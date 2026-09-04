#!/bin/bash

# Stop the whole pipeline when a MATLAB step fails, so the exit code reports it.
# Argument $1: name of the step that failed.
matlab_step_failed() {
    echo -e "\n\n$1 failed, stopping.\n\n"
    # run_matlab.sh is also sourced by scripts that have no TerminateProgram.
    if declare -f TerminateProgram >/dev/null; then
        TerminateProgram "$DebugFlag" 1
    fi
    exit 1
}

# An uncaught error aborts the -r statement list before "exit" runs, and MATLAB
# then reads EOF on stdin and quits 0. The try/catch is what makes the status real.
# Argument $1: name of matlab script
run_matlab() {
    if [[ $compiled_matlab_flag -eq 1 ]]; then
        # run the compiled matlab function
        echo -e "\nRun this command: $MatlabCompiledFunctions/$1 $abs_tmp_dir"
        if [[ $2 == "1" ]]; then
	        read -p "stop before matlab call"
        fi
        "$MatlabCompiledFunctions/$1" "$abs_tmp_dir" || matlab_step_failed "$1"
    else
        # run the matlab script $1
        echo -e "\nRun this command: $matlabp -nodisplay -r \"addpath(genpath('$MatlabFunctionsFolder')); $1('$abs_tmp_dir')\""
        if [[ $2 == "1" ]]; then
	        read -p "stop before matlab call"
        fi
        $matlabp -nodisplay -r "try; addpath(genpath('$MatlabFunctionsFolder')); $1('$abs_tmp_dir'); catch ME; disp(getReport(ME)); exit(1); end; exit(0)" || matlab_step_failed "$1"
    fi
}

# The Julia version reconstructs the data, but the LCModel files are still written
# by MATLAB. Returns 0 (true) if that writer is available. MATLAB reaches the .m files
# of this project either through addpath(genpath("$MatlabFunctionsFolder")), so they can
# sit in any subfolder, or through the script folder itself, so search both.
julia_lcm_writer_available() {
    if [[ $compiled_matlab_flag -eq 1 ]]; then
        [[ -x "$MatlabCompiledFunctions/julia_write_lcm_files" ]]
        return
    fi
    local ScriptDir
    ScriptDir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
    [[ -f "$ScriptDir/julia_write_lcm_files.m" ]] && return 0
    [[ -n $(find -L "$MatlabFunctionsFolder" -name julia_write_lcm_files.m -print -quit 2>/dev/null) ]]
}

# WALINET is a lipid decontamination method, an alternative to the L1 and L2
# regularization rather than a step of its own: -L "WALINET,<model>". The model is
# optional; without it WALINET uses its configured default.
walinet_is_the_lipid_decon() {
    [[ $LipidDecon_flag -eq 1 ]] || return 1
    local Method=${LipidDecon_MethodAndNoOfLoops%%,*}
    [[ ${Method^^} == "WALINET" ]]
}

walinet_lipid_decon_model() {
    local Rest=${LipidDecon_MethodAndNoOfLoops#*,}
    [[ $Rest == "$LipidDecon_MethodAndNoOfLoops" ]] && Rest=""
    echo "$Rest"
}

# Argument $1: CurAv argument for run_julia_reco.jl
# Returns non-zero if nothing was reconstructed, so the caller can use MATLAB instead.
run_julia_reconstruction() {
    local OnlyInMatlab=()
    [[ $TwoDCaipParallelImaging_flag -eq 1 ]] && OnlyInMatlab+=("-r (2D-Caipirinha parallel imaging)")
    [[ $SliceParallelImaging_flag -eq 1 ]] && OnlyInMatlab+=("-R (slice parallel imaging)")
    [[ $NonCartTraj_flag -eq 1 ]] && OnlyInMatlab+=("-s (trajectory file)")
    [[ $TimeInterpolation_flag -eq 1 ]] && OnlyInMatlab+=("-T (time interpolation)")
    [[ $FirstOrderPhaseCorr_flag -eq 1 ]] && OnlyInMatlab+=("-F (first order phase correction)")
    [[ $FirstOrderPhaseModulation_flag -eq 1 ]] && OnlyInMatlab+=("-k (first order phase modulation)")
    [[ $NuisRem_flag -eq 1 ]] && OnlyInMatlab+=("-n (nuisance removal)")

    if [[ ${#OnlyInMatlab[@]} -gt 0 ]]; then
        echo -e "\nThe Julia reconstruction does not implement:"
        for Option in "${OnlyInMatlab[@]}"; do
            echo "    $Option"
        done
        return 1
    fi
    if ! julia_lcm_writer_available; then
        echo -e "\nThe Julia output cannot be used, julia_write_lcm_files was not found."
        return 1
    fi
    if [[ $AlignFreq_flag -eq 1 ]] && [[ $AlignFreq_selected_method != "Patref" ]]; then
        echo -e "\nThe Julia reconstruction implements -A Patref only. $AlignFreq_selected_method runs"
        echo "    inside the MATLAB reconstruction, so it would be prepared here and then"
        echo "    discarded. Use -A Patref, which needs no second acquisition, or drop -S."
        return 1
    fi
    if [[ $WaterReference_flag -eq 1 ]] && [[ ${WaterReference_MethodAndFile%%,*} != "W1" ]]; then
        echo -e "\nThe Julia reconstruction supports W1 water referencing only. W2 fits the"
        echo "    water separately, which needs its own LCModel files written for it."
        return 1
    fi

    # Which pass this is has to be read before the run: the water pass creates
    # WaterReference.mat, so afterwards the two passes look alike.
    local IsWaterPass=0
    [[ $WaterReference_flag -eq 1 ]] && [[ ! -f "$out_path/WaterReference.mat" ]] && IsWaterPass=1

    local ScriptDir
    ScriptDir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
    echo -e "\nRun this command: JULIA_NUM_THREADS=$julia_n_threads julia $ScriptDir/run_julia_reco.jl $abs_tmp_dir $1 $julia_mmap"
    JULIA_NUM_THREADS="$julia_n_threads" julia "$ScriptDir/run_julia_reco.jl" "$abs_tmp_dir" "$1" "$julia_mmap" || return 1

    # The LCModel files are written once, after the last average
    if [[ -n "$NumberOfCSIFiles" ]] && [[ $1 -lt $NumberOfCSIFiles ]]; then
        return 0
    fi
    # -A Patref, at the same point in the chain as MATLAB's alignment: on the
    # finished reconstruction, before anything downstream reads it.
    if [[ ${b0_correction_flag:-0} -eq 1 ]]; then
        local B0Python
        B0Python=$(command -v python3 || command -v python)
        if [[ -z $B0Python ]]; then
            echo -e "\nNeither python3 nor python was found, cannot run -A Patref."
            return 1
        fi
        echo -e "\nRun this command: $B0Python $ScriptDir/apply_b0_patref.py $abs_tmp_dir"
        if ! "$B0Python" "$ScriptDir/apply_b0_patref.py" "$abs_tmp_dir"; then
            # Deliberately not "return 1". That falls back to the MATLAB
            # reconstruction, which would report success for a run whose
            # correction never happened.
            echo -e "\napply_b0_patref.py failed, stopping."
            declare -f matlab_step_failed >/dev/null && matlab_step_failed apply_b0_patref.py
            exit 1
        fi
    fi

    if walinet_is_the_lipid_decon; then
        local WalinetPython WalinetModel
        WalinetPython=$(command -v python3 || command -v python)
        WalinetModel=$(walinet_lipid_decon_model)
        if [[ -z $WalinetPython ]]; then
            echo -e "\nNeither python3 nor python was found, cannot run the WALINET lipid decontamination."
            return 1
        fi
        echo -e "\nRun this command: $WalinetPython $ScriptDir/walinet_clean_csi.py $abs_tmp_dir $WalinetModel"
        if ! "$WalinetPython" "$ScriptDir/walinet_clean_csi.py" "$abs_tmp_dir" "$WalinetModel"; then
            # Deliberately not "return 1". That falls back to the MATLAB
            # reconstruction, which would fit uncleaned spectra and report success
            # for a run that asked for the removal.
            echo -e "\nwalinet_clean_csi.py failed, stopping."
            declare -f matlab_step_failed >/dev/null && matlab_step_failed walinet_clean_csi.py
            exit 1
        fi
    fi

    # A W1 water pass reconstructs no metabolites, so it has no spectra to write.
    # Same condition MRSI_Reconstruction.m applies before its LCM-file block.
    if [[ $IsWaterPass -eq 1 ]]; then
        echo -e "\nWater reference pass: coil weights stored, no LCModel files to write."
        return 0
    fi
    if [[ $compiled_matlab_flag -eq 1 ]]; then
        echo -e "\nRun this command: $MatlabCompiledFunctions/julia_write_lcm_files $abs_tmp_dir"
        "$MatlabCompiledFunctions/julia_write_lcm_files" "$abs_tmp_dir"
    else
        echo -e "\nRun this command: $matlabp -nodisplay -r \"addpath(genpath('$MatlabFunctionsFolder')); julia_write_lcm_files('$abs_tmp_dir')\""
        $matlabp -nodisplay -r "try; addpath(genpath('$MatlabFunctionsFolder')); julia_write_lcm_files('$abs_tmp_dir'); catch ME; disp(getReport(ME)); exit(1); end; exit(0)"
    fi
}

# Argument $1: CurAv argument for MRSI_Reconstruction.m
# The MATLAB reconstruction aligns the frequency inside itself, so no python step
# sees it happen and the deep fitting would correct the same data again. MATLAB
# saves AlignFreq_ShiftMap.mat in every branch of the alignment, whichever method
# ran, so that file is the evidence rather than the flag: a run that asked for
# alignment and did not get it leaves no record and is corrected downstream.
record_matlab_alignment() {
    [[ ${AlignFreq_flag:-0} -eq 1 ]] || return 0
    [[ -f "$out_path/AlignFreq/AlignFreq_ShiftMap.mat" ]] || return 0
    local RecordPython RecordScriptDir
    RecordPython=$(command -v python3 || command -v python)
    [[ -n $RecordPython ]] || return 0
    RecordScriptDir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
    "$RecordPython" "$RecordScriptDir/processing_record.py" "$abs_tmp_dir" b0_corrected=true
}

run_mrsi_reconstruction() {
    if [[ $julia_reconstruction -eq 1 ]]; then
        if run_julia_reconstruction "$1"; then
            return 0
        fi
        # MRSI_Reconstruction adds an existing CombinedCSI.mat as a previous average,
        # so a Julia pass that did not finish must not be left behind for it.
        rm -f "$out_path/CombinedCSI.mat"
        echo -e "Use the MATLAB reconstruction instead.\n"
    fi
    # WALINET decontaminates the Julia output; the MATLAB reconstruction has no
    # equivalent, so reconstructing without it would drop the requested step.
    if walinet_is_the_lipid_decon; then
        echo -e "\nWALINET lipid decontamination is implemented for the Julia reconstruction (-S) only."
        declare -f matlab_step_failed >/dev/null && matlab_step_failed "WALINET lipid decontamination"
        exit 1
    fi
    if [[ $compiled_matlab_flag -eq 1 ]]; then
        # run the compiled matlab function
        echo -e "\nRun this command: $MatlabCompiledFunctions/MRSI_Reconstruction $abs_tmp_dir $1"
        if [[ $2 == "1" ]]; then
	        read -p "stop before matlab call"
        fi
        "$MatlabCompiledFunctions/MRSI_Reconstruction" "$abs_tmp_dir" "$1" || matlab_step_failed MRSI_Reconstruction
    else
        # run the matlab script $1
        echo -e "\nRun this command: $matlabp -nodisplay -r \"addpath(genpath('$MatlabFunctionsFolder')); MRSI_Reconstruction('$abs_tmp_dir', $1)\""
        if [[ $2 == "1" ]]; then
	        read -p "stop before matlab call"
        fi
        $matlabp -nodisplay -r "try; addpath(genpath('$MatlabFunctionsFolder')); MRSI_Reconstruction('$abs_tmp_dir', $1); catch ME; disp(getReport(ME)); exit(1); end; exit(0)" || matlab_step_failed MRSI_Reconstruction
    fi
    record_matlab_alignment
}
