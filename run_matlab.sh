#!/bin/bash

# Argument $1: name of matlab script
run_matlab() {
    if [[ $compiled_matlab_flag -eq 1 ]]; then
        # run the compiled matlab function
        echo -e "\nRun this command: $MatlabCompiledFunctions/$1 $abs_tmp_dir"
        if [[ $2 == "1" ]]; then
	        read -p "stop before matlab call"
        fi
        "$MatlabCompiledFunctions/$1" "$abs_tmp_dir"
    else
        # run the matlab script $1
        echo -e "\nRun this command: $matlabp -nodisplay -r \"addpath(genpath('$MatlabFunctionsFolder')); $1('$abs_tmp_dir')\""
        if [[ $2 == "1" ]]; then
	        read -p "stop before matlab call"
        fi
        $matlabp -nodisplay -r "addpath(genpath('$MatlabFunctionsFolder')); $1('$abs_tmp_dir'); exit"
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

# WALINET removal followed by LCModel instead of the deepmrsi fitters. "-Q off
# <model>" asks for the removal without the deep fitting, so the cleaned spectra
# reach LCModel through the ordinary path.
walinet_before_lcmodel() {
    [[ $deep_learning_flag -eq 1 ]] && [[ $deep_learning_fitting == "off" ]] \
        && [[ -n $deep_learning_walinet_model ]] && [[ $deep_learning_walinet_model != "off" ]]
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
    if walinet_before_lcmodel; then
        local WalinetPython
        WalinetPython=$(command -v python3 || command -v python)
        if [[ -z $WalinetPython ]]; then
            echo -e "\nNeither python3 nor python was found, cannot run the WALINET removal."
            return 1
        fi
        echo -e "\nRun this command: $WalinetPython $ScriptDir/walinet_clean_csi.py $abs_tmp_dir $deep_learning_walinet_model"
        "$WalinetPython" "$ScriptDir/walinet_clean_csi.py" "$abs_tmp_dir" "$deep_learning_walinet_model" || return 1
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
        $matlabp -nodisplay -r "addpath(genpath('$MatlabFunctionsFolder')); julia_write_lcm_files('$abs_tmp_dir'); exit"
    fi
}

# Argument $1: CurAv argument for MRSI_Reconstruction.m
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
    if [[ $compiled_matlab_flag -eq 1 ]]; then
        # run the compiled matlab function
        echo -e "\nRun this command: $MatlabCompiledFunctions/MRSI_Reconstruction $abs_tmp_dir $1"
        if [[ $2 == "1" ]]; then
	        read -p "stop before matlab call"
        fi
        "$MatlabCompiledFunctions/MRSI_Reconstruction" "$abs_tmp_dir" "$1"
    else
        # run the matlab script $1
        echo -e "\nRun this command: $matlabp -nodisplay -r \"addpath(genpath('$MatlabFunctionsFolder')); MRSI_Reconstruction('$abs_tmp_dir', $1)\""
        if [[ $2 == "1" ]]; then
	        read -p "stop before matlab call"
        fi
        $matlabp -nodisplay -r "addpath(genpath('$MatlabFunctionsFolder')); MRSI_Reconstruction('$abs_tmp_dir', $1); exit"
    fi
}
