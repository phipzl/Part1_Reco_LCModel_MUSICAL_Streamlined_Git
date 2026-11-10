#!/bin/bash
# The Julia reconstruction (-S), WALRUS through -L "WALRUS,<model>", the field
# correction -A Patref and the deepmrsi fitters (-l PHIVE, SpatialRegu or Both).
# Sourced by run_matlab.sh, which calls into it only when one of these is given.

# Stop the whole pipeline when a step fails, so the exit code reports it.
# Argument $1: name of the step that failed.
matlab_step_failed() {
    # run_matlab.sh is also sourced by create_mask.sh, a separate process that
    # carries on after a failed step, so only the pipeline itself stops.
    if ! declare -f TerminateProgram >/dev/null; then
        echo -e "\n\n$1 failed.\n\n"
        return 1
    fi
    echo -e "\n\n$1 failed, stopping.\n\n"
    TerminateProgram "$DebugFlag" 1
    exit 1
}

# Called once the options are parsed, before anything is written.
check_julia_deepmrsi_options() {
    # -l takes a method, so a bare -l takes the next option as its method.
    if [[ ${SpectralFittingDontUseLCM_flag:-0} -eq 1 ]] && [[ $SpectralFitting_Method == -* ]]; then
        echo "-l needs a method (LCModel, DeepLearning, None, PHIVE, SpatialRegu or Both), got '$SpectralFitting_Method'."
        matlab_step_failed "-l"
    fi
    # On -S an unknown method would reconstruct and then fit nothing.
    if [[ $julia_reconstruction -eq 1 ]] && [[ ${SpectralFittingDontUseLCM_flag:-0} -eq 1 ]]; then
        case "${SpectralFitting_Method,,}" in
            lcmodel | deeplearning | none | phive | spatialregu | both) ;;
            *)
                echo "-l takes LCModel, None, PHIVE, SpatialRegu or Both, got '$SpectralFitting_Method'."
                matlab_step_failed "-l"
                ;;
        esac
    fi
    # The deep fitters write their maps on csi_template's grid, which -z changes.
    if [[ -n $(deepmrsi_fitting) ]] && [[ ${ZeroFillMetMaps_flag:-0} -eq 1 ]]; then
        echo "-z is not supported with -l $SpectralFitting_Method: the maps are written on the unfilled grid."
        matlab_step_failed "-z"
    fi
    if [[ $julia_reconstruction -ne 1 ]]; then
        matlab_reconstruction_options
    fi
}

# The MATLAB reconstruction with the steps that run after it. What it has to do
# differently: keep the uncombined reference scan, which -A Patref and the deep
# fitters measure the field on and which it drops by default because it is
# several times the size of the CSI.
matlab_reconstruction_options() {
    local Option
    if patref_is_the_alignment || [[ -n $(deepmrsi_fitting) ]]; then
        export KeepFullFIDImage_flag=1
        # The reconstruction resizes the CSI for a target grid but not the reference
        # scan it snapshots before, so the two steps that measure the field on that
        # scan would work across grids. "-I nextpow2" names no target grid.
        if [[ ${InterpolateCSIResolution_flag:-0} -eq 1 && ${InterpolateCSIResolution} == *"["* ]]; then
            echo -e "\n-I with a target grid does not combine with -A Patref or the deep fitters:"
            echo "    they measure the field on the reference scan, which the MATLAB reconstruction"
            echo "    leaves on the acquired grid, and the Julia one does not interpolate at all."
            echo "    Drop -I, or use -A Alignment and LCModel."
            matlab_step_failed "-I with -A Patref or the deep fitters"
        fi
    fi
    # Every -A method other than Patref is applied inside the reconstruction, which
    # writes no record of it, so the fitters would measure the field on the reference
    # scan and correct it a second time.
    if [[ -n $(deepmrsi_fitting) ]] && [[ ${AlignFreq_flag:-0} -eq 1 ]] && ! patref_is_the_alignment; then
        local Method=${AlignFreq_MethodAndPath%%,*}
        [[ -n $Method ]] || Method="(empty)"
        echo -e "\n-A $Method does not combine with -l $SpectralFitting_Method on the MATLAB"
        echo "    reconstruction: it aligns the spectra itself and leaves no record, so the"
        echo "    fitter would correct the field again. Use -A Patref, or fit with LCModel."
        matlab_step_failed "-A $Method with -l $SpectralFitting_Method"
    fi
    patref_is_the_alignment || walrus_is_the_lipid_decon || return 0

    # Both steps run after the reconstruction, so nothing inside it may depend on
    # them: Part1's own Walinet and DeepLearning read the spectra there. With LCModel,
    # the reconstruction also writes the LCModel files there: a W2 water pass writes
    # its water spectra and queues them in a batch line the rewrite replaces, and W1
    # is dropped by that first write but carried by the rewrite, so the same flags
    # would scale the fit or not depending on -A and -L.
    local Refuse=()
    [[ ${NuisRem_flag:-0} -eq 1 ]] && Refuse+=("-n (its nuisance removal runs inside the reconstruction)")
    [[ $SpectralFitting_Method == *DeepLearning* ]] && Refuse+=("-l DeepLearning (it fits inside the reconstruction, before these steps)")
    if [[ $SpectralFitting_Method == "LCModel" && ${WaterReference_flag:-0} -eq 1 ]]; then
        case "$(water_reference_method)" in
            W2) Refuse+=("-w W2 (its water spectra are written and queued inside the reconstruction)") ;;
            *)  Refuse+=("-w W1 (the reconstruction's own LCModel files carry no water reference, the rewritten ones would)") ;;
        esac
    fi
    if [[ ${#Refuse[@]} -gt 0 ]]; then
        echo -e "\nOn the MATLAB reconstruction, -A Patref and -L WALRUS do not combine with:"
        for Option in "${Refuse[@]}"; do echo "    $Option"; done
        echo "    Drop the option, or drop -A Patref and -L WALRUS."
        matlab_step_failed "-A Patref or -L WALRUS on the MATLAB reconstruction"
    fi

    # Spectra the reconstruction is about to sum are aligned with each other first,
    # voxel by voxel onto a reference spectrum, which removes the same field offset
    # the Patref step measures on the reference scan. Correcting it twice is worse
    # than not correcting it, so the combination is refused.
    patref_is_the_alignment || return 0
    if [[ $(grep -oiE "\.dat|\.IMA" <<< "${csi_path:-}" | wc -l) -gt 1 ]]; then
        echo -e "\nSeveral acquisitions given with -c do not combine with -A Patref on the MATLAB"
        echo "    reconstruction: it aligns them with each other before it sums them, which"
        echo "    corrects the same field. Use -A Alignment instead."
        matlab_step_failed "-A Patref on several acquisitions"
    fi
}

walrus_is_the_lipid_decon() {
    [[ $LipidDecon_flag -eq 1 ]] || return 1
    local Method=${LipidDecon_MethodAndNoOfLoops%%,*}
    [[ ${Method^^} == "WALRUS" || ${Method^^} == "WALINET" ]]
}

walrus_lipid_decon_model() {
    local Rest=${LipidDecon_MethodAndNoOfLoops#*,}
    [[ $Rest == "$LipidDecon_MethodAndNoOfLoops" ]] && Rest=""
    echo "$Rest"
}

# The deepmrsi fitter for -l, by the names the online FIRE route uses; empty for
# the other methods.
deepmrsi_fitting() {
    case "${SpectralFitting_Method,,}" in
        phive) echo dlfit ;;
        spatialregu) echo gpufit ;;
        both) echo both ;;
    esac
}

# W1, W2 or empty, classified the way GetPar_CreateTempl_MaskPart1.m does: by the
# first of the two names found anywhere in the -w argument, path included.
water_reference_method() {
    [[ ${WaterReference_flag:-0} -eq 1 ]] || return 0
    case "$WaterReference_MethodAndFile" in
        *W1*) echo W1 ;;
        *W2*) echo W2 ;;
    esac
}

patref_is_the_alignment() {
    [[ ${AlignFreq_flag:-0} -eq 1 ]] || return 1
    local Method=${AlignFreq_MethodAndPath%%,*}
    [[ ${Method,,} == "patref" ]]
}

# InitialParameters.m is MATLAB source. The Julia and Python steps read the same
# content as JSON, written from it so that the two files cannot disagree.
write_initial_parameters_json() {
    local Par="${tmp_dir}/InitialParameters.m"
    local ParJson="${tmp_dir}/InitialParameters.json"
    awk '
BEGIN {
    Quote = sprintf("%c", 39)
    NrOfKeys = 0
}

function JsonEscape(Str,   Out, Pos, Char) {
    Out = ""
    for (Pos = 1; Pos <= length(Str); Pos++) {
        Char = substr(Str, Pos, 1)
        if (Char == "\\")
            Out = Out "\\\\"
        else if (Char == "\"")
            Out = Out "\\\""
        else
            Out = Out Char
    }
    return Out
}

# One line can hold several MATLAB statements ("InPlaneCaipPattern = [...]; VD_Radius = 3;"),
# so cut it at every ";" that is outside brackets and outside a quoted string. Splitting on
# every ";" would tear MATLAB row separators like [1 2; 3 4] apart.
function SplitStatements(Line, Parts,   Pos, Char, Depth, InQuote, Cur, NrOfParts) {
    NrOfParts = 0
    Depth = 0
    InQuote = 0
    Cur = ""
    for (Pos = 1; Pos <= length(Line); Pos++) {
        Char = substr(Line, Pos, 1)
        if (Char == Quote)
            InQuote = !InQuote
        else if (!InQuote && (Char == "[" || Char == "(" || Char == "{"))
            Depth = Depth + 1
        else if (!InQuote && Depth > 0 && (Char == "]" || Char == ")" || Char == "}"))
            Depth = Depth - 1
        if (Char == ";" && !InQuote && Depth == 0) {
            NrOfParts = NrOfParts + 1
            Parts[NrOfParts] = Cur
            Cur = ""
        } else {
            Cur = Cur Char
        }
    }
    NrOfParts = NrOfParts + 1
    Parts[NrOfParts] = Cur
    return NrOfParts
}

# Statements that do not assign to a plain variable name, like the element assignment
# "InPlaneCaipPattern([4 7]) = 1", cannot be expressed in JSON and are skipped here.
function HandleStatement(Line,   EqualPos, Name, Value, CellNr, BracePos, IsQuoted) {
    sub(/^[ \t]+/, "", Line)
    sub(/[ \t]+$/, "", Line)
    if (Line == "" || substr(Line, 1, 1) == "%")
        return

    EqualPos = index(Line, "=")
    if (EqualPos == 0)
        return
    Name = substr(Line, 1, EqualPos - 1)
    Value = substr(Line, EqualPos + 1)
    sub(/[ \t]+$/, "", Name)
    sub(/^[ \t]+/, "", Value)
    sub(/[ \t]+$/, "", Value)
    sub(/;$/, "", Value)
    sub(/[ \t]+$/, "", Value)

    # Cell array entry "Name{Nr} = ...", plain assignment otherwise
    CellNr = 0
    BracePos = index(Name, "{")
    if (BracePos > 0) {
        CellNr = substr(Name, BracePos + 1, index(Name, "}") - BracePos - 1) + 0
        Name = substr(Name, 1, BracePos - 1)
        if (CellNr < 1)
            return
    }
    if (Name !~ /^[A-Za-z_][A-Za-z_0-9]*$/)
        return

    if (!(Name in Seen)) {
        Seen[Name] = 1
        NrOfKeys = NrOfKeys + 1
        Keys[NrOfKeys] = Name
    }

    IsQuoted = (length(Value) >= 2 && substr(Value, 1, 1) == Quote && substr(Value, length(Value), 1) == Quote)
    if (IsQuoted)
        Value = substr(Value, 2, length(Value) - 2)

    if (CellNr > 0) {
        IsCell[Name] = 1
        CellValue[Name, CellNr] = Value
        if (CellNr > CellSize[Name])
            CellSize[Name] = CellNr
    } else {
        SimpleValue[Name] = Value
        # Bare numbers become JSON numbers, everything else a JSON string. So
        # MATLAB arrays and expressions survive as the strings they are.
        IsNumber[Name] = (!IsQuoted && Value ~ /^-?(0|[1-9][0-9]*)([.][0-9]+)?([eE][-+]?[0-9]+)?$/)
    }
}

{
    Line = $0
    gsub(/\r/, "", Line)
    NrOfParts = SplitStatements(Line, Parts)
    for (PartNr = 1; PartNr <= NrOfParts; PartNr++)
        HandleStatement(Parts[PartNr])
}

END {
    print "{"
    for (KeyNr = 1; KeyNr <= NrOfKeys; KeyNr++) {
        Name = Keys[KeyNr]
        Comma = (KeyNr < NrOfKeys) ? "," : ""
        if (Name in IsCell) {
            Entries = ""
            for (CellNr = 1; CellNr <= CellSize[Name]; CellNr++)
                Entries = Entries (CellNr > 1 ? ", " : "") "\"" JsonEscape(CellValue[Name, CellNr]) "\""
            print "    \"" JsonEscape(Name) "\": [" Entries "]" Comma
        } else if (IsNumber[Name]) {
            print "    \"" JsonEscape(Name) "\": " SimpleValue[Name] Comma
        } else {
            print "    \"" JsonEscape(Name) "\": \"" JsonEscape(SimpleValue[Name]) "\"" Comma
        }
    }
    print "}"
}
' "$Par" >"$ParJson"
    chmod 644 "$ParJson"
}

# The Julia version reconstructs the data, but the LCModel files are still written
# by MATLAB. Returns 0 (true) if that writer is available.
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

# Argument $1: CurAv argument for run_julia_reco.jl
# Returns non-zero if the configuration is not supported.
run_julia_reconstruction() {
    local OnlyInMatlab=()
    [[ $TwoDCaipParallelImaging_flag -eq 1 ]] && OnlyInMatlab+=("-r (2D-Caipirinha parallel imaging)")
    [[ $SliceParallelImaging_flag -eq 1 ]] && OnlyInMatlab+=("-R (slice parallel imaging)")
    [[ $NonCartTraj_flag -eq 1 ]] && OnlyInMatlab+=("-s (trajectory file)")
    [[ $TimeInterpolation_flag -eq 1 ]] && OnlyInMatlab+=("-T (time interpolation)")
    [[ $FirstOrderPhaseCorr_flag -eq 1 ]] && OnlyInMatlab+=("-F (first order phase correction)")
    [[ $FirstOrderPhaseModulation_flag -eq 1 ]] && OnlyInMatlab+=("-k (first order phase modulation)")
    [[ $NuisRem_flag -eq 1 ]] && OnlyInMatlab+=("-n (nuisance removal, including Walinet; use -L WALRUS)")
    # "-I nextpow2" too: GetPar then sizes the templates and the mask for the header's final
    # matrix, which this reconstruction never produces.
    [[ ${InterpolateCSIResolution_flag:-0} -eq 1 ]] &&
        OnlyInMatlab+=("-I (this reconstruction stays on the acquired grid)")
    [[ $SpectralFitting_Method == *DeepLearning* ]] && OnlyInMatlab+=("-l DeepLearning (use -l PHIVE)")

    if [[ ${#OnlyInMatlab[@]} -gt 0 ]]; then
        echo -e "\nThe Julia reconstruction does not implement:"
        for Option in "${OnlyInMatlab[@]}"; do
            echo "    $Option"
        done
        return 1
    fi
    # The deepmrsi fitters read CombinedCSI.mat itself and need no LCModel files.
    if [[ -z $(deepmrsi_fitting) ]] && ! julia_lcm_writer_available; then
        echo -e "\nThe Julia output cannot be used, julia_write_lcm_files was not found."
        return 1
    fi
    if [[ $AlignFreq_flag -eq 1 ]] && ! patref_is_the_alignment; then
        echo -e "\nThe Julia reconstruction implements -A Patref only. ${AlignFreq_MethodAndPath%%,*} runs"
        echo "    inside the MATLAB reconstruction. Use -A Patref, which needs no second"
        echo "    acquisition, or drop -S."
        return 1
    fi
    if [[ $WaterReference_flag -eq 1 ]] && [[ $(water_reference_method) != "W1" ]]; then
        echo -e "\nThe Julia reconstruction supports W1 water referencing only. W2 fits the"
        echo "    water separately, which needs its own LCModel files written for it."
        return 1
    fi

    write_initial_parameters_json

    # Which pass this is has to be read before the run: the water pass creates
    # WaterReference.mat, so afterwards the two passes look alike.
    local IsWaterPass=0
    [[ $WaterReference_flag -eq 1 ]] && [[ ! -f "$out_path/WaterReference.mat" ]] && IsWaterPass=1

    # -S "threads,mode": the mode says which reconstruction this run reproduces,
    # ice (the online route, the default) or matlab.
    local Threads=${julia_n_threads%%,*} Mode=ice
    [[ -z $Threads ]] && Threads=auto
    [[ $julia_n_threads == *,* ]] && Mode=${julia_n_threads#*,}
    case "$Mode" in
        ice | matlab) ;;
        *)
            echo -e "\n-S: unknown mode '$Mode'. Use ice or matlab."
            return 1
            ;;
    esac

    local ScriptDir
    ScriptDir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
    echo -e "\nRun this command: JULIA_NUM_THREADS=$Threads julia $ScriptDir/run_julia_reco.jl $abs_tmp_dir $1 $julia_mmap $Mode"
    JULIA_NUM_THREADS="$Threads" julia "$ScriptDir/run_julia_reco.jl" "$abs_tmp_dir" "$1" "$julia_mmap" "$Mode" || return 1

    # The LCModel files are written once, after the last average
    if [[ -n "$NumberOfCSIFiles" ]] && [[ $1 -lt $NumberOfCSIFiles ]]; then
        return 0
    fi

    # A W1 water pass reconstructs no metabolites, so it has no CombinedCSI.mat for the
    # steps below and no spectra to write. Same condition MRSI_Reconstruction.m applies
    # before its LCM-file block.
    if [[ $IsWaterPass -eq 1 ]]; then
        echo -e "\nWater reference pass: coil weights stored, no LCModel files to write."
        return 0
    fi

    run_mindgrab_mask
    run_field_and_lipid_steps

    if [[ -n $(deepmrsi_fitting) ]]; then
        return 0
    fi
    write_lcm_files_from_combined
}

# -m mindgrab: mask the reference scan the reconstruction left in CombinedCSI.mat,
# then let create_mask.sh, which waited for it, read the result as a mask file.
run_mindgrab_mask() {
    [[ ${mindgrab_mask_flag:-0} -eq 1 ]] || return 0
    local Python ScriptDir
    Python=$(command -v python3 || command -v python)
    [[ -n $Python ]] || matlab_step_failed "-m mindgrab (no python3 or python found)"
    ScriptDir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
    echo -e "\nRun this command: $Python $ScriptDir/mindgrab_mask.py $abs_tmp_dir"
    "$Python" "$ScriptDir/mindgrab_mask.py" "$abs_tmp_dir" || matlab_step_failed mindgrab_mask.py
    $rawtomincp -float -clobber -like "$abs_tmp_dir/csi_template.mnc" \
        -input "$abs_tmp_dir/mask_mindgrab.raw" "$mask_file" >/dev/null ||
        matlab_step_failed "the MindGrab mask on the CSI grid"
    (cd "$ScriptDir" && ./create_mask.sh)
    [[ ! -f "$abs_tmp_dir/mask_user_failed" ]] ||
        matlab_step_failed "the MindGrab mask, which holds no brain in the excited volume"
}

# The LCModel inputs, from CombinedCSI.mat as it stands after the steps above.
write_lcm_files_from_combined() {
    if [[ $compiled_matlab_flag -eq 1 ]]; then
        echo -e "\nRun this command: $MatlabCompiledFunctions/julia_write_lcm_files $abs_tmp_dir"
        "$MatlabCompiledFunctions/julia_write_lcm_files" "$abs_tmp_dir" || matlab_step_failed julia_write_lcm_files
    else
        echo -e "\nRun this command: $matlabp -nodisplay -r \"addpath(genpath('$MatlabFunctionsFolder')); julia_write_lcm_files('$abs_tmp_dir')\""
        $matlabp -nodisplay -r "try; addpath(genpath('$MatlabFunctionsFolder')); julia_write_lcm_files('$abs_tmp_dir'); catch ME; disp(getReport(ME)); exit(1); end; exit(0)" || matlab_step_failed julia_write_lcm_files
    fi
}

# -A Patref and WALRUS, on the finished reconstruction, before anything downstream
# reads it. Both routes end here: the steps take CombinedCSI.mat, which the MATLAB
# reconstruction writes in the same layout. A failure stops the run instead of
# returning, which would report the configuration as unsupported.
run_field_and_lipid_steps() {
    patref_is_the_alignment || walrus_is_the_lipid_decon || return 0

    local Python ScriptDir
    Python=$(command -v python3 || command -v python)
    ScriptDir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
    write_initial_parameters_json

    if patref_is_the_alignment; then
        [[ -n $Python ]] || matlab_step_failed "-A Patref (no python3 or python found)"
        echo -e "\nRun this command: $Python $ScriptDir/apply_b0_patref.py $abs_tmp_dir"
        "$Python" "$ScriptDir/apply_b0_patref.py" "$abs_tmp_dir" || matlab_step_failed apply_b0_patref.py
    fi
    if walrus_is_the_lipid_decon; then
        [[ -n $Python ]] || matlab_step_failed "WALRUS (no python3 or python found)"
        local WalrusModel
        WalrusModel=$(walrus_lipid_decon_model)
        echo -e "\nRun this command: $Python $ScriptDir/walrus_clean_csi.py $abs_tmp_dir $WalrusModel"
        "$Python" "$ScriptDir/walrus_clean_csi.py" "$abs_tmp_dir" "$WalrusModel" || matlab_step_failed walrus_clean_csi.py
    fi
}

# The same steps after the MATLAB reconstruction, called once the last average is
# reconstructed. It writes no CombinedCSI.mat for a water reference pass, so there
# is nothing to correct then either.
run_matlab_post_reconstruction() {
    [[ $julia_reconstruction -eq 1 ]] && return 0
    patref_is_the_alignment || walrus_is_the_lipid_decon || [[ -n $(deepmrsi_fitting) ]] || return 0
    [[ -f "$out_path/CombinedCSI.mat" ]] ||
        matlab_step_failed "the MATLAB reconstruction, which wrote no $out_path/CombinedCSI.mat"
    if { patref_is_the_alignment || [[ -n $(deepmrsi_fitting) ]]; } &&
       ! combined_csi_has_reference; then
        echo -e "\nThe MATLAB reconstruction did not keep the reference scan in CombinedCSI.mat."
        echo "    -A Patref and the deep fitters measure the field on it. It is kept when"
        echo "    KeepFullFIDImage is set, which needs the binaries of this Part1 version."
        matlab_step_failed "the reference scan in CombinedCSI.mat"
    fi
    patref_is_the_alignment || walrus_is_the_lipid_decon || return 0
    run_field_and_lipid_steps
    # MRSI_Reconstruction writes the LCModel files inside its own run, before these
    # steps, so those describe the spectra as they were. Write them again from the
    # corrected reconstruction, with the same function the Julia route uses. Its batch
    # script is appended to, so the first one goes: two lines in it fit every voxel
    # twice, at the same time, into the same files.
    if [[ $SpectralFitting_Method == "LCModel" ]]; then
        echo -e "\nThe LCModel files were written before -A and -L; writing them again."
        rm -f "$abs_tmp_dir/lcm_process_core_parallel.sh"
        write_lcm_files_from_combined
    fi
}

# Whether the reconstruction kept image_FullFID, the uncombined reference scan. A
# file that cannot be read is reported as that, not as a missing reference scan.
combined_csi_has_reference() {
    local Python Answer
    Python=$(command -v python3 || command -v python)
    [[ -n $Python ]] || matlab_step_failed "the reference scan check (no python3 or python found)"
    Answer=$("$Python" -c "
import sys, h5py
with h5py.File(sys.argv[1], 'r') as Combined:
    print('yes' if 'image_FullFID' in Combined else 'no')
" "$out_path/CombinedCSI.mat" 2>"$abs_tmp_dir/reference_probe.err") ||
        matlab_step_failed "reading $out_path/CombinedCSI.mat: $(tail -3 "$abs_tmp_dir/reference_probe.err")"
    # A warning on the way out of h5py is not an answer, so only the last line counts.
    [[ ${Answer##*$'\n'} == "yes" ]]
}

# Argument $1: CurAv argument
run_julia_reconstruction_or_stop() {
    run_julia_reconstruction "$1" && return 0
    # No fallback to MATLAB, which needs far more memory for Part1, the reason
    # -S exists. MRSI_Reconstruction adds an existing CombinedCSI.mat as a
    # previous average, so an unfinished Julia pass must not be left behind.
    rm -f "$out_path/CombinedCSI.mat"
    echo -e "\nThe Julia reconstruction does not support this configuration. Change the"
    echo "    setting named above, or drop -S to choose the MATLAB reconstruction."
    matlab_step_failed "the Julia reconstruction"
}

# -l PHIVE, SpatialRegu or Both: fit with deepmrsi instead of LCModel, maps as NIfTI
# in <out>/maps/Orig_<Fitter>.
run_deepmrsi_fit() {
    local OutDir="$out_path" Options=() Python ScriptDir
    write_initial_parameters_json
    # Only the fit: -A and -L ran as Part1's own steps before this one.
    Options+=(--fitting "$(deepmrsi_fitting)")
    Python=$(command -v python3 || command -v python)
    [[ -n $Python ]] || matlab_step_failed "the deepmrsi fit (no python3 or python found)"
    ScriptDir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
    echo -e "\nRun this command: $Python $ScriptDir/run_deepmrsi.py $abs_tmp_dir $OutDir ${Options[*]}"
    "$Python" "$ScriptDir/run_deepmrsi.py" "$abs_tmp_dir" "$OutDir" "${Options[@]}" || matlab_step_failed run_deepmrsi.py
}
