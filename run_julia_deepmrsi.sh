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
    if [[ $julia_reconstruction -ne 1 ]]; then
        refuse_julia_only_options
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

patref_is_the_alignment() {
    [[ ${AlignFreq_flag:-0} -eq 1 ]] && [[ ${AlignFreq_MethodAndPath%%,*} == "Patref" ]]
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

# Options only the Julia route implements, refused on the MATLAB reconstruction
# rather than silently dropped.
refuse_julia_only_options() {
    if walrus_is_the_lipid_decon; then
        echo -e "\nWALRUS through -L is implemented for the Julia reconstruction (-S) only."
        matlab_step_failed "WALRUS lipid decontamination"
    fi
    if patref_is_the_alignment; then
        echo -e "\n-A Patref is implemented for the Julia reconstruction (-S) only."
        matlab_step_failed "-A Patref"
    fi
    # The MATLAB reconstruction does not keep the reference scan in CombinedCSI.mat,
    # which deepmrsi needs.
    if [[ -n $(deepmrsi_fitting) ]]; then
        echo -e "\n-l $SpectralFitting_Method is implemented for the Julia reconstruction (-S) only."
        matlab_step_failed "-l $SpectralFitting_Method"
    fi
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
    if [[ $WaterReference_flag -eq 1 ]] && [[ ${WaterReference_MethodAndFile%%,*} != "W1" ]]; then
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

    local Python
    Python=$(command -v python3 || command -v python)
    # -A Patref, on the finished reconstruction, before anything downstream reads it.
    # A failure here stops the run instead of returning, which would report the
    # configuration as unsupported.
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

    if [[ -n $(deepmrsi_fitting) ]]; then
        return 0
    fi
    if [[ $compiled_matlab_flag -eq 1 ]]; then
        echo -e "\nRun this command: $MatlabCompiledFunctions/julia_write_lcm_files $abs_tmp_dir"
        "$MatlabCompiledFunctions/julia_write_lcm_files" "$abs_tmp_dir" || matlab_step_failed julia_write_lcm_files
    else
        echo -e "\nRun this command: $matlabp -nodisplay -r \"addpath(genpath('$MatlabFunctionsFolder')); julia_write_lcm_files('$abs_tmp_dir')\""
        $matlabp -nodisplay -r "try; addpath(genpath('$MatlabFunctionsFolder')); julia_write_lcm_files('$abs_tmp_dir'); catch ME; disp(getReport(ME)); exit(1); end; exit(0)" || matlab_step_failed julia_write_lcm_files
    fi
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
# in <out>/deepMRSI.
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
