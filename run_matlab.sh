#!/bin/bash

# Argument $1: name of matlab script
run_matlab() {
    if [[ $compiled_matlab_flag -eq 1 ]]; then
        # run the compiled matlab function
        echo -e "\nRun this command: $MatlabCompiledFunctions/$1 $abs_tmp_dir"
        "$MatlabCompiledFunctions/$1" "$abs_tmp_dir"
    else
        # run the matlab script $1
        echo -e "\nRun this command: $matlabp -nodisplay -batch \"addpath(genpath('$MatlabFunctionsFolder')); $1('$abs_tmp_dir')\""
        $matlabp -nodisplay -batch "addpath(genpath('$MatlabFunctionsFolder')); $1('$abs_tmp_dir')"
    fi
}

# Argument $1: CurAv argument for MRSI_Reconstruction.m
run_mrsi_reconstruction() {
    if [[ $compiled_matlab_flag -eq 1 ]]; then
        # run the compiled matlab function
        echo -e "\nRun this command: $MatlabCompiledFunctions/MRSI_Reconstruction $abs_tmp_dir $1"
        "$MatlabCompiledFunctions/MRSI_Reconstruction" "$abs_tmp_dir" "$1"
    else
        # run the matlab script $1
        echo -e "\nRun this command: $matlabp -nodisplay -batch \"addpath(genpath('$MatlabFunctionsFolder')); MRSI_Reconstruction('$abs_tmp_dir', $1)\""
        $matlabp -nodisplay -batch "addpath(genpath('$MatlabFunctionsFolder')); MRSI_Reconstruction('$abs_tmp_dir', $1)"
    fi
}
