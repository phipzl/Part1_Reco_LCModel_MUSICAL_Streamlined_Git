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

# Argument $1: CurAv argument for MRSI_Reconstruction.m
run_mrsi_reconstruction() {
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
