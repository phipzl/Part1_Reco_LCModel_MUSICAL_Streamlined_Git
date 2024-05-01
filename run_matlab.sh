#!/bin/bash

# Argument $1: name of matlab script
# Optional argument $2: additional argument for the matlab script
run_matlab() {
    if [[ $compiled_matlab_flag -eq 1 ]]; then
        # run the compiled matlab function
        "$MatlabCompiledFunctions/$1" "$abs_tmp_dir" "$2"
    else
        # run the matlab script $1
        echo -e "\nRun this command: $matlabp -nodisplay -batch \"addpath(genpath('$MatlabFunctionsFolder')); $1('$abs_tmp_dir', '$2')\""
        $matlabp -nodisplay -batch "addpath(genpath('$MatlabFunctionsFolder')); $1('$abs_tmp_dir', '$2')"
    fi
}
