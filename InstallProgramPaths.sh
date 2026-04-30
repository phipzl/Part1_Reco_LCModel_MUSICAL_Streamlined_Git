#!/bin/bash
##### In this file aliases and paths can be stored so that the program knows all commands.

##### The following programs are necessary:
##### -- OS: Ubuntu 12.04			(12.04.3 LTS, current used Kernel: GNU/Linux 3.2.0-48-generic x86_64)
##### -- Minc 						(program: 2.0.18, libminc: 2.0.18, netcdf: 3.6.3, HDF5: 1.6.6,)
##### -- MATLAB 					(matlab78R2009a)
##### -- BET 						(of FSL package 4.1, 2008)
##### -- tar 						(Any version should work)
##### -- gzip						(Any version should work)
##### -- gunzip						(Any version should work)
##### -- LCModel 					(Version 6.3.1)

# aliases

# MATLAB
export matlabp='/bilbo/usr/local/matlab2022a/bin/matlab'

# Brain extraction tool (bet)
export betp='/usr/local/fsl/bin/bet'

mincpath="/opt/minc"

# Install Minc-Tools or include all the paths so that the script knows rawtominc, mincmath, dcm2mnc etc.
 . ${mincpath}/minc-toolkit-config.sh

# MATLAB Functions Folder
LocalMatDir=$(dirname "${BASH_SOURCE[0]}")
export MatlabFunctionsFolder="$LocalMatDir/MatlabFunctions"

# tmp-folder
tmp_folder="/ceph/mri.meduniwien.ac.at/scratch/radiology/nobackup/tmp_MRSI_processing/Part1"
# tmp_folder=$(pwd)
export tmp_folder

# Gradient Delays (measured at Vienna 7 T scanner, ~2023-10)
export DefaultGradientDelaysForCRTTrajectory="GradDelayPerTempInt_x = [12.42 12.38 10.14]; GradDelayPerTempInt_y = [10.27 10.75 8.99];"

# Per Circle:
# export DefaultGradientDelaysForCRTTrajectory="GradDelayPerAngInt_x = [11.4 11.72 11.8 11.84 13.88 13.88 14.8 17.625 13.86 15.9 11.4 15.165 11.205 10.96 11.49 12.48 11.49 10.5 12.48 10.25 10.25 10.275 10.25 11.55 10.53 11.52 9.6 9.57 10.56 9.3 9.275 9.375]; GradDelayPerAngInt_y = [9.2 9.52 9.64 9.76 11.72 11.76 13.2 15.9 12.24 14.22 9.75 13.635 9.54 9.36 9.81 10.86 9.81 8.85 10.95 8.625 8.675 8.625 8.675 10.38 9.39 10.41 8.4 8.4 9.42 8.175 8.175 8.175];"

# LCModel Path
export LCM_Path="/usr/local/lcmodel/bin/lcmodel"

export RunLCModelOn="lcm"		# Run LCModel on different computer, connecting via ssh. You need a key so that you can automatically connect to this
                            # computer, without needing to type in the password!
                            # BE AWARE THAT THIS COMPUTER HAS TO BE ABLE TO ACCESS THE "LCM_Path", THE BASIS-FILE AND THE "out_path"!
export RunLCModelAs=""  # If you need to be a specific user on the LCModel computer. Leave empty (or dont declare it at all) if not necessary.

export rawtomincp="${mincpath}/bin/rawtominc"

# MATLAB Runtime (runs compiled MATLAB code without license)
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/opt/MATLAB_Runtime_R2021b/v911/runtime/glnxa64:/opt/MATLAB_Runtime_R2021b/v911/bin/glnxa64:/opt/MATLAB_Runtime_R2021b/v911/sys/os/glnxa64:/opt/MATLAB_Runtime_R2021b/v911/sys/opengl/lib/glnxa64
export MatlabCompiledFunctions="Matlab_Compiled"

