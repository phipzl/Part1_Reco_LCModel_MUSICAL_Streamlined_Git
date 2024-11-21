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

# Install Minc-Tools or include all the paths so that the script knows rawtominc, mincmath, dcm2mnc etc.
 . /opt/minc/minc-toolkit-config.sh

# MATLAB Functions Folder
LocalMatDir=$(dirname "${BASH_SOURCE[0]}")
export MatlabFunctionsFolder="$LocalMatDir/MatlabFunctions"

# tmp-folder
tmp_folder="/ceph/mri.meduniwien.ac.at/scratch/radiology/nobackup/tmp_MRSI_processing/Part1"
# tmp_folder=$(pwd)
export tmp_folder

# Gradient Delays (measured at Vienna 7 T scanner, ~2023-10)
export DefaultGradientDelaysForCRTTrajectory="[12.562838, 12.540197, 10.082248]"

# LCModel Path
export LCM_Path="/usr/local/lcmodel/bin/lcmodel"

export RunLCModelOn="lcm"		# Run LCModel on different computer, connecting via ssh. You need a key so that you can automatically connect to this
                            # computer, without needing to type in the password!
                            # BE AWARE THAT THIS COMPUTER HAS TO BE ABLE TO ACCESS THE "LCM_Path", THE BASIS-FILE AND THE "out_path"!
export RunLCModelAs=""  # If you need to be a specific user on the LCModel computer. Leave empty (or dont declare it at all) if not necessary.

# MATLAB Runtime (runs compiled MATLAB code without license)
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/opt/MATLAB_Runtime_R2021b/v911/runtime/glnxa64:/opt/MATLAB_Runtime_R2021b/v911/bin/glnxa64:/opt/MATLAB_Runtime_R2021b/v911/sys/os/glnxa64:/opt/MATLAB_Runtime_R2021b/v911/sys/opengl/lib/glnxa64
export MatlabCompiledFunctions="Matlab_Compiled"

