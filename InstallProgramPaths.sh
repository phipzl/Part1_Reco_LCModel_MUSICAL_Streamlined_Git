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
#export matlabp='/usr/local/matlab78/bin/matlab'
export matlabp='/usr/local/bin/matlab'

# Brain extraction tool (bet)
export betp='/ceph/nchirurg.meduniwien.ac.at/imaging/Software/fsl/bin/bet'

# Install Minc-Tools or include all the paths so that the script knows rawtominc, mincmath, dcm2mnc etc.
 . /opt/minc/1.9.17/minc-toolkit-config.sh

# MATLAB Functions Folder
LocalMatDir=`pwd`
if [ -d $LocalMatDir/MatlabFunctions ]; then
	export MatlabFunctionsFolder="$LocalMatDir/MatlabFunctions"
else
	export MatlabFunctionsFolder="$LocalMatDir/MatlabFunctions" #changed bstrasser to lhingerl for concept
fi
MatlabStartupCommand="Paths = regexp(path,':','split');rmpathss = ~cellfun('isempty',strfind(Paths,'MatlabFunctions')); if(sum(rmpathss) > 0);"
export MatlabStartupCommand="${MatlabStartupCommand} x = strcat(Paths(rmpathss), {':'});x = [x{:}]; rmpath(x); end; clear Paths rmpathss x; addpath(genpath('${MatlabFunctionsFolder}'))"

# tmp-folder
export tmp_folder="/ceph/nchirurg.meduniwien.ac.at/imaging_scratch/tmp_MRSI_processing/Part1"


# Gradient Delays (measured at Vienna 7 T scanner, ~2023-10) 
export DefaultGradientDelaysForCRTTrajectory="[12.562838, 12.540197, 10.082248]"


# LCModel Path
export LCM_Path="/ceph/nchirurg.meduniwien.ac.at/lab/.lcmodel/bin/lcmodel"

export RunLCModelOn="nc1"		# Run LCModel on different computer, connecting via ssh. You need a key so that you can automatically connect to this
									# computer, without needing to type in the password!
									# BE AWARE THAT THIS COMPUTER HAS TO BE ABLE TO ACCESS THE "LCM_Path", THE BASIS-FILE AND THE "out_path"!
export RunLCModelAs=""				# If you need to be a specific user on the LCModel computer. Leave empty (or dont declare it at all) if not necessary.


export rawtomincp="/opt/minc/1.9.17/bin/rawtominc"
