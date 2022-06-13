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
export matlabp='/bilbo/usr/local/matlab2013a/bin/matlab'

# Brain extraction tool (bet)
export betp='/usr/local/fsl/bin/bet'

# Install Minc-Tools or include all the paths so that the script knows rawtominc, mincmath, dcm2mnc etc.
 . /opt/minc/minc-toolkit-config.sh

# MATLAB Functions Folder
LocalMatDir=`pwd`
if [ -d $LocalMatDir/Matlab_Functions ]; then
	export MatlabFunctionsFolder="$LocalMatDir/Matlab_Functions"
else
	export MatlabFunctionsFolder="$LocalMatDir/Matlab_Functions" #changed bstrasser to lhingerl for concept
fi
MatlabStartupCommand="Paths = regexp(path,':','split');rmpathss = ~cellfun('isempty',strfind(Paths,'Matlab_Functions')); if(sum(rmpathss) > 0);"
export MatlabStartupCommand="${MatlabStartupCommand} x = strcat(Paths(rmpathss), {':'});x = [x{:}]; rmpath(x); end; clear Paths rmpathss x; addpath(genpath('${MatlabFunctionsFolder}'))"


# LCModel Path
export LCM_Path="/usr/local/lcmodel/bin/lcmodel"

export RunLCModelOn="lcm"		# Run LCModel on different computer, connecting via ssh. You need a key so that you can automatically connect to this
									# computer, without needing to type in the password!
									# BE AWARE THAT THIS COMPUTER HAS TO BE ABLE TO ACCESS THE "LCM_Path", THE BASIS-FILE AND THE "out_path"!
export RunLCModelAs=""				# If you need to be a specific user on the LCModel computer. Leave empty (or dont declare it at all) if not necessary.
