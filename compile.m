%% matlab -nodisplay -batch "run('compile.m')"
% The matlab version for compiling the code must match the matlab runtime version (MCR) for running it. MCR can be downloaded from the MathWorks website and be used without a license.

OutputDir = 'Matlab_Compiled';

%% Compile GetPar_CreateTempl_MaskPart1.m and MRSI_Reconstruction.m

% Include all MRSIMatlabFunctions files as too many are required for manual selection. This increases the size of the compiled bundle
fileList = dir(fullfile('MatlabFunctions/MRSIMatlabFunctions', '**', '*.m'));
% Include other MatlabFunctions files
fileList = [fileList; dir(fullfile('MatlabFunctions/ToolboxCopies/IndividualFunctions', '**', '*.m'))];
fileList = [fileList; dir(fullfile('MatlabFunctions/MatlabFunctions_3rdParty/catstruct', '**', '*.m'))];
fileList = [fileList; dir(fullfile('MatlabFunctions/MatlabFunctions_3rdParty/fieldnamesr', '**', '*.m'))];
fileList = [fileList; dir(fullfile('MatlabFunctions/MatlabFunctions_3rdParty/mapVBVD', '**', '*.m'))];
fileList = [fileList; dir(fullfile('MatlabFunctions/MatlabFunctions_3rdParty/stringdist', '**', '*.m'))];

% These files need to be excluded since they cause errors during compilation
excludeFiles = {'compute_inverse_chemshift_vector_1_0.m' 'compute_SNR_matrix_0_1.m' 'compute_SNR_matrix_0_3.m' 'compute_source_and_target_points_0_1.m' 'RatioMapsExcludeHighValues_1_0.m' 'read_and_plot_phasemaps_1_2.m'};
fileList = fileList(~ismember({fileList.name}, excludeFiles));
fileNames = {fileList.name};
additionalFiles = fullfile({fileList.folder}, fileNames);
additionalFiles = [additionalFiles {'Create_MincTemplates.m'} {'create_mask_VOI.m'} {'create_magnitude.m'} {'MatlabFunctions/MatlabFunctions_3rdParty/Lipid_Suppression_Toolbox_CodeOnly/utils/get_LipidBasis.m'}];

buildResults = compiler.build.standaloneApplication('GetPar_CreateTempl_MaskPart1.m', 'AdditionalFiles', additionalFiles, 'OutputDir', OutputDir)
buildResults = compiler.build.standaloneApplication('MRSI_Reconstruction.m', 'AdditionalFiles', additionalFiles, 'OutputDir', OutputDir)

%% Compile julia_write_lcm_files.m
% The MATLAB side of the -S route. It was missing here, which is why the binary
% shipped with the pipeline predated every other one: there was no way to rebuild it.
buildResults = compiler.build.standaloneApplication('julia_write_lcm_files.m', 'AdditionalFiles', additionalFiles, 'OutputDir', OutputDir)

%% Compile ExtractBrain_mask.m
buildResults = compiler.build.standaloneApplication('ExtractBrain_mask.m', 'OutputDir', OutputDir)

%% Compile flip_mask.m
buildResults = compiler.build.standaloneApplication('flip_mask.m', 'OutputDir', OutputDir)
