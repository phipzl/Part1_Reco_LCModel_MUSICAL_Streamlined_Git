function GetPar_CreateTempl_MaskPart1(tmp_dir)
%% 0. Preparations

eval(fileread([tmp_dir '/InitialParameters.m']));




%% 1. Save the variables in workspace in the right structures


fprintf('\n\nGather Information\n');
VarsInWorkSpace = who;
VarsInWorkSpace = transpose(VarsInWorkSpace);   % Otherwise you can't loop over the cell.

for var_loopy = VarsInWorkSpace
    
    var_loopy_dummy = var_loopy{:};
    var_loopy_value = eval(var_loopy_dummy);
    
    % Find out if var_loopy is a flag
    if(logical(numel(regexpi(var_loopy_dummy,'_flag'))))
        eval(['Par.Flags.' var_loopy_dummy ' = var_loopy_value;']);
    elseif( logical(numel(regexpi(var_loopy_dummy,'path'))) )
        eval(['Par.Paths.' var_loopy_dummy ' = var_loopy_value;']);
    elseif( logical(numel(regexpi(var_loopy_dummy,'RunLCModel'))) )
        eval(['Par.ServerInfo.' var_loopy_dummy ' = var_loopy_value;']);
    else
        eval(['Par.Settings.' var_loopy_dummy ' = var_loopy_value;']);
    end
    
end

clearvars -except tmp_dir Par ErrorFile
ErrorFile = fopen([tmp_dir '/ErrorFile.sh'],'w+');





%% 2. Get the variables of the LCModel_Control file


%
% % 2.1 SNR Control file
% eval(['run ' Par.Paths.compute_SNR_ControlFile]);
%
% VarsInWorkSpace = who;
% Delete_Partring = logical(cellfun(@numel,strfind(VarsInWorkSpace,'Par')));
% VarsInWorkSpace(Delete_Partring) = [];
% clear Delete_Partring
%
% VarsInWorkSpace = transpose(VarsInWorkSpace);   % Otherwise you can't loop over the cell.
%
% for var_loopy = VarsInWorkSpace
%
%     var_loopy_dummy = var_loopy{:};
%     var_loopy_value = eval(var_loopy_dummy);
%     eval(['Par.ComputeSNR.' var_loopy_dummy ' = var_loopy_value;']);
%
% end




% 2.2 LCModel_Control file
if(isfield(Par.Paths,'LCM_ControlPath'))
    eval(fileread(Par.Paths.LCM_ControlPath));
    Par.LCMControl = ControlWrite;
end


clearvars -except tmp_dir Par ErrorFile

% Write LCM Control file for water dataset
if(isfield(Par.Paths,'LCM_Control_Water_path'))
    eval(fileread(Par.Paths.LCM_Control_Water_path));
    Par.LCMControl_Water = ControlWrite;
end

clearvars -except tmp_dir Par ErrorFile

%% Make some corrections

% Hidden feature: Always use CSI Noise for decorr, if set to -g "UseCSINoise"
if(isfield(Par.Paths,'noisedecorrelation_path') && strcmpi(Par.Paths.noisedecorrelation_path,'UseCSINoise'))
    Par.Flags.noisedecorr_UseCSINoise_flag = 1;
    Par.Paths = rmfield(Par.Paths,'noisedecorrelation_path');
else
    Par.Flags.noisedecorr_UseCSINoise_flag = 0;
end

if(Par.Flags.LipidDecon_flag == 1)
    % L1 or L2 ?
    Par.Flags.LipidDecon_L1_flag = ~isempty(regexp(Par.Settings.LipidDecon_MethodAndNoOfLoops,'L1','ONCE'));
    dummy = regexp(Par.Settings.LipidDecon_MethodAndNoOfLoops,',[0-9]*\.?[0-9]*','match');
    if(Par.Flags.LipidDecon_L1_flag)
        if(isempty(dummy))
            Par.Settings.LipidDecon_NoOfLoops = 10;
        else
            Par.Settings.LipidDecon_NoOfLoops = str2double(dummy{1}(2:end));
        end
    else
        if(isempty(dummy))
            Par.Settings.LipidDecon_L2BetaCorrFactor = 1;
        else
            Par.Settings.LipidDecon_L2BetaCorrFactor = str2double(dummy{1}(2:end));
        end
    end
end


if(Par.Flags.GradientDelay_flag == 1)
    
    if(~isempty(regexp(Par.Settings.GradientDelay,'PerAngInt','ONCE')))
        Tmp = regexp(Par.Settings.GradientDelay,'GradDelayPerAngInt_x = \[\d+.+?\]','match'); % match GradDelayPerAngInt_x = [At least one digit then something at least once, but as few characters as neccessary (specified by ?), and then ]
        Tmp2 = regexp(Tmp{1},'\[.*\]','match');
        Par.Settings.GradientDelayPerAngInt_x = str2num(Tmp2{1});
        Tmp = regexp(Par.Settings.GradientDelay,'GradDelayPerAngInt_y = \[\d+.+?\]','match'); % match GradDelayPerAngInt_x = [At least one digit then something at least once, but as few characters as neccessary (specified by ?), and then ]
        Tmp2 = regexp(Tmp{1},'\[.*\]','match');
        Par.Settings.GradientDelayPerAngInt_y = str2num(Tmp2{1});
    elseif(~isempty(regexp(Par.Settings.GradientDelay,'PerTempInt','ONCE')))
        Tmp = regexp(Par.Settings.GradientDelay,'GradDelayPerTempInt_x = \[\d+.+?\]','match'); % match GradDelayPerAngInt_x = [At least one digit then something at least once, but as few characters as neccessary (specified by ?), and then ]
        Tmp2 = regexp(Tmp{1},'\[.*\]','match');
        Par.Settings.GradientDelayPerTempInt_x = str2num(Tmp2{1});
        Tmp = regexp(Par.Settings.GradientDelay,'GradDelayPerTempInt_y = \[\d+.+?\]','match'); % match GradDelayPerAngInt_x = [At least one digit then something at least once, but as few characters as neccessary (specified by ?), and then ]
        Tmp2 = regexp(Tmp{1},'\[.*\]','match');
        Par.Settings.GradientDelayPerTempInt_y = str2num(Tmp2{1});
    end
end




Par.Flags.ESPIRiT_flag = false;
if(Par.Flags.image_normal_flag)
    if(~isempty(regexpi(Par.Paths.image_normal_path{1},'espirit,','match')))
        Par.Flags.ESPIRiT_flag = true;
        Par.Paths.image_normal_path{1} = regexprep(Par.Paths.image_normal_path{1},'espirit,','','ignorecase');
    end
    % If we want to do ESPIRiT with MUSICAL
    if(isempty(Par.Paths.image_normal_path{1}))
        Par.Flags.image_normal_flag = false;
        Par.Paths = rmfield(Par.Paths,'image_normal_path');
    end
end

Par.Flags.zerofill_to_nextpow2_flag = false;
if(Par.Flags.InterpolateCSIResolution_flag == 1)
    
    Par.Settings.InterpolateCSIResolution_InkSpace = false;
    if(~isempty(regexpi(Par.Settings.InterpolateCSIResolution,'kspace','match')))
        Par.Settings.InterpolateCSIResolution_InkSpace = true;
        Par.Settings.InterpolateCSIResolution_EllipFilter = false;
        if(~isempty(regexpi(Par.Settings.InterpolateCSIResolution,'Ellip','match')))
            Par.Settings.InterpolateCSIResolution_EllipFilter = true;
        end
    end
    
    
    Par.Flags.zerofill_to_nextpow2_flag = ~isempty(regexp(Par.Settings.InterpolateCSIResolution,'nextpow','ONCE'));
    dummy = regexp(Par.Settings.InterpolateCSIResolution,'\[[\d\s,;.]*\]','match');
    try
        Par.Settings.InterpolateCSIResolution = eval(dummy{:});
        Par.Settings.InterpolateCSIResolution = [Par.Settings.InterpolateCSIResolution zeros([1 3-numel(Par.Settings.InterpolateCSIResolution)]) ];	% Guarantee for 3 dims
    catch errie
        Par.Settings.InterpolateCSIResolution = [0 0 0];	% zeros will be replaced later by actual size
    end
else
    Par.Settings.InterpolateCSIResolution = [0 0 0];
end
%% Water Reference Processing
fprintf('\n\nDetermine if and which Water Reference Processing should be done...\n');

if(Par.Flags.WaterReference_flag == 1)
    dummy = char(strtrim(regexp(Par.Settings.WaterReference_MethodAndFile,'[^,]*$','match')));
    % W1 or W2 ?
    if(~isempty(regexp(Par.Settings.WaterReference_MethodAndFile,'W1','ONCE')))
        Par.Settings.WaterReference_method = 'W1';
        Par.Paths.WaterReference_path = dummy;
    elseif(~isempty(regexp(Par.Settings.WaterReference_MethodAndFile,'W2','ONCE')))
        Par.Settings.WaterReference_method = 'W2';
        Par.Paths.WaterReference_path = dummy;
        if(~isfield(Par.Paths,'LCM_Control_Water_path'))
            warning('Missing path to Water Referencing Control File, setting Water flag to zero');
            Par.Flags.WaterReference_flag = 0;
            Par.Settings.WaterReference_method = '0';
        end
    else
        warning(['Not a correct format for water referencing, setting Water flag to zero. ' ...
                 'Expected ''W1,<path>'' or ''W2,<path>'', got ''%s''.'], ...
                 Par.Settings.WaterReference_MethodAndFile);
        Par.Flags.WaterReference_flag = 0;
        Par.Settings.WaterReference_method = '0';
    end
else
    Par.Settings.WaterReference_method = '0';  % This is needed, because many following functions have WaterRefMethod in conditions
end

% Error handling
if(isfield(Par.Paths,'LCM_Control_Water_path') && strcmp(Par.Settings.WaterReference_method, 'W1'))
    Par.Paths=rmfield(Par.Paths,'LCM_Control_Water_path');
    warning('W1 does not use own Control file - Water Control File was removed');
end

%% Create water_spectra folder if W2 is set
if strcmp(Par.Settings.WaterReference_method, 'W2')
    mkdir(sprintf('%s/water_spectra', Par.Paths.out_path));
end


%% 3. Test if all files passed over exist

fields_files = fields(Par.Paths);
[bla,CurComp] = unix('hostname');
TestOnServer = {''};
ServerSSH = '';


if(isfield(Par,'ServerInfo') && ~strcmp(Par.ServerInfo,CurComp))
    TestOnServer = {'LCM_Path','basis_path','out_path'};
    if(strcmpi(Par.Paths.LCM_Path,'lcmodel'))
        TestOnServer(1) = [];
    end
    ServerSSH = [Par.ServerInfo.RunLCModelAs '@' Par.ServerInfo.RunLCModelOn]; ServerSSH(ServerSSH(1) == '@') = []; % Delete @ symbol if it is on first place
end
if(strcmpi(Par.Paths.LCM_Path,'lcmodel'))
    fields_files = setdiff(fields_files,'LCM_Path');
end
TestLocal = setdiff(fields_files,TestOnServer(1:end-1));

% Test for existence
fprintf('\n\nTest if folders exist locally and on server via %s.\nIF PROGRAM STOPS HERE, THE ssh CONNECTION IS NOT WORKING WITH AUTO-LOGIN!\n\n',ServerSSH)
if(~TestForPathExistence(Par.Paths,0,TestLocal,ServerSSH,TestOnServer))
    fprintf(ErrorFile,'ErrorInGetPar_CreateTempl=1; ErrorMessage=''Input path not existing.''');
    fclose(ErrorFile);
    return;
end


[bla, out_dir_attr] = fileattrib(Par.Paths.out_path);
if(~exist(Par.Paths.out_path,'dir') || out_dir_attr.UserWrite == 0 )
    fprintf('\nout_dir\n%s\ndoes not exist or no write rights to it. Aborting . . .\n',Par.Paths.out_path)
    fprintf(ErrorFile,'ErrorInGetPar_CreateTempl=1; ErrorMessage=''Output path not existing.''');
    fclose(ErrorFile);
    return;
end

clear TestLocal bla out_dir_attr;



%% If Folder with DICOMs is passed over, make list of DICOM files

BakCSIPath = Par.Paths.csi_path;
Par.Flags.DICOM_flag = false;
if(numel(Par.Paths.csi_path) == 1)
    Existence = exist(Par.Paths.csi_path{1},'file');
    if(Existence == 7) % Return value is 7 for folder
        

        csi_path_allfiles = dir( fullfile(Par.Paths.csi_path{1},'*.*') );
        csi_path_allfiles = natsortfiles(csi_path_allfiles);
        csi_path_allfiles = {csi_path_allfiles.name}';
        DicomFoundVec = ~cellfun(@isempty,regexpi(csi_path_allfiles,'\.IMA|\.dcm'));
        csi_path_allfiles = csi_path_allfiles(DicomFoundVec);

%         csi_path_allfiles = dir( fullfile(Par.Paths.csi_path{1},'*.IMA') );
        csi_path_allfiles = strcat(Par.Paths.csi_path{1},'/',csi_path_allfiles);
        if(numel(csi_path_allfiles) > 0)
            Par.Paths.csi_path = csi_path_allfiles;
            Par.Flags.DICOM_flag = true;
        else
            fprintf('\nError in GetPar_CreateTempl_MaskPart1.m: Could not find any DICOM files in the folder\n%s. Abort.\n',Par.Paths.csi_path);
            return;
        end
        
    end
end




%% 3. Get the CSI Information

if(isempty(regexp(Par.Paths.csi_path{1},'.*/\w*\.mat','ONCE')))			% Only read header of first csi file. The others must be the same
    
    if(strcmp(Par.Settings.WaterReference_method, 'W2') && ~exist([Par.Paths.out_path '/Parameters_water.mat' ],'file'))
        Par.Paths.csi_path{1} = Par.Paths.WaterReference_path;
        fprintf('\n\nParameters_water.mat are being created.\n\n')
    else
        fprintf('\n\nParameters.mat are being created.\n\n')
    end
    csiPars = read_ascconv(Par.Paths.csi_path{1});		% New function which is adapted to VE11. If some bug occurs, search here.
    % Total Channel Number
    csiPars.total_channel_no = csiPars.total_channel_no_reco; csiPars = rmfield(csiPars,{'total_channel_no_measured','total_channel_no_reco'});
    
    
    % DICOM is always zerofilled to next power of 2.
    if(Par.Flags.zerofill_to_nextpow2_flag || numel(strfind(Par.Paths.csi_path{1}, '.IMA')) > 0)
        csiPars.nFreqEnc = csiPars.nFreqEnc_FinalMatrix;
        csiPars.nPhasEnc = csiPars.nPhasEnc_FinalMatrix;
        if(csiPars.ThreeD_flag)
            csiPars.nPartEnc = csiPars.nSLC_FinalMatrix;    % Only in ThreeD Case we have to change that.
        end
    end
    
    
    % Define Matrix Size
    dummy = [csiPars.nFreqEnc csiPars.nPhasEnc csiPars.nPartEnc*csiPars.nSLC];
    
    Par.Settings.InterpolateCSIResolution(Par.Settings.InterpolateCSIResolution == [0 0 0]) = dummy(Par.Settings.InterpolateCSIResolution == [0 0 0]);
    Par.Settings.InterpolateResolutionRatio = [csiPars.nFreqEnc csiPars.nPhasEnc csiPars.nPartEnc*csiPars.nSLC] ./ Par.Settings.InterpolateCSIResolution;
    
    if(sum(Par.Settings.InterpolateResolutionRatio > 1) >= 1 && ~Par.Settings.InterpolateCSIResolution_InkSpace && ~Par.Flags.AlignFreq_flag)
        fprintf('\n\nWarning: Some spectra will be summed, but no frequency alignment was requested. I set AlignFreq_flag = true.\n\n')
        Par.Flags.AlignFreq_flag = true;
    end
    if(numel(BakCSIPath) > 1 && ~Par.Flags.AlignFreq_flag)  % Need BakCSIPath here, bc we temporarily overwrite this with all dicom files if we read in dicoms
        fprintf('\n\nWarning: Perform averaging of several csi data sets, but no frequency alignment was requested. I set AlignFreq_flag = true.\n\n')
        Par.Flags.AlignFreq_flag = true;
    end
    
    
    
    
    csiPars.nFreqEnc = csiPars.nFreqEnc / Par.Settings.InterpolateResolutionRatio(1);
    csiPars.nPhasEnc = csiPars.nPhasEnc / Par.Settings.InterpolateResolutionRatio(2);
    csiPars.nPartEnc = csiPars.nPartEnc / Par.Settings.InterpolateResolutionRatio(3);
    csiPars.nFreqEnc_nonzf = csiPars.nFreqEnc * Par.Settings.InterpolateResolutionRatio(1);
    csiPars.nPhasEnc_nonzf = csiPars.nPhasEnc * Par.Settings.InterpolateResolutionRatio(2);
    csiPars.nPartEnc_nonzf = csiPars.nPartEnc * Par.Settings.InterpolateResolutionRatio(3);
    csiPars.nFreqEnc_FinalMatrix = csiPars.nFreqEnc_FinalMatrix / Par.Settings.InterpolateResolutionRatio(1);
    csiPars.nPhasEnc_FinalMatrix = csiPars.nPhasEnc_FinalMatrix / Par.Settings.InterpolateResolutionRatio(2);
    csiPars.nSLC_FinalMatrix = csiPars.nSLC_FinalMatrix / Par.Settings.InterpolateResolutionRatio(3);
    
    
    
    
    %csiPars.nFreqEnc = csiPars.nFreqEnc *2;        uncomment if zerofill wanted
    %csiPars.nPhasEnc = csiPars.nPhasEnc *2;
    csiPars = rmfield(csiPars,{'nFreqEnc_FinalMatrix','nPhasEnc_FinalMatrix','nSLC_FinalMatrix' });
    if(~csiPars.ThreeD_flag && numel(strfind(Par.Paths.csi_path{1}, '.IMA')) > 0)
        csiPars.nSLC = 1;
    end
    
    % Stepsize (Voxel Size)
    if(csiPars.nFreqEnc == 1)
        csiPars.FoV_Read = csiPars.VoI_Read;
    end
    if(csiPars.nPhasEnc == 1)
        csiPars.FoV_Phase = csiPars.VoI_Phase;
    end
    csiPars.StepRead = -csiPars.FoV_Read(1) / csiPars.nFreqEnc;		% Coordinate system is reversed in minc with respect to DICOM
    csiPars.StepPhase = -csiPars.FoV_Phase(1) / csiPars.nPhasEnc;
    if(csiPars.nPartEnc == 1)
        csiPars.FoV_Partition(1) = csiPars.VoI_Partition(1);
    end
    csiPars.StepSlice = csiPars.FoV_Partition(1) / (csiPars.nPartEnc * csiPars.nSLC);
    
    
    
    if(~isfield(Par,'Settings') || ~isfield(Par.Settings,'ZeroFillMetMaps'))
        Par.Settings.ZeroFillMetMaps = 1;
    end
    zff = Par.Settings.ZeroFillMetMaps;
    
    
    csiPars.StepRead_zf = csiPars.StepRead / zff;
    csiPars.StepPhase_zf = csiPars.StepPhase / zff;
    if(csiPars.ThreeD_flag)
        csiPars.StepSlice_zf = csiPars.StepSlice / zff; % TODO  this should be implemented later, when zerofilling is done in z direction - csiPars.StepSlice_zf = csiPars.StepSlice / 2;
    else
        csiPars.StepSlice_zf = csiPars.StepSlice;
    end
    
    
    % Rename Position Fields
    [csiPars.POS_X] = csiPars.Pos_Sag;
    [csiPars.POS_Y] = csiPars.Pos_Cor;
    [csiPars.POS_Z] = csiPars.Pos_Tra;
    csiPars = rmfield(csiPars,{'Pos_Sag','Pos_Cor','Pos_Tra'});
    
    % Only for CRT?
    csiPars.InPlaneRotation=csiPars.InPlaneRotation;
    csiPars.InPlaneRotation_VOI=csiPars.InPlaneRotation_VOI;
    
    % Compute direction cosine from x,y and z components of slice normal vector
    [csiPars.PhaseNormalVector, csiPars.ReadNormalVector] = compute_dircos([csiPars.SliceNormalVector_x(1) csiPars.SliceNormalVector_y(1) csiPars.SliceNormalVector_z(1)],csiPars.InPlaneRotation);
    csiPars.SliceNormalVector = [csiPars.SliceNormalVector_x(1) csiPars.SliceNormalVector_y(1) csiPars.SliceNormalVector_z(1)];
    
    MinusVec1 = [-1 -1 1];        % MinusVec are here only for "tuning" signs of the final rotation matrix (RotMat). The values of MinusVecs are
    MinusVec2 = [1 1 -1];         % empirical, thus there might exist a dataset, which will have signs of RotMat uncorrect. However, this setup
    MinusVec3 = [-1 -1 +1];       % works for all tested datasets. The uncorrect signs of direction cosines might be caused by the fact, that the attribute
    % Image Orientation (Patient) tag(0020,0037) is not used in computation of dircos. Instead the Siemens private Tag (0029, 1020)
    % is used
    
    
    csiPars.ReadNormalVector = csiPars.ReadNormalVector .* MinusVec1;
    csiPars.PhaseNormalVector = csiPars.PhaseNormalVector .* MinusVec2;
    csiPars.SliceNormalVector = csiPars.SliceNormalVector .* MinusVec3;
    
    % Create rotation matrix
    RotMat = cat(1,csiPars.ReadNormalVector,csiPars.PhaseNormalVector,csiPars.SliceNormalVector);
    
    % Reverse x- and y- coordinates due to the reversed coordinate system of minc with respect to DICOM
    Pos = [-csiPars.POS_X(1),-csiPars.POS_Y(1),csiPars.POS_Z(1)];
    
    % Convert position from DICOM world coordinates to MINC start values
    % This approach is probably prone to extreme rotation of FOV and is not
    % universal, however for the tested datasets it yielded correct results
    Pos_Minc = RotMat * transpose(Pos);
    % The following line is old, and I think wrong. What I think happened is: Michal fixed the shift
    % in the z-direction by half a voxel with the line below (effectively this line calculates
    % Pos_z - FoV_z/2 + FoV_z/(2*N_z). Then I figured out that the x- and y-positions have to be
    % shifted by half a voxel, and thought by analogy also the z-dimension has to be shifted, not
    % knowing that Michal did that already with the FoVHalf. Now it should be fixed:
    % The FoVHalf is defined "normally" also for z, and the half-voxel shift is done in the 3D-case.
    % For checks: See git commits #1383, #004f, #a9be
    % 	FoVHalf = [csiPars.FoV_Read(1)/2 csiPars.FoV_Phase(1)/2 -csiPars.FoV_Partition(1)/csiPars.nPartEnc*(csiPars.nPartEnc-1)/2];
    FoVHalf = [csiPars.FoV_Read(1)/2 csiPars.FoV_Phase(1)/2 -csiPars.FoV_Partition(1)/2];
    Pos_Minc = transpose(Pos_Minc) + FoVHalf;
    
    % Get from Center of Voxel (DICOM) to corner of voxel (minc) by subtracting half the voxel
    csiPars.POS_X_FirstVoxel = Pos_Minc(1) + csiPars.StepRead/2;         % Be aware that StepRead and StepPhase are reversed and thus the sum is effectively a subtraction.
    csiPars.POS_Y_FirstVoxel = Pos_Minc(2) + csiPars.StepPhase/2;
    csiPars.POS_Z_FirstVoxel = Pos_Minc(3);
    
    
    if(csiPars.ThreeD_flag)
        csiPars.POS_Z_FirstVoxel = csiPars.POS_Z_FirstVoxel + csiPars.StepSlice/2;
    end
    
    
    
    if(Par.Flags.InterpolateCSIResolution_flag == 1)
        % Original CSI Sizes w/o Interpolation
        csiPars.nFreqEnc_BefInterpol = csiPars.nFreqEnc * Par.Settings.InterpolateResolutionRatio(1);
        csiPars.nPhasEnc_BefInterpol = csiPars.nPhasEnc * Par.Settings.InterpolateResolutionRatio(2);
        csiPars.nPartEnc_BefInterpol = csiPars.nPartEnc * Par.Settings.InterpolateResolutionRatio(3);
        csiPars.nSLC_BefInterpol = csiPars.nSLC;
        
        csiPars.StepRead_BefInterpol = -csiPars.FoV_Read(1) / csiPars.nFreqEnc_BefInterpol;
        csiPars.StepPhase_BefInterpol = -csiPars.FoV_Phase(1) / csiPars.nPhasEnc_BefInterpol;
        csiPars.StepSlice_BefInterpol = csiPars.FoV_Partition(1) / (csiPars.nPartEnc_BefInterpol * csiPars.nSLC_BefInterpol);
        
        FoVHalf_BefInterpol = FoVHalf;
        Pos_Minc_BefInterpol = Pos_Minc;
        
        csiPars.POS_X_FirstVoxel_BefInterpol = Pos_Minc_BefInterpol(1) + csiPars.StepRead_BefInterpol/2;
        csiPars.POS_Y_FirstVoxel_BefInterpol = Pos_Minc_BefInterpol(2) + csiPars.StepPhase_BefInterpol/2;
        csiPars.POS_Z_FirstVoxel_BefInterpol = Pos_Minc_BefInterpol(3);
        
        if(csiPars.ThreeD_flag)
            csiPars.POS_Z_FirstVoxel_BefInterpol = csiPars.POS_Z_FirstVoxel_BefInterpol + csiPars.StepSlice_BefInterpol/2;
        end
        
    end
    
    if(Par.Flags.basis_echo_flag == 1) %for spin echo with smaller matrix
        % Original CSI Sizes w/o Interpolation
        csiPars.nFreqEnc_ECHO = csiPars.nFreqEnc /2;
        csiPars.nPhasEnc_ECHO = csiPars.nPhasEnc /2;
        csiPars.nPartEnc_ECHO = csiPars.nPartEnc;
        csiPars.nSLC_ECHO = csiPars.nSLC;
        
        csiPars.StepRead_ECHO = -csiPars.FoV_Read(1) / csiPars.nFreqEnc_ECHO ;
        csiPars.StepPhase_ECHO = -csiPars.FoV_Phase(1) / csiPars.nPhasEnc_ECHO ;
        csiPars.StepSlice_ECHO = csiPars.FoV_Partition(1) / (csiPars.nPartEnc_ECHO * csiPars.nSLC_ECHO);
        
        FoVHalf_ECHO = FoVHalf; FoVHalf_ECHO(3) = -csiPars.FoV_Partition(1)/csiPars.nPartEnc_ECHO*(csiPars.nPartEnc_ECHO-1)/2;
        Pos_Minc_ECHO = Pos_Minc - FoVHalf + FoVHalf_ECHO;
        
        csiPars.POS_X_FirstVoxel_ECHO = Pos_Minc_ECHO(1) + csiPars.StepRead_ECHO/2;
        csiPars.POS_Y_FirstVoxel_ECHO = Pos_Minc_ECHO(2) + csiPars.StepPhase_ECHO/2;
        csiPars.POS_Z_FirstVoxel_ECHO = Pos_Minc_ECHO(3);
        
        %if(csiPars.ThreeD_flag)
        %	csiPars.POS_Z_FirstVoxel_ECHO = csiPars.POS_Z_FirstVoxel_ECHO + csiPars.StepSlice_ECHO/2;
        %end
        
    end
    
    
    % Read Patientname
    if(isfield(Par.Paths,'T1w_path'))
        csiPars.PatName = io_ReadPatientName(Par.Paths.T1w_path);
    else
        csiPars.PatName = io_ReadPatientName(Par.Paths.csi_path{1});
    end
    
    Par.CSI = csiPars; clear csiPars
    
    
    % At this point, a bunch of interesting data is loaded into MATLAB, including B0 [T], w0 [Hz], transmitter voltage, age, sex, weight, flip angle, which we save in the tmp_dir

    PatientData = ReadMeasurementInfoFromMRSI(Par);


    
    fid=fopen([Par.Settings.tmp_dir '/MeasurementInfos.txt'], 'w');
    fprintf(fid, 'Name DoB Age Sex Weight w_0 B_0 U_transmit alpha_E\n%s\n',PatientData);
    fclose(fid);
    
    
    % Check data size and decide if there can by any means be Parallel Imaging enabled
    file_fid = fopen(sprintf('%s', Par.Paths.csi_path{1}),'r');
    headersize = fread(file_fid,1, 'uint32');
    fclose(file_fid);
    size_csidata = dir(Par.Paths.csi_path{1});
    size_csidata = size_csidata.bytes - headersize;
    
    if(isfield(Par.CSI,'WipMemBlockInterpretation') && isfield(Par.CSI.WipMemBlockInterpretation,'Prescan'))
        
        if(isfield(Par.CSI.WipMemBlockInterpretation.Prescan, 'PATREFANDIMASCAN') )
            size_csidata = size_csidata - Par.CSI.WipMemBlockInterpretation.Prescan.PATREFANDIMASCAN.nPhasEnc*2 ...			% Frequency Encoding (sloppily coded...)
                * Par.CSI.WipMemBlockInterpretation.Prescan.PATREFANDIMASCAN.nPhasEnc * Par.CSI.WipMemBlockInterpretation.Prescan.PATREFANDIMASCAN.nSLC ...
                * Par.CSI.WipMemBlockInterpretation.Prescan.PATREFANDIMASCAN.nAverages * 4*2;	% 4: float32, 2: real, imag
        end
        
        if(isfield(Par.CSI.WipMemBlockInterpretation.Prescan, 'NOISEADJSCAN') )
            size_csidata = size_csidata - Par.CSI.WipMemBlockInterpretation.Prescan.NOISEADJSCAN.nReadEnc ...
                * 4000 * 4*2;	% 4000: Sloppily coded size of each vector... 4: float32, 2: real, imag
        end
        
    end
    
    if(Par.CSI.Full_ElliptWeighted_Or_Weighted_Acq == 2)
        EllipFilt = EllipticalFilter(ones([Par.CSI.nFreqEnc Par.CSI.nPhasEnc]),[1 2],[1 1 1 floor(Par.CSI.nFreqEnc/2)-1],1);
        nSpatEncCSI = sum(sum(EllipFilt)); clear EllipFilt;
    else
        nSpatEncCSI = Par.CSI.nFreqEnc*Par.CSI.nPhasEnc;
    end
    size_csidata_predicted = Par.CSI.total_channel_no * nSpatEncCSI * Par.CSI.nPartEnc * Par.CSI.nSLC * Par.CSI.vecSize * 4 * 2;  % 4: float32, 2: real imag
    
    size_csidata_ratio = size_csidata_predicted / size_csidata;
    
    
    
    
    
    
    % Set the Parallel Imaging Flags and values, if info is written in the header of the CSI data
    % 1D Caipi. Only if the size ratio is above 1.7
    if(size_csidata_ratio > 1.7)
        % Check if info is available in wipmemblock
        OneDCaipInfoAvail = isfield(Par.CSI,'WipMemBlockInterpretation') && isfield(Par.CSI.WipMemBlockInterpretation,'OneDCaipi') && isstruct(Par.CSI.WipMemBlockInterpretation.OneDCaipi) ...
            && Par.CSI.WipMemBlockInterpretation.OneDCaipi.SliceParallelImaging_flag;
        % Check if the user gave info. Dont touch user given info.
        if(OneDCaipInfoAvail && ~isfield(Par.Settings,'SliceAliasingPattern') && ~isfield(Par.Settings,'FoVShifts_x') && ~isfield(Par.Settings,'FoVShifts_y') )
            Par.Settings.SliceAliasingPattern = Par.CSI.WipMemBlockInterpretation.OneDCaipi.SliceAliasingPattern;
            Par.Settings.FoVShifts_x = Par.CSI.WipMemBlockInterpretation.OneDCaipi.FoVShifts_x;
            Par.Settings.FoVShifts_y = Par.CSI.WipMemBlockInterpretation.OneDCaipi.FoVShifts_y;
            Par.Flags.SliceParallelImaging_flag = Par.CSI.WipMemBlockInterpretation.OneDCaipi.SliceParallelImaging_flag;
        end
    end
    % Make SliceAliasing Patterns like [1 2], [1 3; 2 4], [1 4; 2 5;3 6] ...
    if(Par.Flags.SliceParallelImaging_flag && ~isfield('Par.Settings','SliceAliasingPattern'))
        for AliLoop = 1:numel(Par.Settings.FoVShifts_x)/2
            Par.Settings.SliceAliasingPattern(AliLoop,:) = [AliLoop AliLoop+numel(Par.Settings.FoVShifts_x)/2];
        end
    end
    
    % Calculate ratio again with new info
    if(Par.Flags.SliceParallelImaging_flag)
        size_csidata_predicted = size_csidata_predicted / size(Par.Settings.SliceAliasingPattern,2);
    end
    size_csidata_ratio = size_csidata_predicted / size_csidata;
    
    
    % 2D Caipi. Only if the size ratio is above 1.4
    if(size_csidata_ratio > 1.4)
        TwoDCaipInfoAvail = isfield(Par.CSI,'WipMemBlockInterpretation') && isfield(Par.CSI.WipMemBlockInterpretation,'TwoDCaipi') && isstruct(Par.CSI.WipMemBlockInterpretation.TwoDCaipi);
        
        TwoDCaipInfoDefective = TwoDCaipInfoAvail && ~isstruct(Par.CSI.WipMemBlockInterpretation.TwoDCaipi) && Par.CSI.WipMemBlockInterpretation.TwoDCaipi == -1;
        if(TwoDCaipInfoDefective && ~isfield(Par.Settings,'InPlaneCaipPattern'))
            fprintf('\nERROR: 2D-CAIPI seems to have been performed, but I could not read the Pattern.')
            fprintf('\nInput the 2D-Pattern manually with option r\n(e.g.: -r ''InPlaneCaipPattern = zeros([5 5]); InPlaneCaipPattern([4 7 9 15 18 21]) = 1; VD_Radius = 1;'').\nStopping here.')
            
            fprintf(ErrorFile,'ErrorInGetPar_CreateTempl=1; ErrorMessage=''TwoDCaipiInfo defective in ascconv header.''');
            fclose(ErrorFile);
            return;
        end
        
        if(TwoDCaipInfoAvail && ~isfield(Par.Settings,'InPlaneCaipPattern') && ~isfield(Par.Settings,'VD_Radius'))
            Par.Settings.InPlaneCaipPattern = Par.CSI.WipMemBlockInterpretation.TwoDCaipi.Skip_Matrix;
            Par.Settings.VD_Radius = Par.CSI.WipMemBlockInterpretation.TwoDCaipi.VD_Radius;
            Par.Flags.TwoDCaipParallelImaging_flag = Par.CSI.WipMemBlockInterpretation.TwoDCaipi.TwoDCaipParallelImaging_flag;
        end
    end
    
    
    
    
    
    if(isfield(Par.Settings,'InPlaneCaipPattern') && sum(sum(Par.Settings.InPlaneCaipPattern)) == numel(Par.Settings.InPlaneCaipPattern) ...
            || isfield(Par.Settings,'VD_Radius') && Par.Settings.VD_Radius > floor(Par.CSI.nFreqEnc/2)-1)		% i.e. if there is no 2D-PI in fact, because all the values are 1
        Par.Flags.TwoDCaipParallelImaging_flag = 0;
    end
    
    if(isfield(Par.Settings,'SliceAliasingPattern') && size(Par.Settings.SliceAliasingPattern,2) == 1)		% i.e. if there is no 2D-PI in fact, because no slices are aliased with each other
        Par.Flags.SliceParallelImaging_flag = 0;
    end
    
    
    
    %% 4. Get the B0map information
    
    if(Par.Flags.AlignFreq_flag && isfield(Par.Paths,'AlignFreq_path') )
        
        % AlignFreq_path is folder
        if(isempty(regexp(Par.Paths.AlignFreq_path{1},'\.mnc','ONCE')))
            CurFile = dir(Par.Paths.AlignFreq_path{1});
            CurFile = {CurFile.name};
            
            B0Map_Par = read_ascconv([Par.Paths.AlignFreq_path{1} '/' CurFile{3}]);
            
            Par.AlignFreq.DeltaTE_ms = diff(B0Map_Par.TEs)/1000;
            if(numel(Par.Paths.AlignFreq_path) > 1)                     % Phase unwrapping should be done --> Need different rescale factor
                Par.AlignFreq.RescaleFactor = 10^6/(2*pi*B0Map_Par.LarmorFreq);
            else
                Par.AlignFreq.RescaleFactor = 10^9/(8192*Par.AlignFreq.DeltaTE_ms*B0Map_Par.LarmorFreq);
            end
            fid = fopen([tmp_dir '/B0Map_Dummy.sh'],'w+');
            fprintf(fid,'AlignFreq_DeltaTE_ms=''%6.4f''; AlignFreq_RescaleFactor=''%9.8f'';',Par.AlignFreq.DeltaTE_ms,Par.AlignFreq.RescaleFactor);
            fclose(fid);
            
            clear fid CurFile
            
        end
        
        % Define size of B0MapTemplate.mnc
        if(isempty(regexpi(Par.Settings.AlignFreq_method,'Align')))     % If Align-method is chosen, use csi_template
            
            % Get Stepsize etc of B0-Map
            if(isempty(regexp(Par.Paths.AlignFreq_path{1},'\.mnc','ONCE')))
                AlignFreqPars.FoV_Read = B0Map_Par.FoV_Read(1); AlignFreqPars.FoV_Phase = B0Map_Par.FoV_Phase(1); AlignFreqPars.FoV_Partition = B0Map_Par.FoV_Partition;
                AlignFreqPars.nFreqEnc = B0Map_Par.nFreqEnc(1); AlignFreqPars.nPhasEnc = B0Map_Par.nPhasEnc(1); AlignFreqPars.nPartEnc = B0Map_Par.nSLC*B0Map_Par.nPartEnc;
                
                
            else         % AlignFreq_path is mnc file
                [deleteme,MincInfo] = unix(['mincinfo ' Par.Paths.AlignFreq_path{1}]);
                Par.AlignFreq.StepRead;
            end
            
            
            % Calculate Stepsize etc for resampled B0-Map
            if(Par.Flags.InterpolateCSIResolution_flag && ~Par.Settings.InterpolateCSIResolution_InkSpace)
                Par.AlignFreq.StepRead = Par.CSI.StepRead_BefInterpol/floor(abs(Par.CSI.StepRead_BefInterpol/(AlignFreqPars.FoV_Read/AlignFreqPars.nFreqEnc)));
                Par.AlignFreq.StepPhase = Par.CSI.StepRead_BefInterpol/floor(abs(Par.CSI.StepPhase_BefInterpol/(AlignFreqPars.FoV_Phase/AlignFreqPars.nPhasEnc)));
                Par.AlignFreq.StepSlice = Par.CSI.StepSlice_BefInterpol/floor(abs(Par.CSI.StepSlice_BefInterpol/(AlignFreqPars.FoV_Partition/AlignFreqPars.nPartEnc)));
                Pos_Minc_Dummy = Pos_Minc_BefInterpol;
            else
                Par.AlignFreq.StepRead = Par.CSI.StepRead/floor(abs(Par.CSI.StepRead/(AlignFreqPars.FoV_Read/AlignFreqPars.nFreqEnc)));
                Par.AlignFreq.StepPhase = Par.CSI.StepRead/floor(abs(Par.CSI.StepPhase/(AlignFreqPars.FoV_Phase/AlignFreqPars.nPhasEnc)));
                Par.AlignFreq.StepSlice = Par.CSI.StepSlice/floor(abs(Par.CSI.StepSlice/(AlignFreqPars.FoV_Partition/AlignFreqPars.nPartEnc)));
                Pos_Minc_Dummy = Pos_Minc;
            end
            
            
            Par.AlignFreq.POS_X_FirstVoxel = Pos_Minc_Dummy(1) + Par.AlignFreq.StepRead/2;
            Par.AlignFreq.POS_Y_FirstVoxel = Pos_Minc_Dummy(2) + Par.AlignFreq.StepPhase/2;
            Par.AlignFreq.POS_Z_FirstVoxel = Pos_Minc_Dummy(3);
            if(Par.CSI.ThreeD_flag)
                Par.AlignFreq.POS_Z_FirstVoxel = Par.AlignFreq.POS_Z_FirstVoxel + Par.CSI.StepSlice_BefInterpol/2;
            end
            
            Par.AlignFreq.nFreqEnc = round(abs(Par.CSI.FoV_Read / Par.AlignFreq.StepRead)); Par.AlignFreq.nPhasEnc = round(abs(Par.CSI.FoV_Phase / Par.AlignFreq.StepPhase));
            Par.AlignFreq.nPartEnc = round(abs(Par.CSI.FoV_Partition / Par.AlignFreq.StepSlice));
            
            Par.AlignFreq.PhaseNormalVector = Par.CSI.PhaseNormalVector; Par.AlignFreq.ReadNormalVector = Par.CSI.ReadNormalVector;
            Par.AlignFreq.SliceNormalVector = Par.CSI.SliceNormalVector;
            
            
        else
            if(Par.Flags.InterpolateCSIResolution_flag && ~Par.Settings.InterpolateCSIResolution_InkSpace)
                Par.AlignFreq.nFreqEnc = Par.CSI.nFreqEnc_BefInterpol;Par.AlignFreq.nPhasEnc = Par.CSI.nPhasEnc_BefInterpol;Par.AlignFreq.nPartEnc = Par.CSI.nPartEnc_BefInterpol;
            else
                Par.AlignFreq.nFreqEnc = Par.CSI.nFreqEnc;Par.AlignFreq.nPhasEnc = Par.CSI.nPhasEnc;Par.AlignFreq.nPartEnc = Par.CSI.nPartEnc;
            end
        end
        
    end
    
    
    %% 5. Get the imaging Information
    
    if(Par.Flags.image_normal_flag || Par.Flags.image_VC_flag || (isfield(Par.CSI,'WipMemBlockInterpretation') && isfield(Par.CSI.WipMemBlockInterpretation,'Prescan') && isfield(Par.CSI.WipMemBlockInterpretation.Prescan,'PATREFANDIMASCAN') && isfield(Par.CSI.WipMemBlockInterpretation.Prescan.PATREFANDIMASCAN,'nPhasEnc') && Par.CSI.WipMemBlockInterpretation.Prescan.PATREFANDIMASCAN.nPhasEnc > 0) )
        if(Par.Flags.image_normal_flag)
            imagingPars = read_ascconv(Par.Paths.image_normal_path{1});
            Par.Image.FoV_Read = imagingPars.FoV_Read;
            Par.Image.FoV_Phase = imagingPars.FoV_Phase;
            Par.Image.FoV_Partition = imagingPars.FoV_Partition;
            Par.Image.nFreqEnc = imagingPars.nFreqEnc;
            Par.Image.nPhasEnc = imagingPars.nPhasEnc;
            Par.Image.nPartEnc = imagingPars.nPartEnc;
            Par.Image.nSLC = imagingPars.nSLC;
        else
            Par.Image.FoV_Read = Par.CSI.FoV_Read;
            Par.Image.FoV_Phase = Par.CSI.FoV_Phase;
            Par.Image.FoV_Partition = Par.CSI.FoV_Partition;
            Par.Image.nFreqEnc = Par.CSI.WipMemBlockInterpretation.Prescan.PATREFANDIMASCAN.nPhasEnc;
            Par.Image.nPhasEnc = Par.CSI.WipMemBlockInterpretation.Prescan.PATREFANDIMASCAN.nPhasEnc;
            Par.Image.nPartEnc = Par.CSI.nPartEnc;
            Par.Image.nSLC = Par.CSI.WipMemBlockInterpretation.Prescan.PATREFANDIMASCAN.nSLC;
        end
        
        
        Par.Image.StepRead = -Par.Image.FoV_Read(1) / Par.Image.nFreqEnc;
        Par.Image.StepPhase = -Par.Image.FoV_Phase(1) / Par.Image.nPhasEnc;
        Par.Image.StepSlice = Par.Image.FoV_Partition(1) / (Par.Image.nPartEnc * Par.Image.nSLC);
        
        
        if(Par.Flags.image_normal_flag)
            if(~Par.Flags.phase_encoding_direction_is_RL_flag)
                Par.Image.FreqDirShift = round(-imagingPars.Pos_Sag/Par.Image.StepRead);
            else
                Par.Image.FreqDirShift = round(imagingPars.Pos_Cor/Par.Image.StepRead);    % Also a "-" necessary?
            end
        else
            Par.Image.FreqDirShift = 0;
        end
        
        clear imagingPars
        
    end
    
    
else
    
    Bak = Par;
    load(Par.Paths.csi_path{1},'Par')
    
    % Use only flags which are 1 in new processing, and 0 before. If its the other way round, produce an error, bc we cannot undo preprocessing.
    fieldy = fieldnames(Bak.Flags); 
    for ii = 1:numel(fieldy)
        TmpFlags = Bak.Flags.(fieldy{ii}) - Par.Flags.(fieldy{ii});
        if(TmpFlags < 0)
            fprintf('\n\nERROR: Re-Processing dataset, which has flag %s enabled, but should be reprocessed with flag disabled.',fieldy{ii})
            clearvars;  % So that program cannot continue in a meaningful way anymore.
            error('\nCannot undo pre-processing steps that were previously done.\n\n\n')
        end
        Par.Flags.(fieldy{ii}) = TmpFlags;
    end
    Par.Flags.mask_flag = Bak.Flags.mask_flag;
    Par.Flags.T1w_flag = Bak.Flags.T1w_flag;
    
    
	%Par.Flags = Bak.Flags;
    Par.Paths = Bak.Paths;
    Par.Settings = Bak.Settings;
    Par.LCMControl = Bak.LCMControl;
end

%% In case we have folder with dicom files
if(exist('BakCSIPath','var'))
    Par.Paths.csi_path = BakCSIPath;
end


%% 6. Print Out & Save Info

fprintf('\n\nI acquired the following information:\n');
Par %#ok
fprintf('\n\n\n')

% Pretty much complicated.
fprintf('\n\n\nPar.Paths\n')
printdummy = cell(2,1);
cnt = 0;
fn = transpose(fieldnames(Par.Paths));
for fn_dumm = fn
    cnt = cnt + 1;
    fn_dummy = fn_dumm{:};
    if(iscell(Par.Paths.(fn_dummy)))
        for cell_loop = 1:numel(Par.Paths.(fn_dummy))
            if(~(cell_loop == 1))
                cnt = cnt+1;
            end
            printdummy(:,cnt) = { fn_dummy, Par.Paths.(fn_dummy){cell_loop} };
        end
    else
        printdummy(:,cnt) = { fn_dummy, Par.Paths.(fn_dummy) };
    end
end
fprintf(['% ' num2str(max(cellfun(@numel,printdummy(1,:)))) 's: %s\n'],printdummy{:});


fprintf('\n\n\nPar.Flags')
Par.Flags
fprintf('\n\n\nPar.Settings')
Par.Settings
if(isfield(Par,'AlignFreq'))
    fprintf('\n\n\nPar.AlignFreq')
    Par.AlignFreq
end
if(isfield(Par,'LCMControl'))
    fprintf('\n\n\nPar.LCMControl')
    Par.LCMControl
end
if(isfield(Par,'ServerInfo'))
    fprintf('\n\n\nPar.ServerInfo')
    Par.ServerInfo
end
fprintf('\n\n\nPar.CSI')
Par.CSI
if(isfield(Par,'Image'))
    fprintf('\n\n\nPar.Image')
    Par.Image
end

clearvars -except tmp_dir Par ErrorFile

if(~exist('Bak','var'))
    if (strcmp(Par.Settings.WaterReference_method, 'W2') && ~exist([Par.Paths.out_path '/Parameters_water.mat' ],'file'))
        save([Par.Paths.out_path '/Parameters_water.mat'],'Par');
        save([tmp_dir '/Parameters_water.mat'],'Par');
    else
        save([Par.Paths.out_path '/Parameters.mat'],'Par');
        save([tmp_dir '/Parameters.mat'],'Par');
    end
end


%% 7. Create Minc Template

fprintf('\n\nCreate Minc Templates\n');
Create_MincTemplates(tmp_dir, Par)





%% 8. Matlab Part of Creating Masks





if(Par.Flags.mask_flag)
    
    if(~isempty(regexpi(Par.Settings.mask_method, 'voi')) || Par.CSI.ThreeD_flag)
        fprintf('\n\nCreate VoI Mask\n');
        create_mask_VOI(tmp_dir)
    end
    
    
else        % PROCESS WHOLE FoV
    
    magnitude_mask = ones([Par.CSI.nFreqEnc Par.CSI.nPhasEnc Par.CSI.nPartEnc*Par.CSI.nSLC]);
    magnitude_fid = fopen([tmp_dir '/mask_brain.raw'],'w');
    fwrite(magnitude_fid,magnitude_mask,'float');
    fclose(magnitude_fid);
    
    if(Par.CSI.ThreeD_flag)
        ZFSize = [Par.CSI.nFreqEnc Par.CSI.nPhasEnc Par.CSI.nPartEnc*Par.CSI.nSLC] * Par.Settings.ZeroFillMetMaps;
    else
        ZFSize = [Par.CSI.nFreqEnc*Par.Settings.ZeroFillMetMaps Par.CSI.nPhasEnc*Par.Settings.ZeroFillMetMaps Par.CSI.nPartEnc*Par.CSI.nSLC];
    end
    
    
    magnitude_mask = ones(ZFSize);
    magnitude_fid = fopen([tmp_dir '/mask_brain_zf.raw'],'w');
    fwrite(magnitude_fid,magnitude_mask,'float');
    fclose(magnitude_fid);
    
    if(Par.Flags.InterpolateCSIResolution_flag == 1)
        magnitude_mask_BefInterpol = ones([Par.CSI.nFreqEnc_BefInterpol Par.CSI.nPhasEnc_BefInterpol Par.CSI.nPartEnc_BefInterpol*Par.CSI.nSLC_BefInterpol]);
        magnitude_fid_BefInterpol = fopen([tmp_dir '/mask_brain_BefInterpol.raw'],'w');
        fwrite(magnitude_fid_BefInterpol,magnitude_mask_BefInterpol,'float');
        fclose(magnitude_fid_BefInterpol);
    end
    
end


if (~Par.Flags.T1w_flag)
    fprintf('\n\nCreate Magnitude of Image, Image_VC or CSI for Masking.\n');
    create_magnitude(tmp_dir)
end



fprintf(ErrorFile,'ErrorInGetPar_CreateTempl=0; ErrorMessage=''''');
fclose(ErrorFile);

end

function Info = ReadMeasurementInfoFromMRSI(Par)

    fprintf('Saving MeasurementInfos.txt to tmp from %s\n', Par.Paths.csi_path{1});
    if(Par.Flags.DICOM_flag)
        fprintf('\nWarning: Input were DICOM files. Not all info can be read from DICOMs. Update the code in GetPar_CreateTempl_MaskPart1.m.')
        
        Test = dicominfo(Par.Paths.csi_path{1});

%         FamName = Test.PatientName.FamilyName;
%         if(isempty(FamName))
%             FamName = 'xxxxx';
%         end
        if(isempty(Test.PatientBirthDate))
            Test.PatientBirthDate = 'Unknown';
        end        
        if(isempty(Test.PatientAge))
            Test.PatientAge = 'Unknown';
        end 
        if(isempty(Test.PatientAge))
            Test.PatientAge = 'Unknown';
        end 
        if(isempty(Test.PatientSex))
            Test.PatientSex = 'Unknown';
        end 
        if(isempty(Test.PatientWeight))
            Test.PatientWeight = 'Unknown';
        end         
        if(isempty(Test.PatientSex))
            Test.PatientSex = 'Unknown';
        end 
        if(isempty(Test.PatientSex))
            Test.PatientSex = 'Unknown';
        end 

        Info=[Par.CSI.PatName, ' ', ...
            num2str(Test.PatientBirthDate), ' ', ...
            num2str(Test.PatientAge), ' ', ...
            num2str(Test.PatientSex), ' ', ...
            num2str(Test.PatientWeight), ' ', ...
            'Unknown', ' ', ...
            'Unknown', ' ', ...
            'Unknown', ' ', ...
            'Unknown'];    
    else
        mapVBVDHdr = read_twix_hdr_standalone(Par.Paths.csi_path{1});
        mapVBVDHdr = mapVBVDHdr{end};
        Info=[mapVBVDHdr.Dicom.tPatientName, ' ', ...
            num2str(mapVBVDHdr.Config.PatientBirthDay), ' ', ...
            num2str(mapVBVDHdr.Dicom.flPatientAge), ' ', ...
            num2str(mapVBVDHdr.Dicom.lPatientSex), ' ', ...
            num2str(mapVBVDHdr.Dicom.flUsedPatientWeight), ' ', ...
            num2str(mapVBVDHdr.Dicom.lFrequency), ' ', ...
            num2str(mapVBVDHdr.Dicom.flMagneticFieldStrength), ' ', ...
            num2str(mapVBVDHdr.Dicom.flTransRefAmpl), ' ', ...
            num2str(mapVBVDHdr.Dicom.adFlipAngleDegree)];
    end
end
