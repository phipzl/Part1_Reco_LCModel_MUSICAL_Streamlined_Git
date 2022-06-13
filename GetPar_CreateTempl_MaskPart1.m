%% 0. Preparations

clear functions; close all;
clearvars -except tmp_dir

run([tmp_dir '/InitialParameters.m'])




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
    eval(['run ' Par.Paths.LCM_ControlPath]);
    Par.LCMControl = ControlWrite;
end


clearvars -except tmp_dir Par ErrorFile

% Write LCM Control file for water dataset
if(isfield(Par.Paths,'LCM_Control_Water_path'))
    eval(['run ' Par.Paths.LCM_Control_Water_path]);
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

% if(Par.Flags.LipidDecon_flag == 1)
% 	% L1 or L2 ?
% 	Par.Flags.LipidDecon_L1_flag = ~isempty(regexp(Par.Settings.LipidDecon_MethodAndNoOfLoops,'L1','ONCE'));
% 	dummy = regexp(Par.Settings.LipidDecon_MethodAndNoOfLoops,',[0-9]*\.?[0-9]*','match');
% 	if(Par.Flags.LipidDecon_L1_flag)
% 		if(isempty(dummy))
% 			Par.Settings.LipidDecon_NoOfLoops = 10;
% 		else
% 			Par.Settings.LipidDecon_NoOfLoops = str2double(dummy{1}(2:end));			
% 		end
% 	else
% 		if(isempty(dummy))
% 			Par.Settings.LipidDecon_L2BetaCorrFactor = 1;
% 		else
% 			Par.Settings.LipidDecon_L2BetaCorrFactor = str2double(dummy{1}(2:end));
% 		end
% 	end
% end

if(Par.Flags.LipidDecon_flag == 1)
	dummy = regexp(Par.Settings.LipidDecon_L2BetaCorrFactor,',','split');
    dummy_factor = dummy(end);
    if(isempty(dummy_factor))
        Par.Settings.LipidDecon_L2BetaCorrFactor = 1;
    else
        Par.Settings.LipidDecon_L2BetaCorrFactor = str2double(dummy_factor{1}(2:end));
    end
end

Par.Flags.zerofill_to_nextpow2_flag = false;
if(Par.Flags.InterpolateCSIResolution_flag == 1)
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
fprintf('\n\nWater Reference Processing\n');

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
        end
    else
        warning('Not a correct format for water referencing, setting Water flag to zero');
        Par.Flags.WaterReference_flag = 0;
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
	ServerSSH = [Par.ServerInfo.RunLCModelAs '@' Par.ServerInfo.RunLCModelOn]; ServerSSH(ServerSSH(1) == '@') = []; % Delete @ symbol if it is on first place
end
TestLocal = setdiff(fields_files,TestOnServer(1:end-1));

% Test for existence
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


%% 3. Get the CSI Information

if(isempty(regexp(Par.Paths.csi_path{1},'.*/\w*\.mat','ONCE')))			% Only read header of first csi file. The others must be the same

    if(strcmp(Par.Settings.WaterReference_method, 'W2') && ~exist([Par.Paths.out_path '/Parameters_water.mat' ],'file'))
        Par.Paths.csi_path{1} = Par.Paths.WaterReference_path;
        fprintf('\n\nParameters_water.mat are being created.\n\n')
    else
        fprintf('\n\nParameters.mat are being created.\n\n')
    end
    %csiPars = read_ascconv(Par.Paths.csi_path{1});		% New function which is adapted to VE11. If some bug occurs, search here.
    inData=mapVBVD(Par.Paths.csi_path{1});
    if(size(inData,2)==2)
        inData=inData{2};
    end
    
	% Total Channel Number
	 if(strcmpi(inData.image.softwareVersion,'vd'))
         csiPars.total_channel_no = numel(inData.hdr.MeasYaps.sCoilSelectMeas.aRxCoilSelectData{1}.asList);
     end
     
     if(strcmpi(inData.image.softwareVersion,'vb'))
        csiPars.total_channel_no = numel(inData.hdr.MeasYaps.asCoilSelectMeas{1}.asList);
     end

    for i = 1 : size(inData.hdr.MeasYaps.alTE,2)
        csiPars.TEs(i)=inData.hdr.MeasYaps.alTE{i};
    end
    
    csiPars.nFreqEnc = inData.hdr.Config.PhaseEncodingLines;
    csiPars.nPhasEnc = inData.hdr.Config.PhaseEncodingLines;
    csiPars.nPartEnc = inData.hdr.Spice.Partitions;
    csiPars.nSLC = 1;
	% Define Matrix Size	
	dummy = [csiPars.nFreqEnc csiPars.nPhasEnc csiPars.nPartEnc*csiPars.nSLC];
	
	Par.Settings.InterpolateCSIResolution(Par.Settings.InterpolateCSIResolution == [0 0 0]) = dummy(Par.Settings.InterpolateCSIResolution == [0 0 0]);
	Par.Settings.InterpolateResolutionRatio = [csiPars.nFreqEnc csiPars.nPhasEnc csiPars.nPartEnc*csiPars.nSLC] ./ Par.Settings.InterpolateCSIResolution;
	
	if(sum(Par.Settings.InterpolateResolutionRatio > 1) >= 1 && ~Par.Flags.AlignFrequency_flag)
		fprintf('\n\nWarning: Some spectra will be summed, but no frequency alignment was requested. I set AlignFrequency_flag = true.\n\n')
		Par.Flags.AlignFrequency_flag = true;
	end
	if(numel(Par.Paths.csi_path) > 1 && ~Par.Flags.AlignFrequency_flag)
		fprintf('\n\nWarning: Perform averaging of several csi data sets, but no frequency alignment was requested. I set AlignFrequency_flag = true.\n\n')
		Par.Flags.AlignFrequency_flag = true;
	end	
	
	

	
	csiPars.nFreqEnc = csiPars.nFreqEnc / Par.Settings.InterpolateResolutionRatio(1);
	csiPars.nPhasEnc = csiPars.nPhasEnc / Par.Settings.InterpolateResolutionRatio(2);
	csiPars.nPartEnc = csiPars.nPartEnc / Par.Settings.InterpolateResolutionRatio(3);
	csiPars.nFreqEnc_nonzf = csiPars.nFreqEnc * Par.Settings.InterpolateResolutionRatio(1);
	csiPars.nPhasEnc_nonzf = csiPars.nPhasEnc * Par.Settings.InterpolateResolutionRatio(2);
	csiPars.nPartEnc_nonzf = csiPars.nPartEnc * Par.Settings.InterpolateResolutionRatio(3);
	csiPars.nFreqEnc_FinalMatrix = csiPars.nFreqEnc / Par.Settings.InterpolateResolutionRatio(1);
	csiPars.nPhasEnc_FinalMatrix = csiPars.nPhasEnc / Par.Settings.InterpolateResolutionRatio(2);
	csiPars.nSLC_FinalMatrix = csiPars.nSLC / Par.Settings.InterpolateResolutionRatio(3);
	
	

	
	% DICOM is always zerofilled to next power of 2.
	if(Par.Flags.zerofill_to_nextpow2_flag || numel(strfind(Par.Paths.csi_path{1}, '.IMA')) > 0)
		csiPars.nFreqEnc = csiPars.nFreqEnc_FinalMatrix;
		csiPars.nPhasEnc = csiPars.nPhasEnc_FinalMatrix;    
		if(csiPars.ThreeD_flag)
			csiPars.nPartEnc = csiPars.nSLC_FinalMatrix;    % Only in ThreeD Case we have to change that.
		end
	end
	%csiPars.nFreqEnc = csiPars.nFreqEnc *2;        uncomment if zerofill wanted
	%csiPars.nPhasEnc = csiPars.nPhasEnc *2;
    csiPars.ThreeD_flag=0;
    if(csiPars.nSLC>1)
        csiPars.ThreeD_flag=1;
    end
	csiPars = rmfield(csiPars,{'nFreqEnc_FinalMatrix','nPhasEnc_FinalMatrix','nSLC_FinalMatrix' });
	if(~csiPars.ThreeD_flag && numel(strfind(Par.Paths.csi_path{1}, '.IMA')) > 0)
		csiPars.nSLC = 1;
    end

    
    csiPars.FoV_Read = inData.hdr.Config.PhaseFoV;
    csiPars.FoV_Phase = inData.hdr.Config.PhaseFoV;
    
    csiPars.FoV_Partition =  inData.hdr.MeasYaps.sSliceArray.asSlice{1}.dThickness;
    csiPars.VoI_Partition =  inData.hdr.Phoenix.sSpecPara.sVoI.dThickness;

    if(csiPars.nSLC==1 && csiPars.nPartEnc==1)
        % Stepsize (Voxel Size)%lukas master hack (acc to stano)
        csiPars.StepRead = -csiPars.FoV_Read(1) / csiPars.nFreqEnc;		% Coordinate system is reversed in minc with respect to DICOM
        csiPars.StepPhase = -csiPars.FoV_Phase(1) / csiPars.nPhasEnc;    
        csiPars.StepSlice = csiPars.VoI_Partition(1) / (csiPars.nPartEnc * csiPars.nSLC); 
    else
        
        % Stepsize (Voxel Size)
        csiPars.StepRead = -csiPars.FoV_Read(1) / csiPars.nFreqEnc;		% Coordinate system is reversed in minc with respect to DICOM
        csiPars.StepPhase = -csiPars.FoV_Phase(1) / csiPars.nPhasEnc;    
        csiPars.StepSlice = csiPars.FoV_Partition(1) / (csiPars.nPartEnc * csiPars.nSLC); 
    end
	
	if(~isfield(Par,'Settings') || ~isfield(Par.Settings,'ZeroFillMetMaps'))
		Par.Settings.ZeroFillMetMaps = 4;
	end
	zff = Par.Settings.ZeroFillMetMaps;

	
	csiPars.StepRead_zf = csiPars.StepRead / zff;
	csiPars.StepPhase_zf = csiPars.StepPhase / zff;
	if(csiPars.ThreeD_flag)
		csiPars.StepSlice_zf = csiPars.StepSlice / zff; % TODO  this should be implemented later, when zerofilling is done in z direction - csiPars.StepSlice_zf = csiPars.StepSlice / 2;
	else
		csiPars.StepSlice_zf = csiPars.StepSlice;
	end

    csiPars.Pos_Tra = inData.image.slicePos(3,1);
    csiPars.Pos_Sag = inData.image.slicePos(1,1);
    csiPars.Pos_Cor = inData.image.slicePos(2,1);
	% Rename Position Fields
	[csiPars.POS_X] = csiPars.Pos_Sag;
	[csiPars.POS_Y] = csiPars.Pos_Cor;
	[csiPars.POS_Z] = csiPars.Pos_Tra;
    csiPars.PosVOI_Sag = csiPars.Pos_Sag;
    csiPars.PosVOI_Cor = csiPars.Pos_Cor;
    csiPars.PosVOI_Tra = csiPars.Pos_Tra;
    csiPars.SliceGap = 0;
    csiPars.SliceThickness = csiPars.FoV_Partition;
    
	csiPars = rmfield(csiPars,{'Pos_Sag','Pos_Cor','Pos_Tra'});

    csiPars.SliceNormalVector_x(1)=0;
    if(isfield(inData.hdr.Phoenix.sSpecPara.sVoI.sNormal,'dSag'))
        csiPars.SliceNormalVector_x(1)=inData.hdr.Phoenix.sSpecPara.sVoI.sNormal.dSag;
    end
    
    csiPars.SliceNormalVector_y(1)=0;
    if(isfield(inData.hdr.Phoenix.sSpecPara.sVoI.sNormal,'dCor'))
        csiPars.SliceNormalVector_y(1)=inData.hdr.Phoenix.sSpecPara.sVoI.sNormal.dCor;
    end
    
    csiPars.SliceNormalVector_z(1)=0;
    if(isfield(inData.hdr.Phoenix.sSpecPara.sVoI.sNormal,'dTra'))
        csiPars.SliceNormalVector_z(1)=inData.hdr.Phoenix.sSpecPara.sVoI.sNormal.dTra;
    end
    
    csiPars.InPlaneRotation=0;
    csiPars.InPlaneRotation_VOI=0;
    if(~isempty(inData.hdr.Spice.VoiInPlaneRot))
       csiPars.InPlaneRotation = inData.hdr.Spice.VoiInPlaneRot; 
       csiPars.InPlaneRotation_VOI = inData.hdr.Spice.VoiInPlaneRot; 
    end
        
    csiPars.InPlaneRotation=csiPars.InPlaneRotation-0.1;
    csiPars.InPlaneRotation_VOI=csiPars.InPlaneRotation_VOI-0.1;

	% Compute direction cosine from x,y and z components of slice normal vector
	[csiPars.PhaseNormalVector, csiPars.ReadNormalVector] = compute_dircos_1_3([csiPars.SliceNormalVector_x(1) csiPars.SliceNormalVector_y(1) csiPars.SliceNormalVector_z(1)],csiPars.InPlaneRotation);
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
	FoVHalf = [csiPars.FoV_Read(1)/2 csiPars.FoV_Phase(1)/2 -csiPars.FoV_Partition(1)/csiPars.nPartEnc*(csiPars.nPartEnc-1)/2];
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
		
		FoVHalf_BefInterpol = FoVHalf; FoVHalf_BefInterpol(3) = -csiPars.FoV_Partition(1)/csiPars.nPartEnc_BefInterpol*(csiPars.nPartEnc_BefInterpol-1)/2;
		Pos_Minc_BefInterpol = Pos_Minc - FoVHalf + FoVHalf_BefInterpol;
		
		csiPars.POS_X_FirstVoxel_BefInterpol = Pos_Minc_BefInterpol(1) + csiPars.StepRead_BefInterpol/2;
		csiPars.POS_Y_FirstVoxel_BefInterpol = Pos_Minc_BefInterpol(2) + csiPars.StepPhase_BefInterpol/2;
		csiPars.POS_Z_FirstVoxel_BefInterpol = Pos_Minc_BefInterpol(3);	
		
		if(csiPars.ThreeD_flag)
			csiPars.POS_Z_FirstVoxel_BefInterpol = csiPars.POS_Z_FirstVoxel_BefInterpol + csiPars.StepSlice_BefInterpol/2;     
		end

	end
	
	

	% Read Patientname
	if(numel(strfind(Par.Paths.csi_path{1}, '.IMA')) > 0)
		[bla, csiPars.PatName] = unix(['dcmdump +P "0010,0010" ' Par.Paths.csi_path{1}] );
		csiPars.PatName = regexp(csiPars.PatName,'\[(?!^).*\^?(?!^).*\]','match');
		if(bla > 0)		% If dcmdump doesnt exist, search for the Patient field manually. One day this probably should be directly implemented in Matlab, without using bash via unix command
			[bla, Border1] = unix(['grep --color=''never'' -o -b -u -P -a ''\x10\x00\x10\x00'' ' Par.Paths.csi_path{1}]);			% Search for byte-offset of field 10001000, which is Patient name
			Border1 = regexp(Border1,'\d+:','match'); Border1 = str2double(Border1{1}(1:end-1));				% Get the byte-offset
			[bla, Border2] = unix(['grep --color=''never'' -o -b -u -P -a ''\x10\x00\x20\x00'' ' Par.Paths.csi_path{1}]);			% Search for byte-offset of field 10002000, which is the following field
			Border2 = regexp(Border2,'\d+:','match'); Border2 = str2double(Border2{1}(1:end-1));
			[bla, csiPars.PatName] = unix(['cat ' Par.Paths.csi_path{1} ' | head -c ' num2str(Border2) ' | tail -c ' num2str(Border2-Border1-8)]);	% Cut everything out between those fields
			csiPars.PatName = cellstr(csiPars.PatName);
		end


		if(isempty(csiPars.PatName))
			csiPars.PatName = cellstr('NoName');		% E.g. if dcmdump does not exist.
		end
		csiPars.PatName = regexprep(csiPars.PatName{:},{'\[','\]','\^'},{'','','_'});
	else
		[bla, csiPars.PatName] = unix(['strings ' Par.Paths.csi_path{1} ' | grep -m 1 "tPatientName"']);		% Example Result: <ParamString."tPatientName"> { "Strasser^Bernhard" }
		csiPars.PatName = regexp(csiPars.PatName,'{ "(?!^).*\^?(?!^).*"  }','match');                       % (?!^).*: Any character 0 or more times, but no caret, \^?: ^ 0 or 1 times. --> { "Strasser^Bernhard" }
		csiPars.PatName = regexp(csiPars.PatName{:},'".*"','match');                                        % --> "Strasser^Bernhard"
		if(isempty(csiPars.PatName))
			csiPars.PatName = cellstr('NoName');		% E.g if dcmdump does not exist.
		end
		csiPars.PatName = regexprep(csiPars.PatName{:},{'"','\^'},{'','_'});                                % Replace " with nothing and caret with underscore. --> Strasser_Bernhard  
	end
	

	csiPars.PatName = regexprep(csiPars.PatName,' ','_');



	Par.CSI = csiPars; clear csiPars





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

% 	if(0)
% 		EllipFilt = EllipticalFilter(ones([Par.CSI.nFreqEnc Par.CSI.nPhasEnc]),[1 2],[1 1 1 floor(Par.CSI.nFreqEnc/2)-1],1);
% 		nSpatEncCSI = sum(sum(EllipFilt)); clear EllipFilt;
% 	else
% 		nSpatEncCSI = Par.CSI.nFreqEnc*Par.CSI.nPhasEnc;
% 	end
% 	size_csidata_predicted = Par.CSI.total_channel_no * nSpatEncCSI * Par.CSI.nPartEnc * Par.CSI.nSLC * Par.CSI.vecSize * 4 * 2;  % 4: float32, 2: real imag
% 
% 	size_csidata_ratio = size_csidata_predicted / size_csidata;






	% Set the Parallel Imaging Flags and values, if info is written in the header of the CSI data
	% 1D Caipi. Only if the size ratio is above 1.7
% 	if(size_csidata_ratio > 1.7)
% 		% Check if info is available in wipmemblock
% 		OneDCaipInfoAvail = isfield(Par.CSI,'WipMemBlockInterpretation') && isfield(Par.CSI.WipMemBlockInterpretation,'OneDCaipi') && isstruct(Par.CSI.WipMemBlockInterpretation.OneDCaipi) ...
% 		&& Par.CSI.WipMemBlockInterpretation.OneDCaipi.SliceParallelImaging_flag;
% 		% Check if the user gave info. Dont touch user given info.
% 		if(OneDCaipInfoAvail && ~isfield(Par.Settings,'SliceAliasingPattern') && ~isfield(Par.Settings,'FoVShifts_x') && ~isfield(Par.Settings,'FoVShifts_y') )
% 			Par.Settings.SliceAliasingPattern = Par.CSI.WipMemBlockInterpretation.OneDCaipi.SliceAliasingPattern;	
% 			Par.Settings.FoVShifts_x = Par.CSI.WipMemBlockInterpretation.OneDCaipi.FoVShifts_x;
% 			Par.Settings.FoVShifts_y = Par.CSI.WipMemBlockInterpretation.OneDCaipi.FoVShifts_y;
% 			Par.Flags.SliceParallelImaging_flag = Par.CSI.WipMemBlockInterpretation.OneDCaipi.SliceParallelImaging_flag;
% 		end
% 	end
% 	% Make SliceAliasing Patterns like [1 2], [1 3; 2 4], [1 4; 2 5;3 6] ...
% 	if(Par.Flags.SliceParallelImaging_flag && ~isfield('Par.Settings','SliceAliasingPattern'))
% 		for AliLoop = 1:numel(Par.Settings.FoVShifts_x)/2
% 			Par.Settings.SliceAliasingPattern(AliLoop,:) = [AliLoop AliLoop+numel(Par.Settings.FoVShifts_x)/2];
% 		end
% 	end
% 
% 	% Calculate ratio again with new info
% 	if(Par.Flags.SliceParallelImaging_flag)
% 		size_csidata_predicted = size_csidata_predicted / size(Par.Settings.SliceAliasingPattern,2);
% 	end
% 	size_csidata_ratio = size_csidata_predicted / size_csidata;
% 
% 
% 	% 2D Caipi. Only if the size ratio is above 1.4
% 	if(size_csidata_ratio > 1.4)
% 		TwoDCaipInfoAvail = isfield(Par.CSI,'WipMemBlockInterpretation') && isfield(Par.CSI.WipMemBlockInterpretation,'TwoDCaipi') && isstruct(Par.CSI.WipMemBlockInterpretation.TwoDCaipi);
% 		
% 		TwoDCaipInfoDefective = ~isstruct(Par.CSI.WipMemBlockInterpretation.TwoDCaipi) && Par.CSI.WipMemBlockInterpretation.TwoDCaipi == -1;
% 		if(TwoDCaipInfoDefective && ~isfield(Par.Settings,'InPlaneCaipPattern'))
% 			fprintf('\nERROR: 2D-CAIPI seems to have been performed, but I could not read the Pattern.')
% 			fprintf('\nInput the 2D-Pattern manually with option r\n(e.g.: -r ''InPlaneCaipPattern = zeros([5 5]); InPlaneCaipPattern([4 7 9 15 18 21]) = 1; VD_Radius = 1;'').\nStopping here.')
% 			
% 			fprintf(ErrorFile,'ErrorInGetPar_CreateTempl=1; ErrorMessage=''TwoDCaipiInfo defective in ascconv header.''');
% 			fclose(ErrorFile);	
% 			return;
% 		end
% 		
% 		
% 		if(TwoDCaipInfoAvail && ~isfield(Par.Settings,'InPlaneCaipPattern') && ~isfield(Par.Settings,'VD_Radius'))
% 			Par.Settings.InPlaneCaipPattern = Par.CSI.WipMemBlockInterpretation.TwoDCaipi.Skip_Matrix;
% 			Par.Settings.VD_Radius = Par.CSI.WipMemBlockInterpretation.TwoDCaipi.VD_Radius;
% 			Par.Flags.TwoDCaipParallelImaging_flag = Par.CSI.WipMemBlockInterpretation.TwoDCaipi.TwoDCaipParallelImaging_flag;
% 		end
% 	end
% 
% 
% 
% 
% 
% 	if(isfield(Par.Settings,'InPlaneCaipPattern') && sum(sum(Par.Settings.InPlaneCaipPattern)) == numel(Par.Settings.InPlaneCaipPattern) ...
% 		|| isfield(Par.Settings,'VD_Radius') && Par.Settings.VD_Radius > floor(Par.CSI.nFreqEnc/2)-1)		% i.e. if there is no 2D-PI in fact, because all the values are 1
% 		Par.Flags.TwoDCaipParallelImaging_flag = 0;
% 	end
% 
% 	if(isfield(Par.Settings,'SliceAliasingPattern') && size(Par.Settings.SliceAliasingPattern,2) == 1)		% i.e. if there is no 2D-PI in fact, because no slices are aliased with each other
% 		Par.Flags.SliceParallelImaging_flag = 0;
% 	end







	%% 4. Get the imaging Information

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
	Par.Flags = Bak.Flags;
	Par.Paths = Bak.Paths;
	Par.Settings = Bak.Settings;
	Par.LCMControl = Bak.LCMControl;
end

%% 5. Print Out & Save Info

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
            if(~cell_loop == 1)
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
if(isfield(Par,'LCMControl'))
    fprintf('\n\n\nPar.LCMControl')
    Par.LCMControl
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

clear inData
%% 6. Create Minc Template

fprintf('\n\nCreate Minc Templates\n');
run ./Create_MincTemplates.m





%% 7. Matlab Part of Creating Masks





if(Par.Flags.mask_flag)
    
    if(~isempty(regexpi(Par.Settings.mask_method, 'voi'))||~isempty(regexpi(Par.Settings.mask_method, 'dreid')))
        fprintf('\n\nCreate VoI Mask\n');    
        run ./create_mask_VOI.m
    end
        
    
else        % PROCESS WHOLE FoV
    
    magnitude_mask = ones([Par.CSI.nFreqEnc Par.CSI.nPhasEnc Par.CSI.nPartEnc*Par.CSI.nSLC]);
    magnitude_fid = fopen([tmp_dir '/mask_brain.raw'],'w');
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
    fprintf('\n\n\n\nWARNING:\nSKIPPING DUE TO PAST ISSUES!!!\nPL20220502\n\n\n\n');
%    run ./create_magnitude.m
end




fprintf(ErrorFile,'ErrorInGetPar_CreateTempl=0; ErrorMessage=''''');
fclose(ErrorFile);		   

