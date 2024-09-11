%% -1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%                 CREATE a MASK TO PROCESS ONLY DATA INSIDE THE VOI              %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function create_mask_VOI(tmp_dir)
%% 0. DEFINITIONS, PREPARATIONS

if ~exist([tmp_dir '/Parameters.mat'],'file')
	fprintf('Parameters.mat not found. Maybe this is water scan.\n');
	return
end

load([tmp_dir '/Parameters.mat'])

fn = transpose(fieldnames(Par.CSI));
for fn_dummy = fn
    eval([fn_dummy{:} ' = Par.CSI.' fn_dummy{:} ';']);
end
clear fn fn_dummy


%magnitude_csi = read_csi_1_4(csi_path,0,1,1);		%eigentlich unnötig CSI einzulesen ... ich weiss alles von Parameter.m
%magnitude_csi = squeeze(magnitude_csi(:,:,:,:,1));

magnitude_mask = zeros([nFreqEnc nPhasEnc nPartEnc*nSLC]);
% VoI_Read = FoV_Read;		% IS THIS NECESSARY???
% VoI_Phase = FoV_Phase;


%% 1. CREATE MASK DIRECTLY based on known VOI size

% Use function to calculate edges of VoI. The VoI is always in the center of the FoV. Therefore, if we have an even number of voxels, the VoI can 
% only fully cover an even number of Voxels. Example: Have 10 voxels, VoI is size of 2 voxels. Center of VoI is at 5.5, covering FoV from 4.5 to 6.5,
% which is one full voxel (the central from 5.0 to 6.0), and two half-voxels (from 4.5 to 5.0 and 6.0 to 6.5).
% Therefore, round it up to even number of voxels in that case.
% For odd number of voxels, the VoI always covers odd number of voxels (bc the VoI-center is at the FoV-center).
ExtraVoxels = 1;
[VOI_X_MIN,VOI_X_MAX] = CalcVoIEdges(nFreqEnc,VoI_Read(1),FoV_Read(1),ExtraVoxels);
[VOI_Y_MIN,VOI_Y_MAX] = CalcVoIEdges(nPhasEnc,VoI_Phase(1),FoV_Phase(1),ExtraVoxels);
[VOI_Z_MIN,VOI_Z_MAX] = CalcVoIEdges(nPartEnc*nSLC,VoI_Partition,FoV_Partition,ExtraVoxels);


if(Par.CSI.ThreeD_flag && ~isempty(regexpi(Par.Settings.mask_method, 'bet')))               % If in fact we want bet, but have a 3D data set, we want to restrict it by VoI_z
    magnitude_mask(:,:,VOI_Z_MIN:VOI_Z_MAX) = 1;    
else
    magnitude_mask(VOI_X_MIN:VOI_X_MAX,VOI_Y_MIN:VOI_Y_MAX,VOI_Z_MIN:VOI_Z_MAX) = 1;
end

% LUKI:
% if (VOI_X_MIN==0 && VOI_Y_MIN==0 && VOI_Z_MIN==0)%there is sum rounding bug...lukas
%   magnitude_mask(VOI_X_MIN+1:VOI_X_MAX-1,VOI_Y_MIN+1:VOI_Y_MAX-1,VOI_Z_MIN+1:VOI_Z_MAX-1) = 1; 
% elseif (VOI_X_MIN==0 && VOI_Y_MIN==0)
%     magnitude_mask(VOI_X_MIN+1:VOI_X_MAX-1,VOI_Y_MIN+1:VOI_Y_MAX-1,VOI_Z_MIN:VOI_Z_MAX) = 1; 
% else
%     magnitude_mask(VOI_X_MIN:VOI_X_MAX,VOI_Y_MIN:VOI_Y_MAX,VOI_Z_MIN:VOI_Z_MAX) = 1;
% end

if(Par.Flags.InterpolateCSIResolution_flag == 1)
	magnitude_mask_BefInterpol = zeros([nFreqEnc_BefInterpol nPhasEnc_BefInterpol nPartEnc_BefInterpol*nSLC_BefInterpol]);
	
    
    [VOI_X_MIN_BefInterpol,VOI_X_MAX_BefInterpol] = CalcVoIEdges(nFreqEnc_BefInterpol,VoI_Read(1),FoV_Read(1),ExtraVoxels);
    [VOI_Y_MIN_BefInterpol,VOI_Y_MAX_BefInterpol] = CalcVoIEdges(nPhasEnc_BefInterpol,VoI_Phase(1),FoV_Phase(1),ExtraVoxels);
    [VOI_Z_MIN_BefInterpol,VOI_Z_MAX_BefInterpol] = CalcVoIEdges(nPartEnc_BefInterpol*nSLC_BefInterpol,VoI_Partition,FoV_Partition,ExtraVoxels);

    
    if(Par.CSI.ThreeD_flag && ~isempty(regexpi(Par.Settings.mask_method, 'bet')))               % If in fact we want bet, but have a 3D data set, we want to restrict it by VoI_z
        magnitude_mask_BefInterpol(:,:,VOI_Z_MIN_BefInterpol:VOI_Z_MAX_BefInterpol) = 1;
    else
	magnitude_mask_BefInterpol(VOI_X_MIN_BefInterpol:VOI_X_MAX_BefInterpol,VOI_Y_MIN_BefInterpol:VOI_Y_MAX_BefInterpol,VOI_Z_MIN_BefInterpol:VOI_Z_MAX_BefInterpol) = 1;
    end
    
	magnitude_fid_BefInterpol = fopen([tmp_dir '/mask_brain_VOI_BefInterpol.raw'],'w');
	fwrite(magnitude_fid_BefInterpol,magnitude_mask_BefInterpol,'float');
	fclose(magnitude_fid_BefInterpol);	
end

%% 2. WRITE DATA AS .RAW-FILES

magnitude_fid = fopen([tmp_dir '/mask_brain_VOI.raw'],'w');
fwrite(magnitude_fid,magnitude_mask,'float');
fclose(magnitude_fid);

end
