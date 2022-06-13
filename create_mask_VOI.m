%% -1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%                 CREATE a MASK TO PROCESS ONLY DATA INSIDE THE VOI              %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 0. DEFINITIONS, PREPARATIONS

%clearvars -except tmp_dir; clear functions; close all;
%load([tmp_dir '/Parameters.mat'])

fn = transpose(fieldnames(Par.CSI));
for fn_dummy = fn
    eval([fn_dummy{:} ' = Par.CSI.' fn_dummy{:} ';']);
end
clear fn fn_dummy


%magnitude_csi = read_csi_1_4(csi_path,0,1,1);		%eigentlich unnötig CSI einzulesen ... ich weiss alles von Parameter.m
%magnitude_csi = squeeze(magnitude_csi(:,:,:,:,1));

magnitude_mask = zeros([nFreqEnc nPhasEnc nPartEnc*nSLC]);
VoI_Read = FoV_Read;
VoI_Phase = FoV_Phase;


%% 1. CREATE MASK DIRECTLY based on known VOI size

VOI_X_MIN = floor((FoV_Read(1)-VoI_Read(1))/(2*FoV_Read(1)/nFreqEnc)+0.5);
VOI_X_MAX = ceil((FoV_Read(1)+VoI_Read(1))/(2*FoV_Read(1)/nFreqEnc)+0.5);
VOI_Y_MIN = floor((FoV_Phase(1)-VoI_Phase(1))/(2*FoV_Phase(1)/nPhasEnc)+0.5);
VOI_Y_MAX = ceil((FoV_Phase(1)+VoI_Phase(1))/(2*FoV_Phase(1)/nPhasEnc)+0.5);
VOI_Z_MIN = floor((FoV_Partition-VoI_Partition)/(2*FoV_Partition/(nPartEnc*nSLC))+0.5);
VOI_Z_MAX = ceil((FoV_Partition+VoI_Partition)/(2*FoV_Partition/(nPartEnc*nSLC))+0.5);

if (VOI_X_MIN==0 && VOI_Y_MIN==0 && VOI_Z_MIN==0)%there is sum rounding bug...lukas
   magnitude_mask(VOI_X_MIN+1:VOI_X_MAX-1,VOI_Y_MIN+1:VOI_Y_MAX-1,VOI_Z_MIN+1:VOI_Z_MAX-1) = 1; 
elseif (VOI_X_MIN==0 && VOI_Y_MIN==0)
     magnitude_mask(VOI_X_MIN+1:VOI_X_MAX-1,VOI_Y_MIN+1:VOI_Y_MAX-1,VOI_Z_MIN:VOI_Z_MAX) = 1; 
else
     magnitude_mask(VOI_X_MIN:VOI_X_MAX,VOI_Y_MIN:VOI_Y_MAX,VOI_Z_MIN:VOI_Z_MAX) = 1;
end

if(Par.Flags.InterpolateCSIResolution_flag == 1)
	magnitude_mask_BefInterpol = zeros([nFreqEnc nPhasEnc nPartEnc*nSLC]);
	VOI_X_MIN_BefInterpol = floor((FoV_Read(1)-VoI_Read(1))/(2*FoV_Read(1)/nFreqEnc_BefInterpol)+0.5);
	VOI_X_MAX_BefInterpol = ceil((FoV_Read(1)+VoI_Read(1))/(2*FoV_Read(1)/nFreqEnc_BefInterpol)+0.5);
	VOI_Y_MIN_BefInterpol = floor((FoV_Phase(1)-VoI_Phase(1))/(2*FoV_Phase(1)/nPhasEnc_BefInterpol)+0.5);
	VOI_Y_MAX_BefInterpol = ceil((FoV_Phase(1)+VoI_Phase(1))/(2*FoV_Phase(1)/nPhasEnc_BefInterpol)+0.5);
	VOI_Z_MIN_BefInterpol = floor((FoV_Partition-VoI_Partition)/(2*FoV_Partition/(nPartEnc_BefInterpol*nSLC_BefInterpol))+0.5);
	VOI_Z_MAX_BefInterpol = ceil((FoV_Partition+VoI_Partition)/(2*FoV_Partition/(nPartEnc_BefInterpol*nSLC_BefInterpol))+0.5);
	
	magnitude_mask_BefInterpol(VOI_X_MIN_BefInterpol:VOI_X_MAX_BefInterpol,VOI_Y_MIN_BefInterpol:VOI_Y_MAX_BefInterpol,VOI_Z_MIN_BefInterpol:VOI_Z_MAX_BefInterpol) = 1;
	magnitude_fid_BefInterpol = fopen([tmp_dir '/mask_brain_VOI_BefInterpol.raw'],'w');
	fwrite(magnitude_fid_BefInterpol,magnitude_mask_BefInterpol,'float');
	fclose(magnitude_fid_BefInterpol);	
end

%% 2. WRITE DATA AS .RAW-FILES

magnitude_fid = fopen([tmp_dir '/mask_brain_VOI.raw'],'w');
fwrite(magnitude_fid,magnitude_mask,'float');
fclose(magnitude_fid);
