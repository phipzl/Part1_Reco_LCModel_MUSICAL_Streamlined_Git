%% -1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%                           PROGRAM TO HELL AND BACK                         %%%%%%%%%%%%%%
%%%%%%%%%%%%%% ENTSTANDEN DURCH DIE KUNST DES PROGRAMMIERENS DURCH KONSEQUENTES ANSTARREN %%%%%%%%%%%%%%
%%%%%%%%%%%%%%                       READ IN MASK DATA AND FLIP IT                        %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 0. DEFINITIONS, PREPARATIONS

clearvars -except tmp_dir; close all;
load([tmp_dir '/Parameters.mat'])



%% 1. READ IN FILE

fid = fopen([tmp_dir '/mask_brain.raw'],'r');
% mask = zeros(ROW,COL,SLC);
% for Slice_no = 1:SLC
%     mask(:,:,Slice_no) = fread(fid, [ROW,COL], 'float');
% end
mask = reshape(fread(fid, 'float'), [Par.CSI.nFreqEnc Par.CSI.nPhasEnc Par.CSI.nPartEnc*Par.CSI.nSLC]);

fclose(fid);


%% 2. FLIP THE MASK DATA

for Slice_no = 1:Par.CSI.nPartEnc*Par.CSI.nSLC
    mask(:,:,Slice_no) = rot90(mask(:,:,Slice_no),-2);
end


%% 3. WRITE DATA AS .RAW-FILES

mask_fid = fopen([tmp_dir '/mask_brain.raw'],'w+');
% for Slice_no = 1:SLC
%     fwrite(mask_fid,mask(:,:,Slice_no),'float');
% end
fwrite(mask_fid,mask,'float');
fclose(mask_fid);



if(Par.Flags.InterpolateCSIResolution_flag == 1)
	fid = fopen([tmp_dir '/mask_brain_BefInterpol.raw'],'r');
	mask_BefInterpol = reshape(fread(fid, 'float'), [Par.CSI.nFreqEnc_BefInterpol Par.CSI.nPhasEnc_BefInterpol Par.CSI.nPartEnc_BefInterpol*Par.CSI.nSLC_BefInterpol]);
	fclose(fid);
	for Slice_no = 1:Par.CSI.nPartEnc_BefInterpol*Par.CSI.nSLC_BefInterpol
		mask_BefInterpol(:,:,Slice_no) = rot90(mask_BefInterpol(:,:,Slice_no),-2);
	end
	mask_fid = fopen([tmp_dir '/mask_brain_BefInterpol.raw'],'w+');
	fwrite(mask_fid,mask_BefInterpol,'float');
	fclose(mask_fid);
end