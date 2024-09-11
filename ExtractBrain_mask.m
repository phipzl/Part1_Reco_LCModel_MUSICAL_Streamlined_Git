%% -1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%                           PROGRAM TO HELL AND BACK                         %%%%%%%%%%%%%%
%%%%%%%%%%%%%% ENTSTANDEN DURCH DIE KUNST DES PROGRAMMIERENS DURCH KONSEQUENTES ANSTARREN %%%%%%%%%%%%%%
%%%%%%%%%%%%%%                       READ IN MASK DATA AND FLIP IT                        %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ExtractBrain_mask(tmp_dir)
%% 0. DEFINITIONS, PREPARATIONS

load([tmp_dir '/Parameters.mat'])





FileList = {[tmp_dir '/mask_brain.raw']};
Sizes = {[Par.CSI.nFreqEnc Par.CSI.nPhasEnc Par.CSI.nPartEnc*Par.CSI.nSLC]};
FileList{2} = [tmp_dir '/mask_brain_zf.raw'];
if(Par.CSI.ThreeD_flag)
    Sizes{2} = [Par.CSI.nFreqEnc Par.CSI.nPhasEnc Par.CSI.nPartEnc*Par.CSI.nSLC] * Par.Settings.ZeroFillMetMaps;
else
    Sizes{2} = [Par.CSI.nFreqEnc*Par.Settings.ZeroFillMetMaps Par.CSI.nPhasEnc*Par.Settings.ZeroFillMetMaps Par.CSI.nPartEnc*Par.CSI.nSLC];
end
if(isfield(Par.CSI,'nPhasEnc_BefInterpol'))
    FileList{3} = [tmp_dir '/mask_brain_BefInterpol.raw'];
    Sizes{3} = [Par.CSI.nFreqEnc_BefInterpol Par.CSI.nPhasEnc_BefInterpol Par.CSI.nPartEnc_BefInterpol*Par.CSI.nSLC_BefInterpol];
end

for ii = 1:numel(Sizes)
	CurFile2 = FileList{ii};
	if(~exist(CurFile2,'file'))
		continue
	end
	
	%% 1. READ IN FILE
	fid = fopen(CurFile2,'r');
	% mask = zeros(ROW,COL,SLC);
	% for Slice_no = 1:SLC
	%     mask(:,:,Slice_no) = fread(fid, [ROW,COL], 'float');
	% end
	mask = fread(fid, 'float');
        mask = reshape(mask, Sizes{ii});
	fclose(fid);


	%% 2. Extract Largest Contiguous Mask Region

	mask2 = zeros(size(mask));
	for slc = 1:size(mask,3)
		dummy = bwconncomp(mask(:,:,slc),4);
		dummy = dummy.PixelIdxList(max(cellfun(@numel,dummy.PixelIdxList)) == cellfun(@numel,dummy.PixelIdxList));
		mask2_dum = mask2(:,:,slc);
        mask2_dum([dummy{:}]) = 1;
        if(any(mask2_dum(:))>0)
            mask2_dum = imfill(mask2_dum);          % Fill in all holes
		mask2(:,:,slc) = mask2_dum;
        end
	end

	%% 3. WRITE DATA AS .RAW-FILES

	mask_fid = fopen(CurFile2,'w+');
	% for Slice_no = 1:SLC
	%     fwrite(mask_fid,mask(:,:,Slice_no),'float');
	% end
	fwrite(mask_fid,mask2,'float');
	fclose(mask_fid);

end

end
