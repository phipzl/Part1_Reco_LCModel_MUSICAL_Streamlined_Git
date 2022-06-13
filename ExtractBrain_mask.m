%% -1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%                           PROGRAM TO HELL AND BACK                         %%%%%%%%%%%%%%
%%%%%%%%%%%%%% ENTSTANDEN DURCH DIE KUNST DES PROGRAMMIERENS DURCH KONSEQUENTES ANSTARREN %%%%%%%%%%%%%%
%%%%%%%%%%%%%%                       READ IN MASK DATA AND FLIP IT                        %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 0. DEFINITIONS, PREPARATIONS

clearvars -except tmp_dir; close all;
load([tmp_dir '/Parameters.mat'])






for CurFile = {[tmp_dir '/mask_brain.raw'],[tmp_dir '/mask_brain_BefInterpol.raw']}
	CurFile2 = Curfile{:};
	if(~exist(CurFile2,'file'))
		continue
	end
	
	%% 1. READ IN FILE
	fid = fopen(CurFile2,'r');
	% mask = zeros(ROW,COL,SLC);
	% for Slice_no = 1:SLC
	%     mask(:,:,Slice_no) = fread(fid, [ROW,COL], 'float');
	% end
	mask = reshape(fread(fid, 'float'), [Par.CSI.nFreqEnc Par.CSI.nPhasEnc Par.CSI.nPartEnc*Par.CSI.nSLC]);

	fclose(fid);


	%% 2. Extract Largest Contiguous Mask Region

	mask2 = zeros(size(mask));
	for slc = 1:size(mask,3)
		dummy = bwconncomp(mask(:,:,slc),4);
		dummy = dummy.PixelIdxList(max(cellfun(@numel,dummy.PixelIdxList)) == cellfun(@numel,dummy.PixelIdxList));
		mask2_dum = mask2(:,:,slc);
		mask2_dum(dummy{:}) = 1;
		mask2(:,:,slc) = mask2_dum;
	end

	%% 3. WRITE DATA AS .RAW-FILES

	mask_fid = fopen(CurFile2,'w+');
	% for Slice_no = 1:SLC
	%     fwrite(mask_fid,mask(:,:,Slice_no),'float');
	% end
	fwrite(mask_fid,mask2,'float');
	fclose(mask_fid);

end

