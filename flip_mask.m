%% -1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%                           PROGRAM TO HELL AND BACK                         %%%%%%%%%%%%%%
%%%%%%%%%%%%%% ENTSTANDEN DURCH DIE KUNST DES PROGRAMMIERENS DURCH KONSEQUENTES ANSTARREN %%%%%%%%%%%%%%
%%%%%%%%%%%%%%                       READ IN MASK DATA AND FLIP IT                        %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function flip_mask(tmp_dir)
%% 0. DEFINITIONS, PREPARATIONS

load([tmp_dir '/Parameters.mat'])



%% 1. READ IN FILE

Files = dir(tmp_dir);
Files = {Files.name};
Files = strcat([tmp_dir '/'],Files);

for ii = 1:numel(Files)

    if(isempty(regexpi(Files{ii},'mask_brain.*\.raw')))
        continue
    end

    fid = fopen(Files{ii},'r');
    CurMask = fread(fid, 'float');
    fclose(fid);

    if(numel(CurMask) == Par.CSI.nFreqEnc*Par.CSI.nPhasEnc*Par.CSI.nPartEnc*Par.CSI.nSLC)
        CurMask = reshape(CurMask,[Par.CSI.nFreqEnc Par.CSI.nPhasEnc Par.CSI.nPartEnc*Par.CSI.nSLC]);
    elseif(numel(CurMask) == Par.CSI.nFreqEnc_BefInterpol*Par.CSI.nPhasEnc_BefInterpol*Par.CSI.nPartEnc_BefInterpol*Par.CSI.nSLC_BefInterpol)
        CurMask = reshape(CurMask,[Par.CSI.nFreqEnc_BefInterpol Par.CSI.nPhasEnc_BefInterpol Par.CSI.nPartEnc_BefInterpol*Par.CSI.nSLC_BefInterpol]);
    elseif(numel(CurMask) == Par.CSI.nFreqEnc*Par.CSI.nPhasEnc*Par.CSI.nPartEnc*Par.CSI.nSLC*(Par.Settings.ZeroFillMetMaps^2))
        CurMask = reshape(CurMask,[Par.CSI.nFreqEnc*Par.Settings.ZeroFillMetMaps Par.CSI.nPhasEnc*Par.Settings.ZeroFillMetMaps Par.CSI.nPartEnc*Par.CSI.nSLC]);
    elseif(numel(CurMask) == Par.CSI.nFreqEnc*Par.CSI.nPhasEnc*Par.CSI.nPartEnc*Par.CSI.nSLC*(Par.Settings.ZeroFillMetMaps^3))
        CurMask = reshape(CurMask,[Par.CSI.nFreqEnc Par.CSI.nPhasEnc Par.CSI.nPartEnc*Par.CSI.nSLC]*Par.Settings.ZeroFillMetMaps);
end

    for Slice_no = 1:size(CurMask,3)
        CurMask(:,:,Slice_no) = rot90(CurMask(:,:,Slice_no),-2);
    end    
    mask_fid = fopen(Files{ii},'w+');
    fwrite(mask_fid,CurMask,'float');
    fclose(mask_fid);

end

end
