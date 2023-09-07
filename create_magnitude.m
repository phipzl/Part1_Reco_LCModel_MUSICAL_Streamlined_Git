%% -1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%                           PROGRAM TO HELL AND BACK                             %%%%%%%%%%%%%%
%%%%%%%%%%%%%% ENTSTANDEN DURCH DIE KUNST DES PROGRAMMIERENS DURCH KONSEQUENTES ANSTARREN     %%%%%%%%%%%%%%
%%%%%%%%%%%%%%   READ IN THE IMAGING OR CSI DATA, SUM UP ALL CHANNELS, WRITE OUTPUT AS .RAW   %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 0. DEFINITIONS, PREPARATIONS

close all;
load([tmp_dir '/Parameters.mat'])

% fn = transpose(fieldnames(Par.CSI));
% for fn_dummy = fn
%     eval([fn_dummy{:} ' = Par.CSI.' fn_dummy{:}]);
% end
% clear fn fn_dummy







%% 1. READ IN FILE

if(isfield(Par.Paths,'NonCartTrajFile_path'))
    NonCartTrajFile_path = Par.Paths.NonCartTrajFile_path{1};
else
    NonCartTrajFile_path = [];
end

if(Par.Flags.image_normal_flag)                                                                                   % If Imaging data was inputted.

    Sett.NonCartReco.Phaseroll_flag = false;
    magnitude = op_ReadAndRecoSiemensData(Par.Paths.image_normal_path{1},NonCartTrajFile_path,[],Sett);           % Read in data
    magnitude = magnitude.Data(:,:,:,1,:);
    if(Par.Flags.image_flip_flag)                                                                             % Same for flip
        magnitude_flip = op_ReadAndRecoSiemensData(Par.Paths.image_flip_path{1},NonCartTrajFile_path,[],Sett);
        magnitude_flip = magnitude_flip.Data(:,:,:,1,:);        
    end
    
    if(Par.Flags.image_flip_flag)                                                                                 % Don't put this if into the above if, because the above if is inside a loop.
        magnitude = ((magnitude + magnitude_flip)/2);
    end
    clear magnitude_flp

elseif(Par.Flags.image_VC_flag)                                                                                   % If ONLY VC-data was inutted for creating the mask
        Sett.NonCartReco.Phaseroll_flag = false;
        magnitude_weighting = (read_image(Par.Paths.image_VC_path{1},[],Sett));

else
    Sett.NonCartReco.Phaseroll_flag = false;
    [MRSI,RefScan] = op_ReadAndRecoSiemensData(Par.Paths.csi_path{1},Par.Paths.NonCartTrajFile_path{1},[],Sett);
    if(~isempty(fieldnames(RefScan)))
        magnitude = RefScan.Data(:,:,:,1,:,:,:);
	else
        magnitude = MRSI.Data(:,:,:,1,:,:,:);        
	end
    clear MRSI RefScan
end 
    
    

%% 2. SOS OF ALL THE CHANNELS

magnitude = HammingFilter(magnitude,[1 2 3],100,'OuterProduct',0);
magnitude = HammingFilter(magnitude,[1 2 3],100,'OuterProduct',0);
magnitude = squeeze(sqrt(sum(abs(magnitude).^2,5)));





%% 3. WRITE DATA AS .RAW-FILES

magnitude_fid = fopen([tmp_dir '/magnitude.raw'],'w');
fwrite(magnitude_fid,magnitude,'float');
fclose(magnitude_fid);






