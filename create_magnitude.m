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


if(Par.Flags.image_normal_flag)                                                                                   % If Imaging data was inputted.

    magnitude = zeros(Par.CSI.total_channel_no,Par.Image.nFreqEnc,Par.Image.nPhasEnc,Par.Image.nPartEnc*Par.Image.nSLC*numel(Par.Paths.image_normal_path) );
    magnitude_flip = zeros(Par.CSI.total_channel_no,Par.Image.nFreqEnc,Par.Image.nPhasEnc,Par.Image.nPartEnc*Par.Image.nSLC*numel(Par.Paths.image_normal_path) );
    for Image_file_no = 1:numel(Par.Paths.image_normal_path)                                                      % If for each slice a different file is inputted

        magnitude_dummy = read_image(Par.Paths.image_normal_path{Image_file_no},[Par.Image.nFreqEnc Par.Image.nPhasEnc]);           % Read in data
        %magnitude_dummy = reorder_multislice_image_1_0(magnitude_dummy,4);                             % For interleaved multislice data         
        if(numel(Par.Paths.image_normal_path) > 1 && size(magnitude_dummy,4) > 1)                                 % If for each slice a different file containint ALL slices was inputted
            magnitude(:,:,:,Image_file_no) = magnitude_dummy(:,:,:,Image_file_no);                      % then cut out from the n'th file the n'th slice
        elseif(numel(Par.Paths.image_normal_path) > 1 && size(magnitude_dummy,4) == 1)                            % If for each slice a single-slice file was inputted
            magnitude(:,:,:,Image_file_no) = magnitude_dummy;                                           % then read this single slice of each file
        elseif(numel(Par.Paths.image_normal_path) == 1)                                                           % If 1 multi-slice file was inputted
            magnitude = magnitude_dummy;                                                                % then just use these data.
        end

        if(Par.Flags.image_flip_flag)                                                                             % Same for flip
            magnitude_flip_dummy = read_image(Par.Paths.image_flip_path{Image_file_no},[Par.Image.nFreqEnc,Par.Image.nPhasEnc],'ZeroFilling',1);
            %magnitude_flip_dummy = reorder_multislice_image_1_0(magnitude_flip_dummy,4);                           
            if(numel(Par.Paths.image_flip_path) > 1 && size(magnitude_flip_dummy,4) > 1)                                                                      
                magnitude_flip(:,:,:,Image_file_no) = magnitude_flip_dummy(:,:,:,Image_file_no);
            elseif(numel(Par.Paths.image_flip_path) > 1 && size(magnitude_dummy,4) == 1)
                magnitude_flip(:,:,:,Image_file_no) = magnitude_flip_dummy;
            elseif(numel(Par.Paths.image_flip_path) == 1)
                magnitude_flip = magnitude_flip_dummy;                                  
            end
        end

    end
    
    if(Par.Flags.image_flip_flag)                                                                                 % Don't put this if into the above if, because the above if is inside a loop.
        magnitude = abs((magnitude + magnitude_flip)/2);
    end
    clear magnitude_flp magnitude_dummy magnitude_flip_dummy


    
elseif(Par.Flags.image_VC_flag)                                                                                   % If ONLY VC-data was inutted for creating the mask

    magnitude = zeros(1,Par.Image.nFreqEnc,Par.Image.nPhasEnc,Par.Image.nPartEnc*Par.Image.nSLC);
	for Image_file_no = 1:numel(Par.Paths.image_VC_path)
        magnitude_weighting_dummy = abs(read_image(Par.Paths.image_VC_path{Image_file_no},[Par.Image.nFreqEnc,Par.Image.nPhasEnc]));
        %magnitude_weighting_dummy = reorder_multislice_image_1_0(magnitude_weighting_dummy,4);        
		if(numel(Par.Paths.image_VC_path) > 1 && size(magnitude_weighting_dummy,4) > 1)                                                                                 
            magnitude(:,:,:,Image_file_no) = magnitude_weighting_dummy(:,:,:,Image_file_no);
        elseif(numel(Par.Paths.image_VC_path) > 1 && size(magnitude_weighting_dummy,4) == 1)
            magnitude(:,:,:,Image_file_no) = magnitude_weighting_dummy;
        elseif(numel(Par.Paths.image_VC_path) == 1)
            magnitude = magnitude_weighting_dummy;               
		end                      
	end

else
	DataSets = 0;
	if(numel(strfind(Par.Paths.csi_path{1}, '.dat')) > 0)                                                                % If the CSI or PRESCAN of CSI data has to be used for creating the mask.
		DataSets = Analyze_mdh(Par.Paths.csi_path{1},1);
	end	
	if(isfield(DataSets,'PATREFANDIMASCAN'))
		iSpace = read_csi( Par.Paths.csi_path{1},0,'PATREFANDIMASCAN' );		
		magnitude = abs(sum(iSpace.PATREFANDIMASCAN{1},7)/2);
	else
		iSpace = read_csi( Par.Paths.csi_path{1},struct('ONLINE',struct('Hamming_flag',true)) );
		magnitude = abs(iSpace.ONLINE{1});
	end
	clear iSpace;
end 
    
    

%% 2. SOS OF ALL THE CHANNELS

magnitude = squeeze(sqrt(sum(magnitude.^2,1)));





%% 3. WRITE DATA AS .RAW-FILES

magnitude_fid = fopen([tmp_dir '/magnitude.raw'],'w');
fwrite(magnitude_fid,magnitude,'float');
fclose(magnitude_fid);






