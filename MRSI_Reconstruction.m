%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%                           PROGRAM TO HELL AND BACK                         %%%%%%%%%%%%%%
%%%%%%%%%%%%%% ENTSTANDEN DURCH DIE KUNST DES PROGRAMMIERENS DURCH KONSEQUENTES ANSTARREN %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -except tmp_dir


if(~exist('CurAvg','var'))		% For those cases someone wants to run this script manually
    CurAvg = 1;
end

CoilCombtic = tic;
precision = 'single';

%pause on

%% Load Parameters

if (exist([tmp_dir '/Parameters_water.mat'],'file'))
    load([tmp_dir '/Parameters_water.mat'])
    Par.Paths.csi_path{CurAvg} = Par.Paths.WaterReference_path;
else
    load([tmp_dir '/Parameters.mat'])
end
if(~isfield(Par.ServerInfo,'RunLCModel_CPUCores'))
	Par.ServerInfo.RunLCModel_CPUCores = 8;
end

FastPIReprocess_flag = false;

%% Read Noise Data

if(~FastPIReprocess_flag)
	if(Par.Flags.noisedecorrelation_flag && isfield(Par.Paths, 'noisedecorrelation_path'))
		fprintf('\nPreparing for Noise Decorrelation.')
        NoiseData = io_ReadAndReshapeSiemensData(Par.Paths.noisedecorrelation_path);
	end	
end


%% Read in CSI file

if(~FastPIReprocess_flag)
    fprintf('\n\nREAD CSI DATA')  
	if(~Par.Flags.NonCartTraj_flag)
        [csi, image, NoiseData] = io_ReadAndReshapeSiemensData(Par.Paths.csi_path(CurAvg));
	else
        [csi,image, NoiseData] = io_ReadAndReshapeSiemensData(Par.Paths.csi_path(CurAvg),Par.Paths.NonCartTrajFile_path{1});     % Spiral sequences need external trajectory file
	end
    
    if(strcmpi(csi.Par.AssumedSequence, 'CSIOrSVS'))
        if(isfield(image,'Data'))
            image.Data = sum(image.Data,7)/2; % Average normal & flip images
        end     
    else
         Par.CSI.vecSize = csi.Par.vecSize;          
    end

end

% Calc NoiseCorrMat
if(isfield(NoiseData,'Data') && numel(NoiseData.Data) > 1 && Par.Flags.noisedecorrelation_flag)
    [NoiseCorrMatStruct,Dummy] = op_CalcNoiseCorrMat(NoiseData);
    NoiseData = Dummy.NoiseData; clear Dummy;
end
if(~isfield(image,'Data'))
    clear image;
end


%% Plot Noise Data

if(~FastPIReprocess_flag)
	
	if(exist('NoiseCorrMatStruct','var'))

		mkdir([Par.Paths.out_path '/noise'])

		NoisecorrMat_fig = figure('visible','off');   
		imagesc(abs(NoiseCorrMatStruct.Data))
		colorbar;
		saveas(NoisecorrMat_fig,sprintf('%s/noise/NoiseCorrMat', Par.Paths.out_path),'epsc2')
		close(NoisecorrMat_fig)


		SubFigRow = ceil(sqrt(Par.CSI.total_channel_no)); SubFigCol = ceil(Par.CSI.total_channel_no / SubFigRow);
		Noise_fig = figure('visible','off');   
		set(gcf, 'Position', get(0,'Screensize')); % Maximize figure

		for cha_no = 1: Par.CSI.total_channel_no
			subplot(SubFigRow,SubFigCol,cha_no);
			plot(real(squeeze(NoiseData.Data(cha_no,:))));
			axis off
			if(cha_no == floor(SubFigCol/2))
				title('Noise Used for Noise Decorrelation','Interpreter','none')
         end
         end

		saveas(Noise_fig,sprintf('%s/noise/Noise', Par.Paths.out_path),'epsc2') 
		close(Noise_fig)
		clear Noise_mat

     end

 end


%% READ IN IMAGING-DATA, IMAGING-FLIP-DATA, IMAGING-VC-DATA, MASK-DATA, PREPARE DATA

if(~FastPIReprocess_flag)
    fprintf('\n\nREAD IMAGING DATA AND MASK')  

	% IMAGE NORMAL
	if(Par.Flags.image_normal_flag)
        if(Par.Flags.NonCartTraj_flag)
            image = io_ReadAndReshapeSiemensData(Par.Paths.image_normal_path{1},Par.Paths.NonCartTrajFile_path{1});
        else
            image = io_ReadAndReshapeSiemensData(Par.Paths.image_normal_path{1}); 
        end
	end


	% IMAGE FLIP
	if(Par.Flags.image_flip_flag)
        if(Par.Flags.NonCartTraj_flag)
            image_flip = io_ReadAndReshapeSiemensData(Par.Paths.image_flip_path{1},Par.Paths.NonCartTrajFile_path{1});
        else
            image_flip = io_ReadAndReshapeSiemensData(Par.Paths.image_flip_path{1}); 
        end
    end    
    
    
	% IMAGE VC
	if(Par.Flags.image_VC_flag)
        if(Par.Flags.NonCartTraj_flag)
            image_VC = io_ReadAndReshapeSiemensData(Par.Paths.image_VC_path{1},Par.Paths.NonCartTrajFile_path{1});
        else
            image_VC = io_ReadAndReshapeSiemensData(Par.Paths.image_VC_path{1}); 
        end
	end    
    


%% MASK
	if(exist([tmp_dir '/mask_brain.raw'],'file'))
		fid_mask = fopen([tmp_dir '/mask_brain.raw'],'r');
		mask = reshape(fread(fid_mask, 'float'), [Par.CSI.nFreqEnc Par.CSI.nPhasEnc Par.CSI.nSLC*Par.CSI.nPartEnc]);
		fclose(fid_mask);
		
		if(Par.Flags.InterpolateCSIResolution_flag == 1)
			fid_mask = fopen([tmp_dir '/mask_brain_BefInterpol.raw'],'r');
			mask_BefInterpol = reshape(fread(fid_mask, 'float'), [Par.CSI.nFreqEnc_BefInterpol Par.CSI.nPhasEnc_BefInterpol Par.CSI.nSLC_BefInterpol*Par.CSI.nPartEnc_BefInterpol]);
			fclose(fid_mask);
		else
			mask_BefInterpol = mask;
		end
		
		if(sum(sum(sum(mask))) == 0)
			mask = ones([Par.CSI.nFreqEnc Par.CSI.nPhasEnc Par.CSI.nPartEnc*Par.CSI.nSLC]);
		end
	else
		mask = ones([Par.CSI.nFreqEnc Par.CSI.nPhasEnc Par.CSI.nPartEnc*Par.CSI.nSLC]);
	end



    %% Noise-Decorrelate Data
    if(exist('NoiseCorrMatStruct','var'))
        if(exist('NoiseData','var'))
            Settings.CreateSyntheticNoise=1;
        end
        csi = op_PerformNoiseDecorrelation(csi,NoiseCorrMatStruct,Settings);
        if(exist('image','var'))
            Settings.CreateSyntheticNoise = 0;
            image = op_PerformNoiseDecorrelation(image,NoiseCorrMatStruct,Settings);
        end
    end


end
    
%% Averaging

if(exist('image','var'))
    image = op_AverageMRData(image);        % Currently does nothing
end
csi = op_AverageMRData(csi);        % Currently does nothing



%% DEBUG MODE: Compute Std and NoiseCorrelation Matrix of Decorrelated Noise and Decorrelated CSI data

% 
% 
% DecorrelationMatrix2 = inv(chol(Noise_CorrMat * 0.5,'lower'));
% Noise_mat_Decorr = DecorrelationMatrix2 * Noise_mat;
% 
% 
% 
% % Gather points at end of FID
% nFIDendpoints = 400;
% TakeNPointsOutOfEnd = 200;
% randy = randperm(nFIDendpoints); randy = randy(1:TakeNPointsOutOfEnd); % Take 'em randomly
% % Take random points at end of FID
% Noise_csi = csi_k(:,:,:,:,end - (nFIDendpoints - 1) : end); Noise_csi = Noise_csi(:,:,:,:,randy);
% 
% 
% % Take only certain k-space points
% [Elliptical_dummy,Elliptical_mask] = EllipticalFilter(squeeze(csi_k(1,:,:,1,1)),[1 2],[1 1 1 size(csi_k,3)/2-1],1);
% [Elliptical_dummy,Elliptical_mask2] = EllipticalFilter(squeeze(csi_k(1,:,:,1,1)),[1 2],[1 1 1 size(csi_k,3)/2-2],1);
% 
% PI_mask = abs(squeeze(csi_k(1,:,:,1,1))); PI_mask(PI_mask > 0) = 1;
% csi_mask = (Elliptical_mask - Elliptical_mask2) .* PI_mask;
% csi_mask = repmat(logical(csi_mask), [1 1 size(csi_k,1)*size(csi_k,4)*TakeNPointsOutOfEnd]);
% csi_mask = reshape(csi_mask, [size(csi_k,2) size(csi_k,3) size(csi_k,1) size(csi_k,4) TakeNPointsOutOfEnd]);
% csi_mask = permute(csi_mask, [3 1 2 4 5]);
% %csi_mask = reshape(repmat(csi_mask, [size(csi_k,1) size(csi_k,4)*TakeNPointsOutOfEnd]), [size(csi_k,1) size(csi_k,2) size(csi_k,3) size(csi_k,4) TakeNPointsOutOfEnd]);
% % Take only csi-points which are farest away from k-space center (so a circle with radius 31 k-space points)
% Noise_mat_csi = Noise_csi(csi_mask);
% Noise_mat_csi = reshape(Noise_mat_csi, [32 numel(Noise_mat_csi)/32]);
% 
% 
% 
% 
% % Check if Noise itself was decorrelated
% CorrMat_OfDecorrelatedNoise = 1/(size(Noise_mat,2)) * (Noise_mat_Decorr * Noise_mat_Decorr');
% 
% % Compute NoiseDecorrelation of CSI data
% CorrMat_OfDecorrelatedCsiNoise = 1/(size(Noise_mat_csi,2)) * (Noise_mat_csi * Noise_mat_csi');
% 
% 
% % Figures
% figure; imagesc(abs(CorrMat_OfDecorrelatedNoise))
% figure; imagesc(abs(CorrMat_OfDecorrelatedCsiNoise))
% 
% figure; plot(squeeze(real(Noise_mat(1,:))))
% figure; plot(squeeze(real(Noise_mat_Decorr(1,:))))
% figure; plot(squeeze(real(Noise_mat_csi(1,:))))
% 
% 
% % Std
% Std_DecorrelatedNoise = cat(2,std(real(Noise_mat_Decorr),0,2), std(imag(Noise_mat_Decorr),0,2), std(real(Noise_mat_csi),0,2), std(imag(Noise_mat_csi),0,2))
% mean(Std_DecorrelatedNoise)



%% Read B0Map

if(Par.Flags.AlignFreq_flag && exist(['./' tmp_dir '/AlignFreq_B0FieldMap.raw'],'file'))
    AlignFreq_B0Fieldmap = read_RawFiles(['./' tmp_dir '/AlignFreq_B0FieldMap.raw'],Par.AlignFreq.nFreqEnc,Par.AlignFreq.nPhasEnc,Par.AlignFreq.nPartEnc,'float32');
    AlignFreq_B0Fieldmap( isnan(AlignFreq_B0Fieldmap) | isinf(AlignFreq_B0Fieldmap) ) = 0;
end



%% Perform Parallel Imaging Simulation and Reconstruction

if(Par.Flags.TwoDCaipParallelImaging_flag || Par.Flags.SliceParallelImaging_flag)
    fprintf('\n\nSIMULATE AND RECONSTRUCT PARALLEL IMAGING')  
    run ./ParallelImagingSimReco.m
end




%% Perform Fourier Transform of image

if(~FastPIReprocess_flag && exist('image','var') && ~(isfield(image.Par,'dicom_flag') && image.Par.dicom_flag))
    if(exist('image','var'))
        if(~isfield(image.Par,'SpatialSpectralEncoding_flag'))
            image.Par.SpatialSpectralEncoding_flag = false;
        end
        if(image.Par.SpatialSpectralEncoding_flag) 
    
            Settings.NonCartReco.DensComp.Method = 'ConcentricRingTrajectory_Theoretical';   
            Settings.NonCartReco.DensComp.ApplyHammingFilter_flag = Par.Flags.hamming_flag;
            Settings.NonCartReco.ConjInkSpace_flag = true; 
            Settings.NonCartReco.Phaseroll_flag = true;
            Settings.NonCartReco.ConjIniSpace_flag = false;    
            Settings.NonCartReco.FlipDim = 1;
            Settings.NonCartReco.CircularSFTFoV_flag = false; 
            Settings.PreWhitenData_flag = 0;
            if(isfield(Par.Settings,'GradientDelayPerAngInt_x') && isfield(Par.Settings,'GradientDelayPerAngInt_y'))
                Settings.ReadInTraj.GradDelayPerAngInt_x_us = Par.Settings.GradientDelayPerAngInt_x;
                Settings.ReadInTraj.GradDelayPerAngInt_y_us = Par.Settings.GradientDelayPerAngInt_y;
            
            elseif(isfield(Par.Settings,'GradientDelayPerTempInt_x') && isfield(Par.Settings,'GradientDelayPerTempInt_y'))
                Settings.ReadInTraj.GradDelayPerTempInt_x_us = Par.Settings.GradientDelayPerTempInt_x;
                Settings.ReadInTraj.GradDelayPerTempInt_y_us = Par.Settings.GradientDelayPerTempInt_y;
            else
                Settings.ReadInTraj.GradDelayPerTempInt_x_us = 0;
                Settings.ReadInTraj.GradDelayPerTempInt_y_us = 0;
            end            
            
            image = op_ReconstructMRData(image,struct(),Settings);
          
        else
          
            % Do z-fft before conj that is sometimes done in in-plane FFT
            if(image.RecoPar.nPartEnc > 1)
                Settings.zCartFFT.ConjFlag = false;    
                Settings.zCartFFT.Ifft_flag = false;    
                Settings.zCartFFT.FlipDim_flag = false;
                Settings.zCartFFT.ApplyAlongDims = [3];
                image = op_FFTOfMRIData_v2(image,Settings.zCartFFT);        
            end
            
            Settings.CartFFT.ConjFlag = true;    
            Settings.CartFFT.Ifft_flag = false;    
            Settings.CartFFT.FlipDim_flag = true;
            Settings.CartFFT.FlipDim = 1;
            Settings.CartFFT.ApplyAlongDims = [1 2];

            image = op_FFTOfMRIData_v2(image,Settings.CartFFT);
            
            image = op_SliceReco(image);


        end

    %         if(Par.Flags.ESPIRiT_flag)      % Directly perform the coil combination here. Maybe later I can just calculate the data for the coil combination...
    %             image = op_CalcSensMaps(image,[],image.RecoPar,struct('Extrapolate_flag',true,'RatioToMax4MaskThreshold',0.18));
    %             image = op_PermuteMRData(image,[5 1 2 3 4]);
    % 
   % end
    end

    if(exist('image_vc','var'))
        if(~isfield(image_vc.Par,'SpatialSpectralEncoding_flag'))
            image_vc.Par.SpatialSpectralEncoding_flag = false;
        end        
        if(image_vc.Par.SpatialSpectralEncoding_flag) 
            image_vc = op_ReconstructMRData(image_vc,struct(),Settings);
        else
            image_vc = op_FFTOfMRIData_v2(image_vc,Settings.CartFFT);
        end
        image_vc = op_SliceReco(image_vc);
        end

    if(exist('image_flip','var'))
        if(~isfield(image_flip.Par,'SpatialSpectralEncoding_flag'))
            image_flip.Par.SpatialSpectralEncoding_flag = false;
        end  
        if(image_flip.Par.SpatialSpectralEncoding_flag) 
            Settings.NonCartReco.FlipDim = 2;
            image_flip = op_ReconstructMRData(image_flip,struct(),Settings);
        else
            Settings.CartFFT.FlipDim = 2;
            image_flip = op_FFTOfMRIData_v2(image_flip,Settings.CartFFT);
        end
        image_flip = op_SliceReco(image_flip);
        image_flip.Data = circshift(image_flip.Data, [-1 0 0 0]);
        image.Data = ( image.Data + image_flip.Data )/2;
    end

end
    
if(size(image.Data,4) > 3)
    image.Data = image.Data(:,:,:,4,:);
else
    image.Data = image.Data(:,:,:,1,:);
end

    


%% Perform Fourier Transform of csi

if(~FastPIReprocess_flag && ~(isfield(csi.Par,'dicom_flag') && csi.Par.dicom_flag))
    if(~isfield(csi.Par,'SpatialSpectralEncoding_flag'))
        csi.Par.SpatialSpectralEncoding_flag = false;
    end 
    if(csi.Par.SpatialSpectralEncoding_flag)
        
        Settings.NonCartReco.DensComp.Method = 'ConcentricRingTrajectory_Theoretical';  
        Settings.NonCartReco.DensComp.ApplyHammingFilter_flag = Par.Flags.hamming_flag;
        Settings.NonCartReco.ConjInkSpace_flag = true; 
        Settings.NonCartReco.CircularSFTFoV_flag = false; 
        Settings.NonCartReco.FlipDim_flag = true;
        Settings.NonCartReco.Phaseroll_flag = true;
        Settings.NonCartReco.ConjIniSpace_flag = false; 
        Settings.PreWhitenData_flag = 0;
        
        if(isfield(Par.Settings,'GradientDelayPerAngInt_x') && isfield(Par.Settings,'GradientDelayPerAngInt_y'))
            Settings.ReadInTraj.GradDelayPerAngInt_x_us = Par.Settings.GradientDelayPerAngInt_x;
            Settings.ReadInTraj.GradDelayPerAngInt_y_us = Par.Settings.GradientDelayPerAngInt_y;

        elseif(isfield(Par.Settings,'GradientDelayPerTempInt_x') && isfield(Par.Settings,'GradientDelayPerTempInt_y'))
            Settings.ReadInTraj.GradDelayPerTempInt_x_us = Par.Settings.GradientDelayPerTempInt_x;
            Settings.ReadInTraj.GradDelayPerTempInt_y_us = Par.Settings.GradientDelayPerTempInt_y;
        else
            Settings.ReadInTraj.GradDelayPerTempInt_x_us = 0;
            Settings.ReadInTraj.GradDelayPerTempInt_y_us = 0;
        end  

        % Reco MRSI Data
        csi = op_ReconstructMRData(csi,struct(),Settings);
        
        
    else
        
        % Do z-fft before conj that is sometimes done in in-plane FFT
        if(csi.RecoPar.nPartEnc > 1)
            Settings.zCartFFT.ConjFlag = false;    
            Settings.zCartFFT.Ifft_flag = false;    
            Settings.zCartFFT.FlipDim_flag = false;
            Settings.zCartFFT.ApplyAlongDims = [3];
            csi = op_FFTOfMRIData_v2(csi,Settings.zCartFFT);        
        end

        Settings.CartFFT.ConjFlag = true;    
        Settings.CartFFT.Ifft_flag = false;    
        Settings.CartFFT.FlipDim_flag = true;
        Settings.CartFFT.FlipDim = 1;
        Settings.CartFFT.ApplyAlongDims = [1 2];
        csi = op_FFTOfMRIData_v2(csi,Settings.CartFFT);
        csi = op_SliceReco(csi);        
    end

end  


%% Apply Hamming Filter to Combined MRSI Data
size_csi = size(csi.Data);
if(Par.Flags.hamming_flag)
	if( ~isfield(csi.Par,'SpatialSpectralEncoding_flag') || ~csi.Par.SpatialSpectralEncoding_flag)      % For SpatialSpectralEncoding we have done it already during reco
    
	    fprintf('\n\nAPPLY HAMMING FILTER')  % If Full_ElliptWeighted_Or_Weighted_Acq = 4, the data is alrdy intrinsically filtered in z-dimension
	    % Was the WeightedAcquisition undone?
	    UndoWeightedAcq_flag = false;
	    if(isfield(csi,'RecoSteps'))
            fieldies = fieldnames(csi.RecoSteps); field = fieldies(~cellfun(@isempty,regexpi(fieldies,'op_AverageMRData')));
            if(~isempty(field))
                UndoWeightedAcq_flag = csi.RecoSteps.(field{1}).UndoWeightedAveraging_flag;
            end
	    end
	    if(size_csi(3) > 1 && Par.CSI.ThreeD_flag && (Par.CSI.Full_ElliptWeighted_Or_Weighted_Acq ~= 4 || UndoWeightedAcq_flag))      
            csi.Data = HammingFilter(csi.Data,[1 2 3],Par.Settings.hamming_factor,'OuterProduct',0);     
        else
            csi.Data = HammingFilter(csi.Data,[1 2],Par.Settings.hamming_factor,'OuterProduct',0);
	    end 
    elseif(isfield(csi.Par,'SpatialSpectralEncoding_flag') && csi.Par.SpatialSpectralEncoding_flag && size_csi(3) > 1 && Par.CSI.ThreeD_flag)
        fprintf('\n\nAPPLY HAMMING FILTER IN z-DIRECTION')  % If Full_ElliptWeighted_Or_Weighted_Acq = 4, the data is alrdy intrinsically filtered in z-dimension
		csi.Data = HammingFilter(csi.Data,[3],Par.Settings.hamming_factor,'OuterProduct',0);   

	end
end


%% Interpolate in kSpace


% NOT YET REIMPLEMENTED!
if(Par.Flags.InterpolateCSIResolution_flag && Par.Settings.InterpolateCSIResolution_InkSpace)
    if(exist('csi','var'))
        csi = ZerofillOrCutkSpace(csi,[size(csi,1) Par.Settings.InterpolateCSIResolution size(csi,5)],1);
        if(Par.Settings.InterpolateCSIResolution_EllipFilter)
            csi = EllipticalFilter(csi, [2 3] ,[1 1 1 Par.Settings.InterpolateCSIResolution(1)/2-1],false);      % So far only 2D and cylindrical 3D-kspaces (circle in kx-ky, everything in kz), and only to quadratic k-spaces
        end
    end
    if(exist('csi_k','var'))
        csi_k = ZerofillOrCutkSpace(csi_k,[size(csi_k,1) Par.Settings.InterpolateCSIResolution size(csi_k,5)],0);
        if(Par.Settings.InterpolateCSIResolution_EllipFilter)
            csi_k = EllipticalFilter(csi_k, [2 3] ,[1 1 1 Par.Settings.InterpolateCSIResolution(1)/2-1],true);      % So far only 2D and cylindrical 3D-kspaces (circle in kx-ky, everything in kz), and only to quadratic k-spaces
     end
end
    if(exist('image','var'))
        image = ZerofillOrCutkSpace(image,[size(image,1) Par.Settings.InterpolateCSIResolution],1);
        if(Par.Settings.InterpolateCSIResolution_EllipFilter)
            image = EllipticalFilter(image, [2 3] ,[1 1 1 Par.Settings.InterpolateCSIResolution(1)/2-1],false);      % So far only 2D and cylindrical 3D-kspaces (circle in kx-ky, everything in kz), and only to quadratic k-spaces
        end
    end
    if(exist('noise_sim','var'))
        noise_sim = ZerofillOrCutkSpace(noise_sim,[size(noise_sim,1) Par.Settings.InterpolateCSIResolution size(noise_sim,5)],0);
        if(Par.Settings.InterpolateCSIResolution_EllipFilter)
            noise_sim = EllipticalFilter(noise_sim, [2 3] ,[1 1 1 Par.Settings.InterpolateCSIResolution(1)/2-1],true);      % So far only 2D and cylindrical 3D-kspaces (circle in kx-ky, everything in kz), and only to quadratic k-spaces
        end
    end
    if(exist('noise_sim_PI','var'))
        noise_sim_PI = ZerofillOrCutkSpace(noise_sim_PI,[size(noise_sim,1) Par.Settings.InterpolateCSIResolution size(noise_sim,5)],1);
        if(Par.Settings.InterpolateCSIResolution_EllipFilter)
            noise_sim_PI = EllipticalFilter(noise_sim_PI, [2 3] ,[1 1 1 Par.Settings.InterpolateCSIResolution(1)/2-1],false);      % So far only 2D and cylindrical 3D-kspaces (circle in kx-ky, everything in kz), and only to quadratic k-spaces
        end
    end
end
size_csi = size(csi.Data);
size_csi(5) = size(csi.Data,5);


%% DEBUG MODE: WRITE PHASEMAPS OF IMAGE_CORR

% if(~FastPIReprocess_flag)
% 
% 	display([char(10) char(10) 'WRITE PHASEMAPS'])
% 
% 	for channel_no = 1:size_csi(1)
% 		for slice_no = 1:size_csi(4)
% 
% 			% CSI PHAMAP
% 			figure('visible','off')
% 			imagesc(-rad2deg(angle(squeeze(csi(channel_no,:,:,slice_no,1)))),[-180 180])
% 			title(sprintf('CSI phamap Slice %d channel %02d',slice_no, channel_no)) 
% 			colorbar
% 			%saveas(gcf,sprintf('%s/phamaps/CSI_Slice%d_channel%02d.fig', Par.Paths.out_path,slice_no,channel_no))
% 			saveas(gcf,sprintf('%s/phamaps/CSI_Slice%d_channel%02d.jpg', Par.Paths.out_path,slice_no,channel_no))
% 			close(gcf)
% 
% 
% 			if(exist('image','var'))
% 
% 
% 				% GRE PHAMAP
% 				figure('visible','off')
% 				imagesc(rad2deg(angle(squeeze(image(channel_no,:,:,slice_no)))),[-180 180])
% 				title(sprintf('GRE phamap Slice %d channel %02d',slice_no,channel_no))      
% 				colorbar
% 				%saveas(gcf,sprintf('%s/phamaps/GRE_Slice%d_channel%02d.fig', Par.Paths.out_path,slice_no,channel_no))
% 				saveas(gcf,sprintf('%s/phamaps/GRE_Slice%d_channel%02d.jpg', Par.Paths.out_path,slice_no,channel_no))
% 				close(gcf)  
% 
% 
% 				if(~Par.Flags.image_VC_flag)          % Don't print this subtraction for sensmap method, because csi and image have different sizes in that case.
% 					% SUB CSI - GRE
% 					figure('visible','off')
% 					imagesc(-rad2deg(angle(squeeze(csi(channel_no,:,:,slice_no,1)))) - rad2deg(angle(squeeze(image(channel_no,:,:,slice_no)))), [-70 70])
% 					title(sprintf('SUB CSI - GRE phamap Slice %d channel %02d',slice_no,channel_no))      
% 					colorbar
% 					%saveas(gcf,sprintf('%s/phamaps/SUB_CSI-GRE_Slice%d_channel%02d.fig', Par.Paths.out_path,slice_no,channel_no))
% 					saveas(gcf,sprintf('%s/phamaps/SUB_CSI-GRE_Slice%d_channel%02d.jpg', Par.Paths.out_path,slice_no,channel_no))
% 					close(gcf)    
% 				end
% 
% 			end
% 
% 		end
% 
% 	end
% end















%% Bilgic Lipid Decontamination
% INSERT PRE-COIL-COMB LIPID DECONTAMINATION HERE!

%% PRE COIL COMBI BILGIC LIPID SUPP



if(Par.Flags.LipidDecon_flag == 1 && ~exist([tmp_dir '/Parameters_water.mat'],'file'))

    csi_bak = csi;
    clear csi;
    csi = csi_bak.Data;
    %save([Par.Paths.out_path '/UnCombinedCSI.mat'],'csi','-v7.3');
    %temp_x=load('/ceph/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Process_Results/MS_3DCRT_Berni_pipeline/test_wLukas/UnCombinedCSI.mat');
    %temp_x=permute(temp_x.csi,[2 3 4 5 1]);
    %csi=temp_x;
    csi=csi*100; % It scales csi data similarly like with Lukas's pipeline --> it makes the L2 work with the same factor
    fprintf('\n\nPerform Bilgic Lipid Decontamination: channelwise & sensitivity-weighted\n')

    mkdir([Par.Paths.out_path '/LipidDecontamination'])
    mkdir([Par.Paths.out_path '/LipidMaskNii'])

    x_dim=size(csi,1);
    y_dim=size(csi,2);
    z_dim=size(csi,3);

    lipid_mask_total=zeros(x_dim,y_dim,z_dim);

    f_range=size(csi,4);
    f_range_start=ceil(f_range*0.82);
    f_range_end=ceil(f_range*0.96);

    Text=sprintf('\n\nVectorsize Points of the csi go from 0 to %d!',f_range);
    disp(Text)
    Text=sprintf('\nRemoving Lipids in the Range from %d to %d:\n',f_range_start,f_range_end);
    disp(Text)

    for nCha = 1:size(csi,5)     

        csi_temp = squeeze_single_dim(csi(:,:,:,:,nCha),5); 
        csi_lip = fftshift(fft(csi_temp(:,:,:,:),[],4),4);

        for a=1:x_dim
            for b=1:y_dim
                for c=1:z_dim
                       ftspectra(a,b,c,:)=abs((fftshift(fft(csi_temp(a,b,c,:)))));
                       ftspectramax(a,b,c)=max(ftspectra(a,b,c,f_range_start:f_range_end));
                end
            end
        end

        max_sig=0.30*max(ftspectramax(:));
        sens=ftspectramax;
        sens(sens<max_sig)=0;
        sens(sens~=0)=1;
        lipid_mask=sens;

        struct_elm=strel('disk',1,0);
        tumor_mask=imerode(mask,struct_elm);
        lipid_mask=lipid_mask-tumor_mask;
        lipid_mask(lipid_mask~=1)=0;

        struct_elm=strel('disk',6,0);
        corner_mask=imdilate(mask,struct_elm);
        corner_mask(corner_mask~=1)=2;
        corner_mask(corner_mask~=2)=0;
        corner_mask(corner_mask~=0)=1;

        lipid_mask=lipid_mask-corner_mask;
        lipid_mask(lipid_mask~=1)=0;

        lipid_mask_total=lipid_mask_total+lipid_mask;

        save(sprintf('%s/LipidMaskNii/lipid_mask_chn_%d.mat',Par.Paths.out_path,nCha),'lipid_mask');  

        mkdir([sprintf('%s/LipidDecontamination/Channel_%d', Par.Paths.out_path,nCha)])

        for Slc = 1:size(csi_temp,3)
            Lipid_fig = figure('visible','on');
            imagesc(squeeze(lipid_mask(:,:,Slc)),[0 1])
            colorbar;      
            saveas(Lipid_fig,sprintf('%s/LipidDecontamination/Channel_%d/lipid_mask_slc_%d', Par.Paths.out_path,nCha,Slc),'epsc2')
            saveas(Lipid_fig,sprintf('%s/LipidDecontamination/Channel_%d/lipid_mask_slc_%d', Par.Paths.out_path,nCha,Slc),'fig')            
            close(Lipid_fig)

            Brain_fig = figure('visible','on');
            imagesc(squeeze(abs(csi(:,:,Slc,1,nCha))),[0 max_sig])
            colorbar;      
            saveas(Brain_fig,sprintf('%s/LipidDecontamination/Channel_%d/brain_mask_slc_%d', Par.Paths.out_path,nCha,Slc),'epsc2')
            saveas(Brain_fig,sprintf('%s/LipidDecontamination/Channel_%d/brain_mask_slc_%d', Par.Paths.out_path,nCha,Slc),'fig')
            close(Brain_fig)
        end

        fprintf('Decontaminating channel %d\n', nCha)    

        for Slc = 1:size(csi_temp,3)

            % Prepare slices: Get slice csi data, Extract Relevant Lipids, Get slice mask
            csi_slc = fftshift(fft(squeeze_single_dim(csi_temp(:,:,Slc,:),3),[],3),3);

                % Get L2 Parameter
                param.beta =  Par.Settings.LipidDecon_L2BetaCorrFactor * 10^-15; 
		
                % Get L2 Lipids
                if(1)%Par.CSI.ThreeD_flag || 1)
                    param.Lipid = [];
                    for SlcAlias = 1:size(csi_temp,3)
                        param.Lipid = cat(1, param.Lipid, get_LipidBasis(squeeze(csi_lip(:,:,SlcAlias,:)),lipid_mask(:,:,SlcAlias)));
                    end
                else
                    param.Lipid = get_LipidBasis(csi_slc,lipid_mask(:,:,Slc));
                end

                param.Lipid = transpose(param.Lipid);

                % Get L2 Brainmask
                if(Par.Flags.InterpolateCSIResolution_flag && ~Par.Settings.InterpolateCSIResolution_InkSpace)
                    param.Bmask = mask_BefInterpol(:,:,Slc);
                else
                    param.Bmask = mask(:,:,Slc);            
                end

            % L2 Regularization
            csi_slc = LipidDecon_L2(csi_slc,param);
            csi(:,:,Slc,:,nCha) = ifft(fftshift(csi_slc,3),[],3);

        end

        clear csi_slc csi_lip sens Slc SlcAlias lipid_mask param

    end

    save(sprintf('%s/LipidMaskNii/lipid_mask_total.mat',Par.Paths.out_path),'lipid_mask_total'); 

    mkdir([sprintf('%s/LipidDecontamination/Total', Par.Paths.out_path)])    

    for Slc = 1:size(csi_temp,3)
            Lipid_fig_total = figure('visible','on');
            imagesc(squeeze(lipid_mask_total(:,:,Slc)),[0 nCha])
            colorbar;      
            saveas(Lipid_fig_total,sprintf('%s/LipidDecontamination/Total/lipid_mask_total_slc_%d', Par.Paths.out_path,Slc),'epsc2')
            saveas(Lipid_fig_total,sprintf('%s/LipidDecontamination/Total/lipid_mask_total_slc_%d', Par.Paths.out_path,Slc),'fig')            
            close(Lipid_fig_total)
    end    

    csi_bak.Data = csi;
    csi = csi_bak;
    clear csi_bak;
    
    
end   


  
    


%% Compute Weights w_n


% spring cleaning
clear image_flip image_dummy image_flip_dummy image_VC_dummy fid_mask file_no slice_no

if(~(isfield(csi.Par,'dicom_flag') && csi.Par.dicom_flag))
    
    if(~FastPIReprocess_flag)
        fprintf('\n\nCompute the Coil Combination Weights w_n')  

        weights.Data = 1;
        if(size_csi(5) > 1 && exist('image','var') && Par.Flags.image_VC_flag)                  % Sensmap Weights
            fprintf('\nCompute a sensitivity map and use it for coil combination . . .\n')  

            weights = image;
            weights.Data = conj(image.Data(:,:,:,1,:)) ./ repmat(image_VC.Data(:,:,:,1,:), [1 1 1 1 size(image.Data,5)]);


        elseif(size_csi(5) > 1 && exist('image','var'))									% MUSICAL weights
            fprintf('\nUse the imaging data for coil combination.\n')  																					% (Siemens does that, so that the spectra increase in frequency from right to left)
            weights = image; weights.Data = conj(weights.Data(:,:,:,1,:));
    %         weights.Mask = mask;

        elseif(size_csi(5) > 1 && ~exist('image','var'))                              % 1st FID point weights
            fprintf('\nUse the fourth FID point for coil combination.\n')  																					% (Siemens does that, so that the spectra increase in frequency from right to left)
            weights = csi;
            weights.Data = conj(csi.Data(:,:,:,4,:));

        elseif(size_csi(5) == 1 && exist('image','var'))                              % Weights and for VC data if image was inputted.
            fprintf('\nPhase csi with the imaging data.\n')  																					% (Siemens does that, so that the spectra increase in frequency from right to left)
            weights = image;
            weights.Data = conj(image.Data(:,:,:,1,:)) ./ abs(image.Data);                                              % Only phase VC data. This makes abs(weights) = 1. Thus also scaling = 1. Hence csi data is only phased.

        end

        clear image image_VC

    end

    % Load the water reference weights and overwrite current ones, because: If
    % imaging weights are used
    if(Par.Flags.WaterReference_flag && exist([Par.Paths.out_path '/WaterReference.mat' ],'file') && exist('weights','var'))	% If water ref alrdy exists, load and use its weights
        load([Par.Paths.out_path '/WaterReference.mat' ],'weights')
    end


    % Weight and Sum CSI Data
    fprintf('\n\nWeight and Sum CSI Data')  

    csi = op_CoilCombineData(csi,weights);



    % Rescale csi
    csi.Data = csi.Data * 10^5;
    if(isfield(csi,'NoiseData'))
        csi.NoiseData = csi.NoiseData * 10^5;
    end
end
    

%% DEBUG MODE: PLOT UNPHASED, PHASED & SUMMED SPECTRA

% for point_pair_index = 1:numel(plot_points)
%     for channel_index = 1:numel(plot_channels)
%      
%         
%         plot_point_vec = plot_points{point_pair_index};
%         
%         
%         size(real(squeeze(csi_unphased(channel_index,plot_point_vec(1),plot_point_vec(2),1,:))))
%         
%         figure
%         plot(real(squeeze(csi_unphased(channel_index,plot_point_vec(1),plot_point_vec(2),1,:))))
%         title(sprintf('UNPHASED, x = %d, y = %d, channel = %d',plot_point_vec(1),plot_point_vec(2),plot_channels(channel_index)))
%         
%         waitforbuttonpress
%         
%         figure
%         plot(real(squeeze(csi_phased(channel_index,plot_point_vec(1),plot_point_vec(2),1,:))))
%         title(sprintf('PHASED, x = %d, y = %d, channel = %d',plot_point_vec(1),plot_point_vec(2),plot_channels(channel_index)))
%         
%         waitforbuttonpress
%         
%         figure
%         plot(real(squeeze(csi_summed(plot_point_vec(1),plot_point_vec(2),1,:))))
%         title(sprintf('SUMMED, x = %d, y = %d',plot_point_vec(1),plot_point_vec(2)))
%         
%         waitforbuttonpress
%         close all;
%         
%         
%     end
% end
% 
% clear csi_unphased csi_phased csi_summed


% % % RESHAPE CSI?
% % reshape csi
% csi = reshape(squeeze(csi),[size_csi(2) size_csi(3) size_csi(4) size_csi(5)]);




%%
% Hack for Korbinian: Load Data

% load('Data4KorbinianBeforeLipidDecon.mat');

% % L2 Decontamination
% Par.Flags.LipidDecon_flag = 1; Par.Flags.LipidDecon_L1_flag = false; 
% Par.Settings.LipidDecon_MethodAndNoOfLoops = 'L2,0.1'; Par.Settings.LipidDecon_L2BetaCorrFactor = 0.1;


% % L1 Decontamination
% Par.Flags.LipidDecon_flag = 1; Par.Flags.LipidDecon_L1_flag = true; 
% Par.Settings.LipidDecon_MethodAndNoOfLoops = 'L1,5'; Par.Settings.LipidDecon_NoOfLoops = 5;



%% 14. Bilgic Lipid Decontamination

% % HackForNow:
% csiBak = csi;
% csi = csi.Data;
% 
% if(Par.Flags.LipidDecon_flag == 1 && ~exist([tmp_dir '/Parameters_water.mat'],'file'))
% 	fprintf('\n\nPerform Bilgic Lipid Decontamination.')
% 
% % MASK
% 	if(exist([tmp_dir '/mask_lipid.raw'],'file'))
% 		fid_mask_lipid = fopen([tmp_dir '/mask_lipid.raw'],'r');
%         mask_lipid = reshape(fread(fid_mask_lipid, 'float'), [Par.CSI.nFreqEnc Par.CSI.nPhasEnc Par.CSI.nSLC*Par.CSI.nPartEnc]);
%         fclose(fid_mask_lipid);
% 
% 		if(Par.Flags.InterpolateCSIResolution_flag && ~Par.Settings.InterpolateCSIResolution_InkSpace)  % If we do image-domain interpol --> haven't done interpol yet
% 			fid_mask_lipid = fopen([tmp_dir '/mask_lipid_BefInterpol.raw'],'r');
% 			mask_lipid_BefInterpol = reshape(fread(fid_mask_lipid, 'float'), [Par.CSI.nFreqEnc_BefInterpol Par.CSI.nPhasEnc_BefInterpol Par.CSI.nSLC_BefInterpol*Par.CSI.nPartEnc_BefInterpol]);
% 			fclose(fid_mask_lipid);				
%             mask_lipid = mask_lipid_BefInterpol;
% %         else
% % 			mask_lipid_BefInterpol = mask_lipid;
% 
% 		end
% 
% 		if(sum(sum(sum(mask_lipid))) == 0)
% 			mask_lipid = ones([Par.CSI.nFreqEnc Par.CSI.nPhasEnc Par.CSI.nPartEnc*Par.CSI.nSLC]);
% 		end
%         else
% 		mask_lipid = ones([Par.CSI.nFreqEnc Par.CSI.nPhasEnc Par.CSI.nPartEnc*Par.CSI.nSLC]);
% 
% 	end
% 
% 	% Lipid Mask
% 	lipid_mask = mask_lipid;
% 	%lipid_mask = lipid_mask .* ~mask_BefInterpol;
% 
% 	mkdir([Par.Paths.out_path '/LipidDecontamination'])					% mod stuff gives 1 if true, and 2 if false.
% 	fiddy = fopen(sprintf('%s/LipidDecontamination/UsedL%dNorm.txt',Par.Paths.out_path,Par.Flags.LipidDecon_L1_flag + mod(2,Par.Flags.LipidDecon_L1_flag)),'w+'); 
% 	fclose(fiddy); 
% 
% 	for Slc = 1:size(csi,3)
% %		Scaling_fig = figure('visible','off');
% %		imagesc(1./squeeze(scaling_inv(:,:,Slc,1)),[0 5/scaling_inv(ceil(size(scaling_inv,1)/2),ceil(size(scaling_inv,2)/2),1,1) ])
% %		colorbar;
% %		saveas(Scaling_fig,sprintf('%s/LipidDecontamination/scaling_slc%d', Par.Paths.out_path,Slc),'epsc2')
% %		saveas(Scaling_fig,sprintf('%s/LipidDecontamination/scaling_slc%d', Par.Paths.out_path,Slc),'fig')
% %		close(Scaling_fig)
% 
% 		Lipid_fig = figure('visible','off');
% 		imagesc(squeeze(lipid_mask(:,:,Slc)),[0 1])
%		colorbar;
% 		saveas(Lipid_fig,sprintf('%s/LipidDecontamination/lipid_mask_slc%d', Par.Paths.out_path,Slc),'epsc2')
% 		saveas(Lipid_fig,sprintf('%s/LipidDecontamination/lipid_mask_slc%d', Par.Paths.out_path,Slc),'fig')
% 		close(Lipid_fig)	
% 	end
% 
% 	param.beta = NaN;
% 	csi_lip = fftshift(fft(squeeze(csi(:,:,:,:)),[],4),4);
% 	for Slc = 1:size(csi,3)
% 		
% 		% Prepare slices: Get slice csi data, Extract Relevant Lipids, Get slice mask
% 		fprintf('\n\nDecontaminating slice %d', Slc)
% 		csi_slc = squeeze(csi(:,:,Slc,:));
% 		
% 		csi_slc = fftshift(fft(csi_slc,[],3),3);
% 		
% 
%         if(Par.Flags.SliceParallelImaging_flag == 1)
% 						
% 			[y, x] = find(Slc == Par.Settings.SliceAliasingPattern);
% 			AccessSlc = Par.Settings.SliceAliasingPattern(y,:);
%             param.Lipid = [];
% 			for SlcAlias = AccessSlc
% 				if(SlcAlias ~= Slc)
% 					LipAliasMask = lipid_mask(:,:,SlcAlias);
% 					LipAliasMask = reshape(LipAliasMask, [1 size(LipAliasMask)]);
% 					LipAliasMask1_k = FFTOfMRIData(LipAliasMask,0,[2 3],1);
% 					FoVShift = cat(2,Par.Settings.FoVShifts_x(SlcAlias)-Par.Settings.FoVShifts_x(Slc),Par.Settings.FoVShifts_y(SlcAlias)-Par.Settings.FoVShifts_y(Slc));
% 					LipAliasMask1 = kSpace_FoVShift(LipAliasMask1_k,FoVShift);
% 					LipAliasMask1 = FFTOfMRIData(LipAliasMask1,0,[2 3],0);
% 					LipAliasMask_From2To1 = (LipAliasMask1) .* reshape(mask(:,:,Slc), [1 size(mask,1) size(mask,2)]);
% 					LipAliasMask_From2To1_k = FFTOfMRIData(LipAliasMask_From2To1,0,[2 3],1);
% 					LipAliasMask_From2To1_Final = kSpace_FoVShift(LipAliasMask_From2To1_k,-FoVShift);
% 					LipAliasMask_From2To1_Final = abs(squeeze(FFTOfMRIData(LipAliasMask_From2To1_Final,0,[2 3],0)));
% 					LipAliasMask_From2To1_Final(LipAliasMask_From2To1_Final <= 0.5) = 0; 
% 					LipAliasMask_From2To1_Final(LipAliasMask_From2To1_Final > 0.5) = 1; 
% 				else
% 					LipAliasMask_From2To1_Final = lipid_mask(:,:,SlcAlias);
%             end
% 				param.Lipid = cat(1,   param.Lipid,get_LipidBasis(squeeze(csi_lip(:,:,SlcAlias,:)), LipAliasMask_From2To1_Final)   );
	% end
% 			clear x y AccessSlc SlcAlias;
% %         elseif(Par.CSI.ThreeD_flag)
% %             param.Lipid = [];
% %             for SlcAlias = 1:size(csi,3)
% %                 param.Lipid = cat(1,   param.Lipid,get_LipidBasis(squeeze(csi_lip(:,:,SlcAlias,:)), lipid_mask(:,:,SlcAlias))   );            
% %             end
%         else
% 			param.Lipid = get_LipidBasis(csi_slc,lipid_mask(:,:,Slc));
%         end
	% 
%         if(Par.Flags.InterpolateCSIResolution_flag && ~Par.Settings.InterpolateCSIResolution_InkSpace)
%             param.Bmask = mask_BefInterpol(:,:,Slc);
%         else
%             param.Bmask = mask(:,:,Slc);            
%         end
% 		
% 		% L1 or L2
%         if(Par.Flags.LipidDecon_L1_flag)
% 			
% 			
% 			% L1 
% 			if(Slc == 1)
% 				param_old = param; load('./Bilgic_param_template.mat'); param.Bmask=param_old.Bmask;param.Lipid=param_old.Lipid; param.beta=param_old.beta;
% 				size_csi_bilgic = size(csi);
% 				% Elliptical Filter Mask (zf necessary???)
% 				EllipticalWMask = repmat(EllipticalFilter(ones(size(csi,1),size(csi,2)),[1 2],[1 1 1 size(csi,1)/2+1],1),[1 1 size_csi_bilgic(4)]);
% 				param.FT = FT_v2(EllipticalWMask);	
% 			end
% 			param.data =  param.FT*csi_slc;
% 
% 			for t = 1:Par.Settings.LipidDecon_NoOfLoops
% 				fprintf('\nDecontamination Loop %d of %d. Norm Results:\n',t,Par.Settings.LipidDecon_NoOfLoops)
% 				csi_slc = lipid_suppression(csi_slc,param);
% 			end
% 		else
% 			
% 			
% 			% L2
% 			param.beta =  Par.Settings.LipidDecon_L2BetaCorrFactor * 10^-15; 
% 			param.Lipid = transpose(param.Lipid);
% 			csi_slc = LipidDecon_L2(csi_slc,param);
%         end
% 		
% 		
% 		
% 		
% 		
% 		csi(:,:,Slc,:) = ifft(fftshift(csi_slc,3),[],3);
% 	end
% 
% 	%gh save lipid mask!
% 	param_beta = param.beta;
% 	save([Par.Paths.out_path '/LipidDecontamination/Lipid_mask.mat'], 'lipid_mask','param_beta');
% 
% 	clear csi_slc param EllipticalWMask size_csi_bilgic Scaling_fig Lipid_fig scaling_inv2 lipid_mask Slc t param_beta
% 
% 	% % For Visualization
% 	% % Lipid Decon
% 	% csi_LipidDecon = csi;
% 	% csi_fft = fftshift(fft(csi_LipidDecon,[],3),3);
% 	% for t = 1:10
% 	%     csi_LipidDecon = lipid_suppression(csi_LipidDecon,param);
% 	%     csi_LipidDecon_fft = fftshift(fft(csi_LipidDecon,[],3),3);
% 	%     figure;
% 	%     imagesc( [ mask.*squeeze(sum(abs(csi_fft(:,:,1199:1559)),3))  mask.*squeeze(sum(abs(csi_LipidDecon(:,:,1199:1559)),3))] ), ...
% 	%     axis image, colorbar, drawnow
% 	%     title(['t = ' num2str(t)] )
% 	% end
% 
% 	% lipid_noLipidDecon = abs(fftshift(fft(csi,[],3),3));
% 	% lipid_noLipidDecon = squeeze(sum(lipid_noLipidDecon(:,:,1199:1559),3));
% 	% 
% 	% lipid_LipidDecon = abs(fftshift(fft(csi_LipidDecon,[],3),3));
% 	% lipid_LipidDecon = squeeze(sum(lipid_LipidDecon(:,:,1199:1559),3));
% 
% end
% 
% csiBak.Data = csi; csi = csiBak; clear csiBak;


%% Frequency Alignment


if(Par.Flags.AlignFreq_flag)
	fprintf('\n\nPERFORM FREQUENCY ALIGNMENT')  
    if(exist('AlignFreq_B0Fieldmap','var'))
        % Convert from ppm to Hz
        AlignFreq_B0Fieldmap_Hz = -AlignFreq_B0Fieldmap*Par.CSI.LarmorFreq/10^6;
        % time vector
        nS  = Par.CSI.vecSize;
        dwt = Par.CSI.Dwelltimes(1);
        t   = (0:nS-1)*dwt/10^9;
        csi.Data = odMRSIrecon(csi.Data,AlignFreq_B0Fieldmap_Hz,t);
        save(sprintf('%s/AlignFreq/AlignFreq_ShiftMap.mat', Par.Paths.out_path), 'AlignFreq_B0Fieldmap')
        clear nS dwt

        
    else
        InData.csi = csi.Data; 
        if(Par.Flags.InterpolateCSIResolution_flag && ~Par.Settings.InterpolateCSIResolution_InkSpace)  % In case we do image-domain interpolation --> Havent done interpol yet
            InData.mask = mask_BefInterpol;                                                             % --> use BefInterpol
        else
            InData.mask = mask;        
        end
        if(exist([Par.Paths.out_path '/scalings/AlignFreq_ShiftMap.mat'],'file'))
            load([Par.Paths.out_path '/scalings/AlignFreq_ShiftMap.mat'],'FreqAlignRefSpec') 	% Make sure that we align all the averages not just within each other, but also with each other.
        else																		% Therefore use a reference spectrum from prev
            FreqAlignRefSpec = [];
        end

        FreqAlignSettings.LarmorFreq = Par.CSI.LarmorFreq;
        FreqAlignSettings.Dwelltime = Par.CSI.Dwelltimes(1);
        FreqAlignSettings.vecsize = size(csi,4);
        FreqAlignSettings.PeakSearchPPM = 4.65;
        FreqAlignSettings.PeakSearchRangePPM = 0.2;
        [csi.Data, AlignFreq_ShiftMap,FreqAlignRefSpec] = FrequencyAlignment(InData,FreqAlignSettings,4,2,FreqAlignRefSpec);
        
        % Print figures
        for Slc = 1:size(csi.Data,3)
            Scaling_fig = figure('visible','off');
            imagesc(AlignFreq_ShiftMap(:,:,Slc))
            colorbar;
            saveas(Scaling_fig,sprintf('%s/AlignFreq/AlignFreq_ShiftMap_slc%d', Par.Paths.out_path,Slc),'epsc2')
            saveas(Scaling_fig,sprintf('%s/AlignFreq/AlignFreq_ShiftMap_slc%d', Par.Paths.out_path,Slc),'fig')
            close(Scaling_fig)
        end	
        save(sprintf('%s/AlignFreq/AlignFreq_ShiftMap.mat', Par.Paths.out_path), 'AlignFreq_ShiftMap','FreqAlignRefSpec')

        % Free memory
        InData.csi = 0;
    end
	
end



%% Nuisance Removal

if(Par.Flags.NuisRem_flag == 1)
    run ./NuisanceRemoval_HSVD.m
end


% %% Apply Hamming Filter to Combined MRSI Data
% 
% if(Par.Flags.hamming_flag)
% 	if( ~isfield(csi.Par,'SpatialSpectralEncoding_flag') || ~csi.Par.SpatialSpectralEncoding_flag)      % For SpatialSpectralEncoding we have done it already during reco
%     
% 	    fprintf('\n\nAPPLY HAMMING FILTER')  % If Full_ElliptWeighted_Or_Weighted_Acq = 4, the data is alrdy intrinsically filtered in z-dimension
% 	    % Was the WeightedAcquisition undone?
% 	    UndoWeightedAcq_flag = false;
% 	    if(isfield(csi,'RecoSteps'))
%             fieldies = fieldnames(csi.RecoSteps); field = fieldies(~cellfun(@isempty,regexpi(fieldies,'op_AverageMRData')));
%             if(~isempty(field))
%                 UndoWeightedAcq_flag = csi.RecoSteps.(field{1}).UndoWeightedAveraging_flag;
%             end
% 	    end
% 	    if(size_csi(3) > 1 && Par.CSI.ThreeD_flag && (Par.CSI.Full_ElliptWeighted_Or_Weighted_Acq ~= 4 || UndoWeightedAcq_flag))      
%             csi.Data = HammingFilter(csi.Data,[1 2 3],Par.Settings.hamming_factor,'OuterProduct',0);     
%         else
%             csi.Data = HammingFilter(csi.Data,[1 2],Par.Settings.hamming_factor,'OuterProduct',0);
% 	    end 
%     elseif(isfield(csi.Par,'SpatialSpectralEncoding_flag') && csi.Par.SpatialSpectralEncoding_flag && size_csi(3) > 1 && Par.CSI.ThreeD_flag)
% 		csi.Data = HammingFilter(csi.Data,[3],Par.Settings.hamming_factor,'OuterProduct',0);   
% 
% 	end
% end


%% DEBUG MODE: SHOW UNHAMMINGED AND HAMMINGED MAGNITUDE

% csi_hamminged = squeeze(csi(:,:,1,1));
% 
% figure
% imagesc(abs(csi_unhamminged))
% title('UNHAMMINGED')
% 
% figure
% imagesc(abs(csi_hamminged))
% title('HAMMINGED')
% 
% 
% clear csi_hamminged csi_unhamminged






%% Perform Interpolation

% In image domain
if(Par.Flags.InterpolateCSIResolution_flag && ~Par.Settings.InterpolateCSIResolution_InkSpace)
	csi.Data = InterpolateCSIImages(csi.Data,0,mask_BefInterpol,Par.Settings.InterpolateResolutionRatio);
    csi.RecoPar.DataSize = size(csi.Data);
end


%% Time Domain Interpolation

% Interpolation Time Domain
% In time domain __________________________________________________________________________________
if(Par.Flags.TimeInterpolation_flag)
	display([ char(10) 'PERFORM INTERPOLATION IN TIME-DOMAIN' char(10)])
	csi.Data = TimeInterpolation(csi.Data,Par.Settings.TimeInterpolationFactor(1),Par.Settings.TimeInterpolationFactor(2),Par.Settings.TimeInterpolationFactor(3));
	display(['TruncateFactor = ', num2str(Par.Settings.TimeInterpolationFactor(1)), ', ZerofillFactor = ', num2str(Par.Settings.TimeInterpolationFactor(2)), '.'])
	if(Par.Settings.TimeInterpolationFactor(3) == 1)
	display(['Filling to original size.'])	
	end

end
		display(['Done.']) 		


%% Apply exponential Time-Domain Filter to Combined MRSI Data

if (Par.Flags.exponential_filter_Hz_flag)
    fprintf('\nAPPLY EXP TIME-DOMAIN FILTER\n')
	csi = ExponentialFilter(csi, Par.CSI.Dwelltimes,Par.Settings.exponential_filter_Hz,4);
end





%% Perform same things on noise as on CSI
% Do that for the non-PI and the PI noise seperately, in order to calculate the g-factor map


if(exist('noise_sim','var'))
    fprintf('\n\nPerform Same Calculations on Simulated Noise for SNR Calculation')

    % Hack hardcode number of Pseudo Replicas 
    NoOfPseudoReplicas = 100; 
    % Parallel Imaging Performed in ParallelImagingSimReco.m


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%   NON-PI NOISE   %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Spatial FFT
    noise_sim = ifftshift(ifftshift(noise_sim,2),3);
    noise_sim = fft(fft(noise_sim,[],2),[],3);
    if(Par.CSI.ThreeD_flag && size(noise_sim,3) > 1)
        noise_sim = fft(ifftshift(noise_sim,4),[],4);
    end
    noise_sim = conj(noise_sim);
    noise_sim = fftshift(fftshift(noise_sim,2),3);
    if(Par.CSI.ThreeD_flag && size(noise_sim,3) > 1)
        noise_sim = fftshift(noise_sim,4);
    end
    noise_sim = flip(noise_sim,2);

	
	
	
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%    BOTH NOISES   %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%	
	% Coil Comb, Scaling, Hamming Filtering, Exponential Filter
	LoopOverNoises = {'noise_sim'};
	if(Par.Flags.TwoDCaipParallelImaging_flag || Par.Flags.SliceParallelImaging_flag)	
		LoopOverNoises{2} = 'noise_sim_PI';
	end
	
	
	for CurNoise = LoopOverNoises
		% Process current noise
		CurNoiseData = eval(CurNoise{:});
		
		% Coil Combination
		if(exist('weights','var')) 
			CurNoiseData = CurNoiseData .* repmat(weights,[1 1 1 1 NoOfPseudoReplicas]);
			CurNoiseData = sum(CurNoiseData,1); 
		end
		CurNoiseData = reshape(squeeze(CurNoiseData),[size(CurNoiseData,2) size(CurNoiseData,3) size(CurNoiseData,4) NoOfPseudoReplicas]);

		% Scaling
		if(exist('scaling_inv','var'))
			CurNoiseData = CurNoiseData ./ scaling_inv(:,:,:,1); 
		end
		
		
		% Lipid Decon (Can this be simulated?)

		
		
		% Frequency Alignment
		if(Par.Flags.AlignFreq_flag)
            if(exist('AlignFreq_B0Fieldmap','var'))
                CurNoiseData = odMRSIrecon(CurNoiseData,AlignFreq_B0Fieldmap_Hz,t);
            else
                InData.csi = CurNoiseData;
                CurNoiseData = FrequencyAlignment(InData,AlignFreq_ShiftMap,4,2);       % Only apply frequency shifts, don't calculate them again!
            end
		end

		
		
		% Hamming Filtering
		if(Par.Flags.hamming_flag)
			if(size(CurNoiseData,3) > 1 && Par.CSI.ThreeD_flag)
				CurNoiseData = HammingFilter(CurNoiseData,[1 2 3],Par.Settings.hamming_factor,'OuterProduct',0);    
			else
				CurNoiseData = HammingFilter(CurNoiseData,[1 2],Par.Settings.hamming_factor,'OuterProduct',0);
			end 
		end
		
		
		% Interpolation
        % In Image Domain
		if(Par.Flags.InterpolateCSIResolution_flag && ~Par.Settings.InterpolateCSIResolution_InkSpace)
			CurNoiseData = InterpolateCSIImages(CurNoiseData,0,mask_BefInterpol,Par.Settings.InterpolateResolutionRatio);
		end		
		
		% Interpolation Time Domain
		if(Par.Flags.TimeInterpolation_flag)
			display(['Perform FID ZF on Noise'])
			Par.Settings.TimeInterpolationFactor(1)
			Par.Settings.TimeInterpolationFactor(2)
			Par.Settings.TimeInterpolationFactor(3)  			
			CurNoiseData = TimeInterpolation(CurNoiseData,Par.Settings.TimeInterpolationFactor(1),Par.Settings.TimeInterpolationFactor(2),Par.Settings.TimeInterpolationFactor(3));
		end
		
		
		% Exponential Filter
		if (Par.Flags.exponential_filter_Hz_flag)
			CurNoiseData = ExponentialFilter(CurNoiseData,Par.CSI.Dwelltimes,Par.Settings.exponential_filter_Hz,4);
		end		
		
		
		% Write data back 
		eval([CurNoise{:} ' = CurNoiseData;']);		
	end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%     PI NOISE      %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(Par.Flags.TwoDCaipParallelImaging_flag || Par.Flags.SliceParallelImaging_flag)
        % In case of PI: The SNR has to be calculated based on the PI-noise
        noise_sim_spectral = noise_sim_PI;
        noise_sim_spectral = fftshift(fft(noise_sim_spectral,[],5),5) / sqrt(2*size(noise_sim_spectral,5)); % Why is this sqrt(2) necessary?  BUG: 5 --> 4 [*]
        noise_sim_time = noise_sim_PI * 0.1370;                                                             % But without it, the spectrum is differently scaled than the LCModel spectrum.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      % COMMENT ON [*]: By performing an fft and dividing by sqrt(size), the std is not changed. Therefore it doesnt matter if the bug is fixed or not. 
    %%%%%    NON-PI NOISE    %%%%%      % But the sqrt(2) has an effect. What it is for, I don't know.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	else
        % In case of no PI: SNR calculated based on normal noise
        noise_sim_spectral = noise_sim;
        noise_sim_spectral = fftshift(fft(noise_sim_spectral,[],5),5) / sqrt(2*size(noise_sim_spectral,5)); % BUG: 5 --> 4
        noise_sim_time = noise_sim * 0.1370;   % 0.1492 (originally), 0.137 (now): To make a NAA peak that has amplitude 1 in MATLAB have the concentration of 1 in LCModel.
    end



    % Compute NoiseScaling Matrix
    if(Par.Flags.TwoDCaipParallelImaging_flag)
        NoiseScalingMatrix_PI = mean( cat(4, std(real(noise_sim_PI),0,4) , std(imag(noise_sim_PI),0,4)), 4);
        NoiseScalingMatrix_Fully = mean( cat(4, std(real(noise_sim),0,4) , std(imag(noise_sim),0,4)), 4);	% Use info of real and imaginary part, since both should have the same properties.  
        g_FactorMap = NoiseScalingMatrix_PI ./ (NoiseScalingMatrix_Fully * sqrt(R));						% See Robson et al., "Comprehensive Quantification of Signal-to-Noise Ratio...",
        g_FactorMap = g_FactorMap .* mask;																% MRM 60:895-907, 2008, Equation [5]. Assumption: Signal is same in Accelerated & Non-A.
    end
    NoiseScalingMatrix_spectral = mean( cat(4, std(real(noise_sim_spectral),0,4) , std(imag(noise_sim_spectral),0,4)), 4);
    NoiseScalingMatrix_time = mean( cat(4, std(real(noise_sim_time),0,4) , std(imag(noise_sim_time),0,4)), 4);

    
    
    % COMPUTE Internal NoiseScalingMatrix
%     noise_sim2 = csi(:,:,:,end-150:end) / sqrt(2);
%     noise_sim3 = fftshift(fft(noise_sim2,[],4),4) / sqrt(size(noise_sim2,4));
%     NoiseScalingMatrix_spectral2 = mean( cat(4, std(real(noise_sim2),0,4) , std(imag(noise_sim2),0,4)), 4);
%     NoiseScalingMatrix_spectral3 = mean( cat(4, std(real(noise_sim3),0,4) , std(imag(noise_sim3),0,4)), 4);
    

    clear scaling_inv noise_sim noise_sim_PI noise_sim_time noise_sim_spectral



end
%% Hack for new Data structure 
if(isfield(csi,'NoiseData'))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%     PI NOISE      %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(Par.Flags.TwoDCaipParallelImaging_flag || Par.Flags.SliceParallelImaging_flag)
        % In case of PI: The SNR has to be calculated based on the PI-noise
        noise_sim_spectral = noise_sim_PI;
        noise_sim_spectral = fftshift(fft(noise_sim_spectral,[],5),5) / sqrt(2*size(noise_sim_spectral,5)); % Why is this sqrt(2) necessary?  BUG: 5 --> 4 [*]
        noise_sim_time = noise_sim_PI * 0.1370;                                                             % But without it, the spectrum is differently scaled than the LCModel spectrum.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      % COMMENT ON [*]: By performing an fft and dividing by sqrt(size), the std is not changed. Therefore it doesnt matter if the bug is fixed or not. 
    %%%%%    NON-PI NOISE    %%%%%      % But the sqrt(2) has an effect. What it is for, I don't know.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	else
        % In case of no PI: SNR calculated based on normal noise
        noise_sim_spectral = csi.NoiseData;
        noise_sim_spectral = fftshift(fft(noise_sim_spectral,[],5),5) / sqrt(2*size(noise_sim_spectral,5)); % BUG: 5 --> 4
        noise_sim_time = csi.NoiseData * 0.1370;   % 0.1492 (originally), 0.137 (now): To make a NAA peak that has amplitude 1 in MATLAB have the concentration of 1 in LCModel.
    end



    % Compute NoiseScaling Matrix
    if(Par.Flags.TwoDCaipParallelImaging_flag)
        NoiseScalingMatrix_PI = mean( cat(4, std(real(noise_sim_PI),0,4) , std(imag(noise_sim_PI),0,4)), 4);
        NoiseScalingMatrix_Fully = mean( cat(4, std(real(csi.NoiseData),0,4) , std(imag(csi.NoiseData),0,4)), 4);	% Use info of real and imaginary part, since both should have the same properties.  
        g_FactorMap = NoiseScalingMatrix_PI ./ (NoiseScalingMatrix_Fully * sqrt(R));						% See Robson et al., "Comprehensive Quantification of Signal-to-Noise Ratio...",
        g_FactorMap = g_FactorMap .* mask;																% MRM 60:895-907, 2008, Equation [5]. Assumption: Signal is same in Accelerated & Non-A.
    end
    NoiseScalingMatrix_spectral = mean( cat(4, std(real(noise_sim_spectral),0,4) , std(imag(noise_sim_spectral),0,4)), 4);
    NoiseScalingMatrix_time = mean( cat(4, std(real(noise_sim_time),0,4) , std(imag(noise_sim_time),0,4)), 4);

    
    
    % COMPUTE Internal NoiseScalingMatrix
%     noise_sim2 = csi(:,:,:,end-150:end) / sqrt(2);
%     noise_sim3 = fftshift(fft(noise_sim2,[],4),4) / sqrt(size(noise_sim2,4));
%     NoiseScalingMatrix_spectral2 = mean( cat(4, std(real(noise_sim2),0,4) , std(imag(noise_sim2),0,4)), 4);
%     NoiseScalingMatrix_spectral3 = mean( cat(4, std(real(noise_sim3),0,4) , std(imag(noise_sim3),0,4)), 4);
    

    clear scaling_inv csi.NoiseData noise_sim_PI noise_sim_time noise_sim_spectral
end 


%% Correct for 1storderphase
%     Rationale: Only do 1stPhCorr if: - FID sequence (--> TE < 10000 us = 10ms) and - The user asks for it or it is the water reference scan
    if( (Par.Flags.FirstOrderPhaseCorr_flag || (strcmp(Par.Settings.WaterReference_method, 'W1') && ~exist([Par.Paths.out_path '/WaterReference.mat' ],'file'))) && min(Par.CSI.TEs) < 10000 )
           size_csi=size(csi.Data);
           fprintf('\n\nDoFirstOrderPhaseCorr')
           omega_step = 2*pi/((size_csi(4)-1)*Par.CSI.Dwelltimes(1)/10^9);
           if(~mod(size_csi(4),2))
               
               FreqVec = repmat(-ceil(size_csi(4)/2) : (ceil(size_csi(4)/2)-1), [numel(Par.CSI.TEs(1)) 1]);
           else
               
               FreqVec = repmat(-ceil((size_csi(4)-1)/2) : (ceil((size_csi(4)-1)/2)), [numel(Par.CSI.TEs(1)) 1]);
           end
            TEs = repmat(Par.CSI.TEs(1)/1000, [1 size(FreqVec,2)]);
            FirstOrderCorrVec = exp(-1i * (TEs/1000)*omega_step .* FreqVec);
            clear TEs FreqVec
            FirstOrderCorrVec = reshape(FirstOrderCorrVec, [1 1 size(FirstOrderCorrVec,1) size(FirstOrderCorrVec,2)]);

            csi_spec = fftshift(fft(csi.Data,[],4),4);
            csi_spec = csi_spec .* repmat(FirstOrderCorrVec, [size_csi(1) size_csi(2) size_csi(3) 1]);
            csi.Data = ifft(ifftshift(csi_spec,4),[],4);

    end



%% Manipulate 1storderphase for AD-"change" (GH 2019)
if(Par.Flags.FirstOrderPhaseModulation_flag)
	display([ char(10) char(10) 'Truncate FID to simulate different AD'])
	Cut_time=Par.Settings.FID_Truncation_in_ms   %should be in ms here
	ReadInInfo.Par.Dwelltime
	Cut_length=floor(Cut_time/(ReadInInfo.Par.Dwelltime*2/10^6)) %just calc it as ms %*2 bc displayed BW is half of the right value (verified vy manual calc)
	Cut_length=Cut_length+1; %account for starting at 1

    %csi_spec = fftshift(fft(csi,[],4),4);
    csi.Data = csi.Data(:,:,:,Cut_length:size(csi.Data,4));
    %csi = ifft(ifftshift(csi_spec,4),[],4);
end


% % % % % USE MEANINGFUL IF-CONDITION HERE
if(false)
    csi.Data=circshift(csi.Data,[1 -1 0 0]);
    csi.Data=Zerofilling_Spectral(csi.Data,[size(csi.Data,1) size(csi.Data,2) size(csi.Data,3) size(csi.Data,4)*2]);
% 	mask=zeros(size(mask));
%     mask(31,23,25)=1;
end



%% Handle Water Reference

    % Initialize WatRef variables
    IsWatRef = false;            % Flag telling if the current processed data set "csi" is a water-reference (either W1- or W2-method).
    IsW2WatRef = false;            % Flag telling if it is W1- or W2-method

    if(Par.Flags.WaterReference_flag)

        if(~exist([Par.Paths.out_path '/WaterReference.mat' ],'file') )
            IsWatRef = true;

            % W1 & W2 --> If water does not exist yet, save it
            WaterReferenceCSI = csi;
            save([Par.Paths.out_path '/WaterReference.mat' ], 'WaterReferenceCSI'); clear WaterReferenceCSI
            if(exist('weights','var'))
                save([Par.Paths.out_path '/WaterReference.mat' ],'weights','-append')
            end

            % W2 --> Prepare appropriate LCModel processing
            if(strcmp(Par.Settings.WaterReference_method, 'W2'))
                display([ char(10) 'Water spectra are out_path' char(10)])
                out_dir_spectra = [ Par.Paths.out_path '/water_spectra' ];
                if(isfield(Par.Paths,'LCM_Control_Water_path') )
                    display([ char(10) 'Control file is water file' char(10)])
                    ControlInfo = Par.Paths.LCM_Control_Water_path;
                end
                delete(sprintf('%s/Parameters_water.mat', tmp_dir))            % No need for it anymore
                IsW2WatRef = true;
            end

        % W1 --> if water does exist, load it
        elseif(strcmp(Par.Settings.WaterReference_method, 'W1'))
            display([ char(10) 'Load water reference' char(10)])
            load([Par.Paths.out_path '/WaterReference.mat' ],'WaterReferenceCSI')
            Data.watref = WaterReferenceCSI;
        end

    end


%% 22. Perform averaging

if(~IsWatRef)		% No averaging for watref
	IsLastAvg = CurAvg == numel(Par.Paths.csi_path);
	if(exist([Par.Paths.out_path '/CombinedCSI.mat'],'file'))
		csi_bak = csi.Data;
		load([Par.Paths.out_path '/CombinedCSI.mat'],'csi')
		csi.Data = csi.Data + csi_bak;
	end
	if(IsLastAvg)
		csi.Data = csi.Data/numel(Par.Paths.csi_path);
		if(exist('NoiseScalingMatrix_spectral','var'))
			NoiseScalingMatrix_spectral = NoiseScalingMatrix_spectral * sqrt(numel(Par.Paths.csi_path));
		end
		if(exist('NoiseScalingMatrix_time','var'))
			NoiseScalingMatrix_time = NoiseScalingMatrix_time * sqrt(numel(Par.Paths.csi_path));
		end
		clear csi_bak
	end
end

%% Write the files necessary for LCM-processing

% Dont do for W1-WatRef
if( IsW2WatRef || (~IsWatRef && IsLastAvg) )  % i.e. either we have W2-Watref, or no Watref and last Avg

	display([ char(10) char(10) 'WRITE LCM-FILES'])

	% Dont overwrite if variables were alrdy set previously
	if(~exist('out_dir_spectra','var'))	
		out_dir_spectra = [ Par.Paths.out_path '/spectra' ];
	end
	if(~exist('ControlInfo','var'))
		if(isfield(Par.Paths,'LCM_ControlPath'))
			ControlInfo = Par.Paths.LCM_ControlPath;
		else
			ControlInfo = 0;
		end
	end
	
	Paths.out_dir = out_dir_spectra;
	Paths.basis_file = Par.Paths.basis_path;                            % Path to Basis set that should be used ('.basis')
	Paths.LCM_ProgramPath = Par.Paths.LCM_Path;                         % Path of the LCModel program
	Paths.batchdir = tmp_dir;
	MetaInfo.DimNames = {'x','y','z'};
	MetaInfo.pat_name = Par.CSI.PatName;                                % Name of Patient which determines the naming of the output spectra
	MetaInfo.LarmorFreq = Par.CSI.LarmorFreq;
	MetaInfo.dwelltime = Par.CSI.Dwelltimes(1);
    
    %BOW - define path containing priors
	if(Par.Flags.priors_flag)
        fprintf('\nUse priors from OFF spectra.\nWARNING: THIS OPTION WORKS ONLY FOR USING MEGA-OFF PRIOR KNOWLEDGE FOR MEGA-DIFF-FITTING\n' )
		Paths.priors_dir = Par.Paths.priors_path;
	end	
    
	Data.csi = csi.Data;
	Write_LCM_files(Data,Paths,MetaInfo,ControlInfo,mask,Par.ServerInfo.RunLCModel_CPUCores) 

end


%% Save the processed MRSI data for SNR-Computation after LCModel processing

fprintf('\n\nThe MRSI Pre-Processing and LCModel preparations took %10.6f s.\n',toc(CoilCombtic))
clearvars -except Par csi weights image image_VC mask g_FactorMap NoiseScalingMatrix_spectral NoiseScalingMatrix_time IsWatRef

if(~IsWatRef)		% No need to save water reference
	save([Par.Paths.out_path '/CombinedCSI.mat'], '-v7.3')
end



