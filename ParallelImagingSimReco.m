%% 0. Preparations

% Housekeeping


% Definitions


% Minimum number of SrcPts in Reco.
MinKernelSrcPts = 20;


% How much of the imaging-kspace should be used as ACS data (e.g. 128x128 or 64x64)
ACS_size = 128;




%% 1. Read ACS-Data

% IMAGE NORMAL
if(Par.Flags.image_normal_flag)
	for file_no = 1:numel(Par.Paths.image_normal_path)
        [image_dummy, image_dummy_k] = ...
        read_image(Par.Paths.image_normal_path{file_no},0,'ZeroFilling',0,Par.Image.FreqDirShift(1),0,0,Noise_CorrMat,Par.Settings.phase_encod_dir_is);
          
        if(numel(Par.Paths.image_normal_path) > 1 && size(image_dummy_k,4) > 1)                                 % If e.g. 4 images were acquired, each with different echo time and with 4 slices.
            image_k(:,:,:,file_no) = image_dummy_k(:,:,:,file_no);
        elseif(numel(Par.Paths.image_normal_path) > 1 && size(image_dummy_k,4) == 1)                            % If e.g. 4 images were acquired, each with different echo time but only 1 slice.
            image_k(:,:,:,file_no) = image_dummy_k;  
        elseif(numel(Par.Paths.image_normal_path) == 1)                                                       % If either 1 image with 4 Slices (each slice with different echo time) or single slice was acquired.
            image_k = image_dummy_k;               
        end
	end
    image_k = image_k(:,1:2:end,:,:);
end

if(~exist('image_k','var'))
    display(['Problem: You need to provide imaging data for usage as Auto Calibration Signals,' char(10) 'if GRAPPA/CAIPIRINHA Reconstruction should be performed.'])
    quit force
end


% IMAGE FLIP
if(Par.Flags.image_flip_flag)
    for file_no = 1:numel(Par.Paths.image_flip_path)
        [image_flip_dummy, image_flip_dummy_k] = ...
        read_image(Par.Paths.image_flip_path{file_no},0,'ZeroFilling',1,Par.Image.FreqDirShift(1),0,0,Noise_CorrMat,Par.Settings.phase_encod_dir_is);
        
        if(numel(Par.Paths.image_flip_path) > 1 && size(image_flip_dummy_k,4) > 1) 
            image_flip_k(:,:,:,file_no) = image_flip_dummy_k(:,:,:,file_no);
        elseif(numel(Par.Paths.image_flip_path) > 1 && size(image_flip_dummy_k,4) == 1) 
            image_flip_k(:,:,:,file_no) = image_flip_dummy_k;  
        elseif(numel(Par.Paths.image_flip_path) == 1)
            image_flip_k = image_flip_dummy_k;               
        end        
    end
    image_flip_k = image_flip_k(:,1:2:end,:,:);
    image_k = ( image_k + image_flip_k )/2;
end



ACS = image_k(:,(end/2 + 1) - ACS_size/2 : (end/2 + 1) + ACS_size/2 - 1, (end/2 + 1) - ACS_size/2 : (end/2 + 1) + ACS_size/2 - 1,:);

clear image_dummy image_dummy_k file_no image_flip_dummy image_flip_dummy_k image_flip_k image_k





%% 2. Define 


% Flag Defining if ParallelImaging should be simulated by skipping part of k-space or if the measured data is already undersampled.
% These should be replaced one day by the header.
%TwoDCaipParallelImaging_SimulationFlag = 1;
TwoDCaipParallelImaging_SimulationFlag = Par.Flags.TwoDCaipParallelImaging_flag && numel(csi_k(1,:,:,1,1,1,1)) / sum(sum(~(csi_k(1,:,:,1,1,1,1) == 0))) < 1.5;
%SliceParallelImaging_SimulationFlag = 0;
WholeSliceIsZero = false;
for slice = 1:size(csi_k,4)
	if(WholeSliceIsZero)
		break;
	end
	WholeSliceIsZero = WholeSliceIsZero && (sum(sum(~(csi_k(1,:,:,slice,1) == 0))) == 0);		% ~(csi_k(1,:,:,slice,1) == 0					--> Elements that are non-zero
end																								% sum(sum(~(csi_k(1,:,:,slice,1) == 0)))		--> number of non-zero elements
SliceParallelImaging_SimulationFlag = (size(ACS,4) == size(csi_k,4)) && ~WholeSliceIsZero ...	% (sum(sum(~(csi_k(1,:,:,slice,1) == 0))) == 0)	--> check if number of non-zero els is zero
&& Par.Flags.SliceParallelImaging_flag;

% Make SliceAliasing Patterns like [1 2], [1 3; 2 4], [1 4; 2 5;3 6] ...
if(Par.Flags.SliceParallelImaging_flag && ~isfield('Par.Settings','SliceAliasingPattern'))
	for AliLoop = 1:numel(Par.Settings.FoVShifts_x)/2
        Par.Settings.SliceAliasingPattern(AliLoop,:) = [AliLoop AliLoop+numel(Par.Settings.FoVShifts_x)/2];
	end
end
if(Par.Flags.SliceParallelImaging_flag)
	R_Slice = size(Par.Settings.SliceAliasingPattern,2);
else
	R_Slice = 1;
end


% Compute Mask of Measured Fully Sampled kPoints
if(TwoDCaipParallelImaging_SimulationFlag)
    MeasuredkPoints_FullySampled = logical(squeeze(sum(abs(csi_k(:,:,:,:,1)),1))); 
	MeasuredkPoints_FullySampled_NoEllip = ones([size(csi_k,2) size(csi_k,3) size(csi_k,4)]);
else
    MeasuredkPoints_FullySampled_NoEllip = ones([size(csi_k,2) size(csi_k,3) size(csi_k,4)]);
    if(Par.CSI.Full_ElliptWeighted_Or_Weighted_Acq == 2)
        MeasuredkPoints_FullySampled = EllipticalFilter(MeasuredkPoints_FullySampled_NoEllip,[1 2],[1 1 1 size(csi_k,2)/2-1]);
    end
end

if(~SliceParallelImaging_SimulationFlag)
    MeasuredkPoints_FullySampled = repmat(MeasuredkPoints_FullySampled, [1 1 R_Slice]);
    MeasuredkPoints_FullySampled_NoEllip = repmat(MeasuredkPoints_FullySampled_NoEllip, [1 1 R_Slice]);	
end


% Save measured data for later usage
csi_k_meas = csi_k(:,:,:,:,1);
csi_meas = flipdim(fftshift(fftshift(conj(fft(fft(ifftshift(ifftshift(csi_k_meas,2),3),[],2),[],3)),2),3),2);






%% 2. SliceSimReco: Perform Slice Aliasing

if(Par.Flags.SliceParallelImaging_flag)

    ACS_Slice = ACS;
    if(Par.Flags.TwoDCaipParallelImaging_flag)
        ACS = zeros([size(ACS_Slice,1) size(ACS_Slice,2) size(ACS_Slice,3) size(ACS_Slice,4)/R_Slice]);
    end
    if(SliceParallelImaging_SimulationFlag)
        csi_k_dummy = zeros([size(csi_k,1) size(csi_k,2) size(csi_k,3) size(csi_k,4)/R_Slice size(csi_k,5)]);
    end
    
	for AliLoop = 1:size(Par.Settings.SliceAliasingPattern,1)
        CurSliceAliased = Par.Settings.SliceAliasingPattern(AliLoop,:);

        % For being more memory efficient, do this channel by channel, and parfor
        for cha = 1:size(csi_k,1)
            % Alias CSI Data  ( S I M U L A T I O N )
            if(SliceParallelImaging_SimulationFlag)
                csi_k_dummy2 = csi_k(cha,:,:,CurSliceAliased,:);
                csi_k_dummy(cha,:,:,AliLoop,:) = sum(kSpace_FoVShift(csi_k_dummy2,cat(2,Par.Settings.FoVShifts_x(CurSliceAliased)',Par.Settings.FoVShifts_y(CurSliceAliased)')),4);
            end
            % Alias ACS data
            if(Par.Flags.TwoDCaipParallelImaging_flag)
                ACS(cha,:,:,AliLoop,:) = sum(kSpace_FoVShift(ACS_Slice(cha,:,:,CurSliceAliased,:),cat(2,Par.Settings.FoVShifts_x(CurSliceAliased)',Par.Settings.FoVShifts_y(CurSliceAliased)')),4);
            end
        end
        
	end
    
    if(exist('csi_k_dummy','var'))
        csi_k = csi_k_dummy;
        clear csi_k_dummy csi_k_dummy2
    end

    
end

% Compute Mask of Slice Undersampled kPoints
MeasuredkPoints_SliceUndersampled = logical(squeeze(sum(abs(csi_k(:,:,:,:,1)),1))); 
size_csi = size(csi_k);




%% 3. InPlaneSimReco: Skip MRSI-Data
if(Par.Flags.TwoDCaipParallelImaging_flag)
    csi_k_skipped = feval(precision,Skip_kPoints(csi_k,Par.Settings.InPlaneCaipPattern));        % Skip Parts of the kSpace & convert to single if necessary

    % Compute the Mask of Measured Undersampled kPoints
    MeasuredkPoints_BothUndersampled = logical(squeeze(sum(abs(csi_k_skipped(:,:,:,:,1)),1))); 
    MeasuredkPoints_InPlaneUndersampled = repmat(MeasuredkPoints_BothUndersampled,[1 1 R_Slice]);
end
    
% Save first FID point of aliased images
if(SliceParallelImaging_SimulationFlag || TwoDCaipParallelImaging_SimulationFlag)
	csi_k_skipped_ForPlotting = csi_k_skipped(:,:,:,:,1,1);
	csi_skipped_ForPlotting = flipdim(fftshift(fftshift(conj(fft(fft(ifftshift(ifftshift(csi_k_skipped_ForPlotting,2),3),[],2),[],3)),2),3),2);
	csi_skipped_ForPlotting = sqrt(sum(abs(csi_skipped_ForPlotting).^2));
	csi_k_skipped_ForPlotting = csi_k_skipped_ForPlotting(1,:,:,:);
end


    
%% 4. InPlaneSimReco: Save the VD-Data

if(Par.Flags.TwoDCaipParallelImaging_flag)

    % Define kSpaceCenter_CSI
    kSpaceCenter_CSI = [floor(size_csi(2)/2)+1 floor(size_csi(3)/2)+1 0]; 

    % Create the VD Mask by Creating a Circle Around kSpace Center
    [x_grid,y_grid,z_grid] = ndgrid(1:size_csi(2), 1:size_csi(3),0);
    ExtraMeasuredkPoints_VD = (x_grid - kSpaceCenter_CSI(1)).^2 + (y_grid - kSpaceCenter_CSI(2)).^2 + (z_grid - kSpaceCenter_CSI(3)).^2 <= Par.Settings.VD_Radius^2;
    ExtraMeasuredkPoints_VD = repmat(ExtraMeasuredkPoints_VD, [1 1 size_csi(4)]); 

    % Compute mask of VD measured points
    MeasuredkPoints_VD = or(ExtraMeasuredkPoints_VD,MeasuredkPoints_BothUndersampled);

    % Save Fully Sampled csi_k Within VD-mask.  
    VD_ExtraData = csi_k(repmat(reshape(ExtraMeasuredkPoints_VD, [1 size_csi(2) size_csi(3) size_csi(4) 1]), [size_csi(1) 1 1 1 size_csi(5)]));

    clear csi_k kSpaceCenter_CSI x_grid y_grid z_grid


end





%% 5. InPlaneSimReco: InPlane Reconstruction

if(Par.Flags.TwoDCaipParallelImaging_flag)
    [csi_k_skipped,weights2DCaip,kernelsize2DCaip,SrcRelativeTarg2DCaip] = opencaipirinha_MRSI(csi_k_skipped,ACS,Par.Settings.InPlaneCaipPattern,false,precision,MinKernelSrcPts); 
    csi_k = csi_k_skipped;
    clear csi_k_skipped
end





%% 6. InPlaneSimReco: Elliptical Weight the Reconstructed Data
% The reconstruction assumes no elliptical weighting. So some kPoints outside the "Elliptical Weighting"-Circle were reconstructed and have meaningless values.

if(Par.Flags.TwoDCaipParallelImaging_flag)
    
    if(Par.CSI.Full_ElliptWeighted_Or_Weighted_Acq == 2)
        csi_k = EllipticalFilter(csi_k,[2 3], [1 1 1 Par.CSI.nFreqEnc_nonzf/2-1],1);
    end
    clear mask_FullySampled_cha cha

end



%% 7. InPlaneSimReco: Get the VD-Data Back

if(Par.Flags.TwoDCaipParallelImaging_flag)
    csi_k(repmat(reshape(ExtraMeasuredkPoints_VD, [1 size_csi(2) size_csi(3) size_csi(4) 1]), [size_csi(1) 1 1 1 size_csi(5)])) = VD_ExtraData;
    clear VD_ExtraData

else % These variables are needed later on. If TwoDCaip..._flag = 0, then the Slice_flag must be 1.
    
    MeasuredkPoints_InPlaneUndersampled = MeasuredkPoints_SliceUndersampled;
    MeasuredkPoints_BothUndersampled = MeasuredkPoints_SliceUndersampled;    
    MeasuredkPoints_VD = MeasuredkPoints_SliceUndersampled;
    
end



%% 8. Both: Compute the Reduction Factors

% Compute the InPlane Outer Reduction Factor
ORF_InPlane = sum(sum(MeasuredkPoints_FullySampled(:,:,1))) / sum(sum(MeasuredkPoints_InPlaneUndersampled(:,:,1)));
ORF_InPlane_NoEllip = sum(sum(MeasuredkPoints_FullySampled_NoEllip(:,:,1))) / sum(sum(MeasuredkPoints_InPlaneUndersampled(:,:,1)));

% The InPlane Reduction Factor R
R_InPlane = sum(sum(MeasuredkPoints_FullySampled(:,:,1))) / sum(sum(MeasuredkPoints_VD(:,:,1)));
R_InPlane_NoEllip = sum(sum(MeasuredkPoints_FullySampled_NoEllip(:,:,1))) / sum(sum(MeasuredkPoints_VD(:,:,1)));

% The Overall Reduction Factor
R = sum(MeasuredkPoints_FullySampled(:)) / sum(MeasuredkPoints_VD(:));
R_NoEllip = sum(MeasuredkPoints_FullySampled_NoEllip(:)) / sum(MeasuredkPoints_VD(:));



%% 9. SliceSimReco:  Perform Slice Unaliasing


if(Par.Flags.SliceParallelImaging_flag)    
	[csi_k,weights1DCaip] = openslicecaipirinha_MRSI(csi_k,ACS_Slice,cat(2,Par.Settings.FoVShifts_x',Par.Settings.FoVShifts_y'),false,precision); 
	size_csi = size(csi_k);    
end







%% 10. Both: Perform FFT

csi = feval(precision,csi_k);	% convert to single if necessary
csi_k = csi_k(1,:,:,:,1);	% Dont need that anymore, only for plotting
toc_sum = 0;

for cha = 1:size(csi,1)
    
    tic
    fprintf('\nFouriertransforming channel %02d\t...', cha)

    csi_channel = csi(cha,:,:,:,:,:);                           % This extra assignment proved to be faster than using always csi(cha,:,:,:,:,:). Seems that indexing is rather slow.
    
    csi_channel = ifftshift(ifftshift(csi_channel,2),3);
    csi_channel = fft(fft(csi_channel,[],2),[],3);
    csi_channel = conj(csi_channel);                            % the chem shift gets higher from right to left --> conj reverses that
    csi_channel = fftshift(fftshift(csi_channel,2),3);
    csi_channel = flipdim(csi_channel,2);                       % THIS FLIPS LEFT AND RIGHT IN SPATIAL DOMAIN BECAUSE PHYSICIANS WANT TO SEE IMAGES FLIPPED

    csi(cha,:,:,:,:,:) = csi_channel; 

    toc_sum = toc_sum + toc;
    fprintf('\ttook\t%10.6f seconds', toc)       

end

fprintf('\nOverall FFT Process\t\t...\ttook\t%10.6f seconds', toc_sum)
clear csi_channel toc_sum cha







%% 11. Both: Perform Same Stuff on Noise Data

if(exist('noise_sim','var'))
	noise_sim_PI = noise_sim;

    
    % Perform Slice Aliasing
    if(Par.Flags.SliceParallelImaging_flag)

        noise_sim_dummy = feval(precision,zeros([size(noise_sim,1) size(noise_sim,2) size(noise_sim,3) size(noise_sim,4)/R_Slice size(noise_sim,5)]));

        for AliLoop = 1:size(Par.Settings.SliceAliasingPattern,1)
            CurSliceAliased = Par.Settings.SliceAliasingPattern(AliLoop,:);
            % Alias Noise Data  (This has to be done in simulation, and in non-simulation, because also if the slices are measured aliased, the same is performed.)
            noise_sim_dummy(:,:,:,AliLoop,:) = sum(kSpace_FoVShift(noise_sim(:,:,:,CurSliceAliased,:),cat(2,Par.Settings.FoVShifts_x(CurSliceAliased)',Par.Settings.FoVShifts_y(CurSliceAliased)')),4);
        end
        
        noise_sim_PI = noise_sim_dummy;
        clear noise_sim_dummy
    end
    

    if(Par.Flags.TwoDCaipParallelImaging_flag)
        % Skip Parts of the kSpace; Has to be done in simulation AND in non-sim, because the noise data was always simulated
        noise_sim_PI_skipped = Skip_kPoints(noise_sim_PI,Par.Settings.InPlaneCaipPattern);    

        % Save VD kPoints
        VD_ExtraData_noise = noise_sim_PI(repmat(reshape(ExtraMeasuredkPoints_VD, [1 size(noise_sim_PI,2) size(noise_sim_PI,3) size(noise_sim_PI,4) 1]), [size(noise_sim_PI,1) 1 1 1 NoOfPseudoReplicas]));

        % Reconstruct data 
        noise_sim_PI = opencaipirinha_MRSI(noise_sim_PI_skipped,ACS,Par.Settings.InPlaneCaipPattern,true,precision,MinKernelSrcPts,weights2DCaip,kernelsize2DCaip,SrcRelativeTarg2DCaip);
        clear noise_sim_PI_skipped weights2DCaip kernelsize2DCaip SrcRelativeTarg2DCaip


        % Set to zeros those points that are not necessary (outside of circle)
        if(Par.CSI.Full_ElliptWeighted_Or_Weighted_Acq == 2)
            noise_sim_PI = EllipticalFilter(noise_sim_PI,[2 3], [1 1 1 size(noise_sim_PI,2)/2-1],1);
        end

        % Get VD Data back
        noise_sim_PI(repmat(reshape(ExtraMeasuredkPoints_VD, [1 size(noise_sim_PI,2) size(noise_sim_PI,3) size(noise_sim_PI,4) 1]), [size(noise_sim_PI,1) 1 1 1 NoOfPseudoReplicas])) = VD_ExtraData_noise;
    end


    
    % Perform Slice Unaliasing
    if(Par.Flags.SliceParallelImaging_flag)
        noise_sim_PI = openslicecaipirinha_MRSI(noise_sim_PI,weights1DCaip,cat(2,Par.Settings.FoVShifts_x',Par.Settings.FoVShifts_y'), true,precision); 
    end
    
    
    % FFT 
    noise_sim_PI = ifftshift(ifftshift(noise_sim_PI,2),3);
    noise_sim_PI = fft(fft(noise_sim_PI,[],2),[],3);
    noise_sim_PI = conj(noise_sim_PI);
    noise_sim_PI = fftshift(fftshift(noise_sim_PI,2),3);
    noise_sim_PI = flipdim(noise_sim_PI,2); 
    
end




%% 12. Both: Write Some Infos

mkdir([Par.Paths.out_path '/ParallelImaging']);

% Write a Text-File With Some Info
PI_fid = fopen(sprintf('%s/ParallelImaging/PI_Info.txt',Par.Paths.out_path),'w');
fprintf(PI_fid,'PARALLEL IMAGING INFO\n\n');
fprintf(PI_fid,'SliceParallelImaging_SimulationFlag = %d, TwoDCaipParallelImaging_SimulationFlag = %d\n',SliceParallelImaging_SimulationFlag,TwoDCaipParallelImaging_SimulationFlag);
fprintf(PI_fid,'In-Plane Outer Reduction Factor,\t\t\t\t\tORF_InPlane\t\t\t=\t%6.4f\n',ORF_InPlane);
fprintf(PI_fid,'In-Plane Reduction Factor,\t\t\t\t\t\t\tR_InPlane\t\t\t=\t%6.4f\n',R_InPlane);
fprintf(PI_fid,'In-Plane Outer Reduction Factor vs. fully Sampled,\tORF_InPlane_Fully\t=\t%6.4f\n',ORF_InPlane_NoEllip);
fprintf(PI_fid,'In-Plane Reduction Factor vs. fully Sampled,\t\tR_InPlane_Fully\t\t=\t%6.4f\n',R_InPlane_NoEllip);
fprintf(PI_fid,'Slice Reduction Factor,\t\t\t\t\t\t\t\tR_Slice\t\t\t\t=\t%6.4f\n',R_Slice);
fprintf(PI_fid,'Overall Reduction Factor,\t\t\t\t\t\t\tR\t\t\t\t\t=\t%6.4f\n',R);
fprintf(PI_fid,'Overall Reduction Factor vs. fully Sampled,\t\t\tR_Fully\t\t\t\t=\t%6.4f\n',R_NoEllip);

if(Par.Flags.TwoDCaipParallelImaging_flag); fprintf(PI_fid,'\nVD_Radius = %d\n',int16(Par.Settings.VD_Radius)); end

if(Par.Flags.SliceParallelImaging_flag)
	fprintf(PI_fid,'SliceAliasingPattern = \t'); 
	for i=1:size(Par.Settings.SliceAliasingPattern,1);
		fprintf(PI_fid,'Slicegroup %d: ',int16(i)); for j=1:size(Par.Settings.SliceAliasingPattern,2); fprintf(PI_fid,'%d ',int16(Par.Settings.SliceAliasingPattern(i,j))); end; 
		fprintf(PI_fid,'\n'); if(i<size(Par.Settings.SliceAliasingPattern,1)); fprintf(PI_fid,'\t\t\t\t\t\t'); end
	end

	fprintf(PI_fid,'FoVShifts_x = [ ');
	for i=1:numel(Par.Settings.FoVShifts_x); fprintf(PI_fid,'%f ',Par.Settings.FoVShifts_x(i)); end;
	fprintf(PI_fid,']\n');

	fprintf(PI_fid,'FoVShifts_y = [ ');
	for i=1:numel(Par.Settings.FoVShifts_y); fprintf(PI_fid,'%f ',Par.Settings.FoVShifts_y(i)); end;
	fprintf(PI_fid,']\n');
end

fprintf(PI_fid,'\nGeneral Settings:\n');
fprintf(PI_fid,'MinKernelSrcPts=%f\n',MinKernelSrcPts);
fprintf(PI_fid,'ACS_size=%f\n',ACS_size);

fclose(PI_fid);

% Data Itself
for Slc = 1:size(csi_k,4)
    
    Reco_k_fig = figure('visible','off');
    imagesc(squeeze(abs(csi_k(1,:,:,Slc,1))))
	colorbar;	
    saveas(Reco_k_fig,sprintf('%s/ParallelImaging/kSpaceReco_cha1_Slice%d', Par.Paths.out_path,Slc),'epsc2')
    close(Reco_k_fig)
    
    Reco_Image_fig = figure('visible','off');
    imagesc(squeeze(sqrt(sum(abs(csi(:,:,:,Slc,1)).^2))))
	colorbar;
    saveas(Reco_Image_fig,sprintf('%s/ParallelImaging/ImageReco_SoS_Slice%d', Par.Paths.out_path,Slc),'epsc2')
    close(Reco_Image_fig)
	
	if(Slc <= size(csi_k_meas,4))
		Measured_k_fig = figure('visible','off');
		imagesc(squeeze(abs(csi_k_meas(1,:,:,Slc,1))))
		colorbar;		
		saveas(Measured_k_fig,sprintf('%s/ParallelImaging/kSpaceMeas_cha1_Slice%d', Par.Paths.out_path,Slc),'epsc2')
		close(Measured_k_fig) 	

		Measured_Image_fig = figure('visible','off');
		imagesc(squeeze(sqrt(sum(abs(csi_meas(:,:,:,Slc,1)).^2))))
		colorbar;
		saveas(Measured_Image_fig,sprintf('%s/ParallelImaging/ImageMeas_SoS_Slice%d', Par.Paths.out_path,Slc),'epsc2')
		close(Measured_Image_fig)
	end
	
	
    
    if(SliceParallelImaging_SimulationFlag || TwoDCaipParallelImaging_SimulationFlag)
        
        Diff_Measured_Reco_k_fig = figure('visible','off');
        imagesc(squeeze(abs(csi_k_meas(1,:,:,Slc,1) - csi_k(1,:,:,Slc,1))))
		colorbar;
        saveas(Diff_Measured_Reco_k_fig,sprintf('%s/ParallelImaging/kSpace_Meas_-_Reco_cha1_Slice%d', Par.Paths.out_path,Slc),'epsc2')	
        close(Diff_Measured_Reco_k_fig)
        
        Diff_Measured_Reco_Image_fig = figure('visible','off');
        imagesc(squeeze(sqrt(sum(abs(csi_meas(:,:,:,Slc,1) - csi(:,:,:,Slc,1)).^2))))
		colorbar;
        saveas(Diff_Measured_Reco_Image_fig,sprintf('%s/ParallelImaging/Image_Meas_-_Reco_SoS_Slice%d', Par.Paths.out_path,Slc),'epsc2')
		close(Diff_Measured_Reco_Image_fig)
    
		if(Slc <= size(csi_skipped_ForPlotting,4))
			csi_skipped_ForPlotting_fig = figure('visible','off');
			imagesc(squeeze(abs(csi_skipped_ForPlotting(1,:,:,Slc,1) )))
			colorbar;
			saveas(csi_skipped_ForPlotting_fig,sprintf('%s/ParallelImaging/Image_Undersampled_SoS_Slice%d', Par.Paths.out_path,Slc),'epsc2')	
			close(csi_skipped_ForPlotting_fig)

			csi_k_skipped_ForPlotting_fig = figure('visible','off');
			imagesc(squeeze(abs(csi_k_skipped_ForPlotting(1,:,:,Slc,1) )))
			colorbar;
			saveas(csi_k_skipped_ForPlotting_fig,sprintf('%s/ParallelImaging/kSpace_Undersampled_cha1_Slice%d', Par.Paths.out_path,Slc),'epsc2')	
			close(csi_k_skipped_ForPlotting_fig)
		end
		
		
    end
    
end

% kSpace Patterns
if(Par.Flags.TwoDCaipParallelImaging_flag)
    
    UndersamplingPattern_fig = figure('visible','off');
    imagesc(Par.Settings.InPlaneCaipPattern)
    colormap(gray)
    saveas(UndersamplingPattern_fig,sprintf('%s/ParallelImaging/UndersamplingPattern', Par.Paths.out_path),'epsc2')
    close(UndersamplingPattern_fig)
	
	Measured_fig = figure('visible','off');
	imagesc(squeeze(logical(abs(csi_k_meas(1,:,:,1,1)))))
	colormap(gray)
	saveas(Measured_fig,sprintf('%s/ParallelImaging/MeasuredMask', Par.Paths.out_path),'epsc2')
	close(Measured_fig)	
	
	VDSampled_fig = figure('visible','off');
    imagesc(MeasuredkPoints_VD(:,:,1))
    colormap(gray)
    saveas(VDSampled_fig,sprintf('%s/ParallelImaging/VDSampledMask', Par.Paths.out_path),'epsc2')
    close(VDSampled_fig)
	
	UnderSampled_fig = figure('visible','off');
	imagesc(MeasuredkPoints_BothUndersampled(:,:,1))
	colormap(gray)
	saveas(UnderSampled_fig,sprintf('%s/ParallelImaging/UnderSampledMask', Par.Paths.out_path),'epsc2')
	close(UnderSampled_fig)
    
	if(~TwoDCaipParallelImaging_SimulationFlag)			% In the simulated case this is identical to the MeasuredMask
		FullySampled_fig = figure('visible','off');
		imagesc(MeasuredkPoints_FullySampled(:,:,1))
		colormap(gray)
		saveas(FullySampled_fig,sprintf('%s/ParallelImaging/FullySampledMask', Par.Paths.out_path),'epsc2')
		close(FullySampled_fig)
	end
	
end

clear PI_fid UndersamplingPattern_fig FullySampled_fig UnderSampled_fig VDSampled_fig Slc Measured_Image_fig Reco_k_fig Reco_Image_fig Measured_k_fig Diff_Measured_Reco_k_fig Diff_Measured_Reco_Image_fig
clear csi_k_meas csi_meas








%% +. Postparations

clear R_x R_y ORF_InPlane InPlaneCaipPattern VD_Radius ACS_size TwoDCaipParallelImaging_SimulationFlag ACS MeasuredkPoints_FullySampled MeasuredkPoints_FullySampled_NoEllip
clear MeasuredkPoints_UnderSampled MeasuredkPoints_VD Skip_kPoints_Spatial ExtraMeasuredkPoints_VD
clear ExtraMeasuredkPoints_VD VD_ExtraData_noise
clear weights2DCaip kernelsize2DCaip SrcRelativeTarg2DCaip weights1DCaip
