%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%                           PROGRAM TO HELL AND BACK                         %%%%%%%%%%%%%%
%%%%%%%%%%%%%% ENTSTANDEN DURCH DIE KUNST DES PROGRAMMIERENS DURCH KONSEQUENTES ANSTARREN %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -except tmp_dir
load([tmp_dir '/Parameters.mat'])
Par.CSI.total_channel_no;
global totCoili;

global NoOfTempInterleaves;
global NoiseCorrMat_post;
FastPIReprocess_flag = 0;



if(~exist('CurAvg','var'))		% For those cases someone wants to run this script manually
    CurAvg = 1;
end
CurAvg;
CoilCombtic = tic;
%pause on

%% 1. Load Parameters

if (exist([tmp_dir '/Parameters_water.mat'],'file'))
    load([tmp_dir '/Parameters_water.mat'])
    Par.Paths.csi_path{CurAvg} = Par.Paths.WaterReference_path;
else
    load([tmp_dir '/Parameters.mat'])
end


Noise_CorrMat = 1;


%% 3a. Read in CSI file

% Hack
Concept_flag = true;
display([char(10) char(10) 'READ CSI DATA'])    % char(10) = newline

SettingsInfo.ONLINE.NoFFT_flag = true; 
[Noise,SettingsInfo,ReadInInfo,kSpace] = read_csi_new(Par.Paths.csi_path{CurAvg},tmp_dir,SettingsInfo);
csi_k = feval('single',kSpace.ONLINE); kSpace.ONLINE = 0;
size_csi = size(csi_k);	

image_k = feval('single',kSpace.PATREFSCAN);
image_k = reshape(image_k, [size(image_k,1) size(image_k,2) size(image_k,3) size(image_k,4)*size(image_k,5) size(image_k,6)]);		           

size_csi = [size_csi(1) size_csi(2) size_csi(3) size_csi(4)*size_csi(5) size_csi(6)];
csi_k = reshape(csi_k, [size(csi_k,1) size(csi_k,2) size(csi_k,3) size(csi_k,4)*size(csi_k,5) size(csi_k,6)]);		           


% Get Prescan Noise
if(isfield(Noise,'ONLINE') && numel(Noise.ONLINE) > 1)
    Noise_mat = Noise.ONLINE;
end
Noise_CorrMat = SettingsInfo.ONLINE.NoiseCorrMat;



clear Noise iSpace kSpace;


size_csi=size(csi_k); 
totCoili=size(csi_k,1);


%% 3c Perform Noise Decorrelation
 if(exist('NoiseCorrMat_post','var') && totCoili >1 && Par.Flags.mask_flag)
     fprintf('\nPerform Noise Decorrelation!.')
     if(size(csi_k,2) > size(csi_k,3))
         for kreiserl = 1:size(csi_k,2) %better dont change this index, because this is the most memory efficient looping
             csi_k(:,kreiserl,:,:,:) = PerformNoiseDecorrelation(csi_k(:,kreiserl,:,:,:),NoiseCorrMat_post);    
             image_k(:,kreiserl,:,:,:) = PerformNoiseDecorrelation(image_k(:,kreiserl,:,:,:),NoiseCorrMat_post);    
         end
     else
         for kreiserl = 1:size(csi_k,3) %bLipidDecon_MethodAndNoOfLoopsetter dont change this index, because this is the most memory efficient looping
             csi_k(:,:,kreiserl,:,:) = PerformNoiseDecorrelation(csi_k(:,:,kreiserl,:,:),NoiseCorrMat_post);    
             image_k(:,:,kreiserl,:,:) = PerformNoiseDecorrelation(image_k(:,:,kreiserl,:,:),NoiseCorrMat_post);    
         end
     end
 end



%% 4. READ IN IMAGING-DATA, IMAGING-FLIP-DATA, IMAGING-VC-DATA, MASK-DATA, PREPARE DATA

if(~FastPIReprocess_flag)
	display([char(10) char(10) 'READ IMAGING DATA AND MASK'])

    % 	MASK
	if(exist([tmp_dir '/mask_brain.raw'],'file'))
		fid_mask = fopen([tmp_dir '/mask_brain.raw'],'r');
		mask = reshape(fread(fid_mask, 'float'), [ReadInInfo.Par.ReadLines ReadInInfo.Par.PhaseLines ReadInInfo.Par.SliceLines]);
		fclose(fid_mask);
		
		if(Par.Flags.InterpolateCSIResolution_flag == 1)
			fid_mask = fopen([tmp_dir '/mask_brain_BefInterpol.raw'],'r');
			mask_BefInterpol = reshape(fread(fid_mask, 'float'), [Par.CSI.nFreqEnc_BefInterpol Par.CSI.nPhasEnc_BefInterpol Par.CSI.nSLC_BefInterpol*Par.CSI.nPartEnc_BefInterpol]);
			fclose(fid_mask);
		else
			mask_BefInterpol = mask;
		end
		
		if(sum(sum(sum(mask))) == 0)
			mask = ones([ReadInInfo.Par.ReadLines ReadInInfo.Par.ReadLines ReadInInfo.Par.SliceLines]);
		end
	else
		mask = ones([ReadInInfo.Par.ReadLines ReadInInfo.Par.ReadLines ReadInInfo.Par.SliceLines]);
	end

end


%% 5. Simulate noise
%Since the data is already prewhitened, the noise can be simulated as random numbers with std = 1.


if(numel(Noise_CorrMat) > 1 || (numel(Noise_CorrMat) == 1 && Noise_CorrMat ~= 1 && Noise_CorrMat ~= 0))		% For VC data Noise_CorrMat holds only one element,
    display([char(10) char(10) 'Simulate Noise For SNR Calculation.'])										% but still we might want to calculate the SNR...
    NoOfPseudoReplicas = 200;
    % noise_sim = complex(randn([size_csi(1) size(mask_BefInterpol,1) size(mask_BefInterpol,2) size(mask_BefInterpol,3) NoOfPseudoReplicas]), ...
    % randn([size_csi(1) size(mask_BefInterpol,1) size(mask_BefInterpol,2) size(mask_BefInterpol,3) NoOfPseudoReplicas]));
    size_csi_ForNoise = size_csi; size_csi_ForNoise(5) = NoOfPseudoReplicas; size_csi_ForNoise(6) = ReadInInfo.Par.nAve;
    noise_sim = complex(randn(size_csi_ForNoise), randn(size_csi_ForNoise));
	noise_sim = feval('single',noise_sim);		% convert ot single if necessary
    if(Par.Flags.zerofill_to_nextpow2_flag || numel(strfind(Par.Paths.csi_path{CurAvg}, '.IMA')) > 0)
        noise_sim = ZerofillOrCutkSpace(noise_sim,[size_csi(1) size(mask_BefInterpol,1) size(mask_BefInterpol,2) size(mask_BefInterpol,3) NoOfPseudoReplicas],1);        
    end
    
    % Sum Averages
    if(size_csi_ForNoise(6) > 1)
        noise_sim = sum(noise_sim,6)/size(noise_sim,6);
    end
    
    % Hadamard Decoding
    if(size_csi(4) > 1 && ~Par.CSI.ThreeD_flag)
        display('Perform Hadamard Decoding')
        noise_sim = hadamard_decoding(noise_sim,4);
    end
end

isStackofRings = 0;
if(ReadInInfo.Par.SliceLines*ReadInInfo.Par.PhaseLines/2 == size(csi_k,2))
    isStackofRings = 1;
end
% display([char(10) char(10) 'PERFORM PI SIM'])


%% 7. Perform Regridding for Concept
if(Concept_flag)

    ReadInInfo.Par.SliceNormalVector_x = Par.CSI.SliceNormalVector_x;   %<<<<------ these 3 lines
    ReadInInfo.Par.SliceNormalVector_y = Par.CSI.SliceNormalVector_y;
    ReadInInfo.Par.SliceNormalVector_z = Par.CSI.SliceNormalVector_z;   
    ReadInInfo.Par.InPlaneRotation=Par.CSI.InPlaneRotation;
    nSlices=ReadInInfo.Par.SliceLines;

    image_new=zeros([size(image_k,1) ReadInInfo.Par.ReadLines ReadInInfo.Par.PhaseLines nSlices 1],'single');	   

    totCoili_parfor=totCoili;
    mask_flag_parfor=Par.Flags.mask_flag;

    tic
    %%B0MAP (out of musical scann)
%     csi_k=abs(csi_k).*exp(1i*angle(csi_k)).*exp(1i*repmat(angle(image_k(:,:,:,:,4)-image_k(:,:,:,:,3)),[1 1 1 1 size(csi_k,5)]));
    
   % if (matlabpool('size')~=12) 
   %     matlabpool(12)
   % end 
    display([char(10) char(10) 'PERFORM 2D GRIDDING FOR CONCEPT'])    % char(10) = newline 
    for slicenmb=1:nSlices
        [csi_k_p,image_k_p,noise_sim_k_p,NumberOfLoopsPerSlice,NoOfTempInterleaves_vec]=Initialize3DRollercoaster(slicenmb,csi_k,image_k,noise_sim,nSlices,ReadInInfo,isStackofRings);

%         NoOfTempInterleaves_parfor=unique(ReadInInfo.Data(5,:)); for VB17  old
%          NoOfTempInterleaves_vec=repmat(NoOfTempInterleaves_parfor,[size(csi_k_p,2) 1]);
          
        nc=max(NumberOfLoopsPerSlice); % this will be overwritten in the reco soon
          
        csi_k_backup=csi_k_p;         %for rollercoaster
        conj_flag = false;       %also the flag for no-musical
        post_grd_dc_flag = true;    %Do post grd weight calculation during musical: samples=1.
        [~,true_weights]=RollercoasterRecon_fast(totCoili_parfor,slicenmb,NumberOfLoopsPerSlice,ReadInInfo,nSlices,post_grd_dc_flag,conj_flag,NoOfTempInterleaves_vec,csi_k_p); % POST-gridding weights (true_weights)
        post_grd_dc_flag = false;
        [output_k,~]=RollercoasterRecon_fast(totCoili_parfor,slicenmb,NumberOfLoopsPerSlice,ReadInInfo,nSlices,post_grd_dc_flag,conj_flag,NoOfTempInterleaves_vec,image_k_p,true_weights);%apply post grd weights here
        image_new(:,:,:,slicenmb,:)=output_k;
        conj_flag = true;
        if(totCoili_parfor>1 && mask_flag_parfor)
            [output_k,~]=RollercoasterRecon_fast(totCoili_parfor,slicenmb,NumberOfLoopsPerSlice,ReadInInfo,nSlices,post_grd_dc_flag,conj_flag,NoOfTempInterleaves_vec,noise_sim_k_p,true_weights); %NOISE
            noise_sim_new(:,:,:,:,:,slicenmb)=output_k;
        end
        csi_k_p=csi_k_backup; %this is for rollercoasting now and has output csi
        [output_k,~]=RollercoasterRecon_fast(totCoili_parfor,slicenmb,NumberOfLoopsPerSlice,ReadInInfo,nSlices,post_grd_dc_flag,conj_flag,NoOfTempInterleaves_vec,csi_k_p,true_weights);%apply post grd weights here
        csi(:,:,:,:,:,slicenmb)=output_k;
    end
  %  clear output_k
  % if(totCoili_parfor>1 && mask_flag_parfor && matlabpool('size')==12)
  %      noise_sim_new=permute(noise_sim_new, [1 2 3 6 5 4]);
   % end
   % if(matlabpool('size')==12)%parfor makes some strange permutations
   %     csi=permute(csi, [1 2 3 6 5 4]);
   % end
   % toc

 	clear output_k
        if(totCoili_parfor>1 && mask_flag_parfor)
            noise_sim_new=permute(noise_sim_new, [1 2 3 6 5 4]);
        end

            csi=permute(csi, [1 2 3 6 5 4]);                         %<<<<----- no if (matlabpool) etc..

        toc 
    
    clear csi_k image_k noise_sim_k_p
    
    if(totCoili>1 && Par.Flags.mask_flag)
    noise_sim=noise_sim_new;
    end
    image=image_new(:,:,:,:,:);

    clear csi_new noise_sim_new image_new
    if(nSlices>1)
        %hamming in z
%              dummyfilter=chebwin(size(csi,4));
        dummyfilter=HammingFilter(ones(size(csi,4),1),[1 2]);
        for i=1:size(csi,4)
            csi(:,:,:,i,:)=csi(:,:,:,i,:)*dummyfilter(i,1);
            if(totCoili>1 && Par.Flags.mask_flag)
            noise_sim(:,:,:,i,:)=noise_sim(:,:,:,i,:)*dummyfilter(i,1);
            end
            image(:,:,:,i,:)=image(:,:,:,i,:)*dummyfilter(i,1);
        end  


        csi = ifftshift(csi,4); %Note that ifftshift(fft(fftshift is not the same as  fftshift(fft(ifftshift for odd N
        csi = fft(csi,[],4);
        csi = fftshift(csi,4);
        if(totCoili>1 && Par.Flags.mask_flag)
        noise_sim = ifftshift(noise_sim,4);
        noise_sim = fft(noise_sim,[],4);
        noise_sim = fftshift(noise_sim,4);
        end
        image = ifftshift(image,4);
        image = fft(image,[],4);
        image = fftshift(image,4);
    end

     csi=conj(csi);            % the chem shift gets higher from right to left --> conj reverses that
     csi = flipdim(csi,4);     % THIS FLIPS LEFT AND RIGHT IN SPATIAL DOMAIN BECAUSE PHYSICIANS WANT TO SEE IMAGES FLIPPED
     image = flipdim(image,4);
     if(totCoili>1 && Par.Flags.mask_flag)
     noise_sim=conj(noise_sim);   
     noise_sim = flipdim(noise_sim,4);
     end
end

size_csi = size(csi);


%% 7. Bilgic Lipid Decontamination

if(Par.Flags.LipidDecon_flag == 1 && ~exist([tmp_dir '/Parameters_water.mat'],'file'))
    
	fprintf('\n\nPerform Bilgic Lipid Decontamination: channelwise & sensitivity-weighted\n')
	
    mkdir([Par.Paths.out_path '/LipidDecontamination'])
    mkdir([Par.Paths.out_path '/LipidMaskNii'])
   
    x_dim=size(csi,2);
    y_dim=size(csi,3);
    z_dim=size(csi,4);
     
    lipid_mask_total=zeros(x_dim,y_dim,z_dim);
    
    	f_range=size(csi,5);
	f_range_start=ceil(f_range*0.82*(1-(180000/(ReadInInfo.Par.Dwelltime)-1)/2));
	f_range_end=ceil(f_range*0.96*(1-(180000/(ReadInInfo.Par.Dwelltime)-1)/2)); 
        
    Text=sprintf('\n\nVectorsize Points of the csi go from 0 to %d!',f_range);
    disp(Text)
    Text=sprintf('\nRemoving Lipids in the Range from %d to %d:\n',f_range_start,f_range_end);
    disp(Text)
    
    for nCha = 1:size(csi,1)     
        
        csi_temp = squeeze(csi(nCha,:,:,:,:)); 
        csi_lip = fftshift(fft(squeeze(csi_temp(:,:,:,:)),[],4),4);
        
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
            imagesc(squeeze(abs(csi(nCha,:,:,Slc))),[0 max_sig])
            colorbar;      
            saveas(Brain_fig,sprintf('%s/LipidDecontamination/Channel_%d/brain_mask_slc_%d', Par.Paths.out_path,nCha,Slc),'epsc2')
            saveas(Brain_fig,sprintf('%s/LipidDecontamination/Channel_%d/brain_mask_slc_%d', Par.Paths.out_path,nCha,Slc),'fig')
            close(Brain_fig)
        end
      
        fprintf('Decontaminating channel %d\n', nCha)    
        
        for Slc = 1:size(csi_temp,3)

            % Prepare slices: Get slice csi data, Extract Relevant Lipids, Get slice mask
            csi_slc = fftshift(fft(squeeze(csi_temp(:,:,Slc,:)),[],3),3);

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
            csi(nCha,:,:,Slc,:) = ifft(fftshift(csi_slc,3),[],3);
            
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
    
end  
%% 8. Compute Weights w_n


if(~FastPIReprocess_flag)

	display([ char(10) char(10) 'Compute the Coil Combination Weights w_n'])



	if(size_csi(1) > 1 && exist('image','var') && Par.Flags.image_VC_flag)                  % Sensmap Weights
		display([char(10) 'Compute a sensitivity map and use it for coil combination . . .' char(10)])       
		csi_size = size(csi);
		image_resized = zeros(csi_size(1:4));
		image_VC_resized = zeros([1 csi_size(2:4)]);
		for slice_no = 1: size(image,4)
			for channel_no = 1:size(image,1)
				image_resized(channel_no,:,:,slice_no) = imresize(imresize(squeeze(image(channel_no,:,:,slice_no)), [size(image,2)/4 size(image,3)/4]), [size_csi(2) size_csi(3)]);
			end
			image_VC_resized(1,:,:,slice_no) = imresize(imresize(squeeze(image_VC(1,:,:,slice_no)), [size(image_VC,2)/4 size(image_VC,3)/4]), [size_csi(2) size_csi(3)]);
		end
		weights = image_resized ./ repmat(image_VC_resized, [32 1 1 1]);




	elseif(totCoili > 1 && exist('image','var'))									% MUSICAL weights
		display([char(10) 'Use the imaging data for coil combination.' char(10)])   
		weights = image;                                                            % No need for conjugating, because the csi data is already conjugated 
																					% (Siemens does that, so that the spectra increase in frequency from right to left)



	elseif(size_csi(1) > 1 && ~exist('image','var'))                              % 1st FID point weights
		display([char(10) 'Use the first FID point for coil combination.' char(10)])
		weights = conj(csi(:,:,:,:,1));




	elseif(totCoili == 1 && exist('image','var'))                              % Weights and for VC data if image was inputted.
		display([char(10) 'Phase csi with the imaging data.' char(10)])    
		weights = image;                                              % Only phase VC data. This makes abs(weights) = 1. Thus also scaling = 1. Hence csi data is only phased.

	end

	clear image image_VC

end

% Load the water reference weights and overwrite current ones, because: If
% imaging weights are used
if(Par.Flags.WaterReference_flag && exist([Par.Paths.out_path '/WaterReference.mat' ],'file') && exist('weights','var'))	% If water ref alrdy exists, load and use its weights
	load([Par.Paths.out_path '/WaterReference.mat' ],'weights')
end



%% 10. Weight and Sum CSI Data
display([ char(10) char(10) 'Weight and Sum CSI Data'])
csi = feval('single',csi);		% convert to single if necessary
if(exist('weights','var'))
	
	% To be more memory efficient...
	for t = 1:size_csi(5)
		csi_t = csi(:,:,:,:,t);
		csi_t = csi_t .* weights;
		csi(:,:,:,:,t) = csi_t;
	end
	csi = sum(csi,1);    
    
end


% reshape csi
csi = reshape(squeeze(csi),[size_csi(2) size_csi(3) size_csi(4) size_csi(5)]);



%% 10. Compute Inverse Scaling Matrix 

if(~FastPIReprocess_flag)
    
    display([ char(10) char(10) 'Compute the Scaling Matrix . . .'])
    % scaling-matrix. Has only effect on displaying the result (e.g. if for inhomogeneous reception profile is corrected for)
    
    if(exist('weights','var'))
        if( Par.Flags.WaterReference_flag ) % If water reference is used, use the same scaling for both csi and water dataset
            %scaling_inv = squeeze_single_dim(sqrt(sum(abs(weights).^2,1)),1);
	    scaling_inv = squeeze_single_dim(sum(abs(weights).^2,1),1);
        else
            scaling_inv = squeeze_single_dim(sum(abs(weights).^2,1),1);
        end
        scaling_inv = repmat(scaling_inv/(10^5), [1 1 1 size_csi(5)]);       % 10^5: To make concentrations in a nicer range.
        
        % Print figures
        for Slc = 1:size(scaling_inv,3)
            Scaling_fig = figure('visible','off');
            imagesc(1./squeeze(scaling_inv(:,:,Slc,1)),[0 5/scaling_inv(ceil(size(scaling_inv,1)/2),ceil(size(scaling_inv,2)/2),1,1) ])
            colorbar;
            saveas(Scaling_fig,sprintf('%s/scalings/scaling_slc%d', Par.Paths.out_path,Slc),'epsc2')
            saveas(Scaling_fig,sprintf('%s/scalings/scaling_slc%d', Par.Paths.out_path,Slc),'fig')
            close(Scaling_fig)
        end
        
    end
end




%% 11. Scale MRSI DATA
display([ char(10) char(10) 'Scale MRSI Data'])

% B1mapi=dicomread('/net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Measurement_Data/ComparisonRollerVsSpiral/AbstractMeasNew/ConceptB1/b1/B1CONCEPT_TEST.MR.PHYSIKER_PMOS.0007.0001.2017.09.05.13.59.59.359375.95736104.IMA');
% B1mapi=B1mapi';
% B1mapi=im2double(B1mapi);
% scaling_inv=repmat(B1mapi,[1 1 1 size(scaling_inv,4)]);
% 
if(exist('scaling_inv','var'))
                                                % 0.1492: To make a NAA peak that has amplitude 1 in MATLAB have the concentration of 1 in LCModel.
	csi = csi ./ (scaling_inv);          % Additional scaling of csi by 0.1492. Noise is not scaled by this factor, because the noise has std 1 in MATLAB. But the noise is not
                                         
end






%% 13. Perform Interpolation

% Interpolation Image Domain
if(Par.Flags.InterpolateCSIResolution_flag)
	csi = InterpolateCSIImages(csi,0,mask_BefInterpol,Par.Settings.InterpolateResolutionRatio);
end

% Interpolation Time Domain
% In time domain __________________________________________________________________________________
if(Par.Flags.TimeInterpolation_flag)
	display([ char(10) 'PERFORM INTERPOLATION IN TIME-DOMAIN' char(10)])
	csi = TimeInterpolation(csi,Par.Settings.TimeInterpolationFactor(1),Par.Settings.TimeInterpolationFactor(2),Par.Settings.TimeInterpolationFactor(3));
	display(['TruncateFactor = ', num2str(Par.Settings.TimeInterpolationFactor(1)), ', ZerofillFactor = ', num2str(Par.Settings.TimeInterpolationFactor(2)), '.'])
	if(Par.Settings.TimeInterpolationFactor(3) == 1)
	display(['Filling to original size.'])	
	end

end
		display(['Done.']) 		


%% 13. Apply exponential Time-Domain Filter to Combined MRSI Data

if (Par.Flags.exponential_filter_Hz_flag)
    display([ char(10) 'APPLY EXP TIME-DOMAIN FILTER' char(10)])   
	csi = ExponentialFilter(csi, Par.CSI.Dwelltimes,Par.Settings.exponential_filter_Hz,4);
end





%% 14. Perform same things on noise as on CSI
% Do that for the non-PI and the PI noise seperately, in order to calculate the g-factor map


if(exist('noise_sim','var'))
    display([ char(10) char(10) 'Perform Same Calculations on Simulated Noise for SNR Calculation'])   


    % Parallel Imaging Performed in ParallelImagingSimReco.m


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%   NON-PI NOISE   %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Spatial FFT
%     noise_sim = ifftshift(ifftshift(noise_sim,2),3);
%     noise_sim = fft(fft(noise_sim,[],2),[],3);
%     if(Par.CSI.ThreeD_flag && size(noise_sim,3) > 1)
%         noise_sim = fft(ifftshift(noise_sim,4),[],4);
%     end
%     noise_sim = conj(noise_sim);
%     noise_sim = fftshift(fftshift(noise_sim,2),3);
%     if(Par.CSI.ThreeD_flag && size(noise_sim,3) > 1)
%         noise_sim = fftshift(noise_sim,4);
%     end
%     noise_sim = flipdim(noise_sim,2);

	
	
	
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
			CurNoiseData = CurNoiseData ./ repmat(scaling_inv(:,:,:,1),[1 1 1 size(CurNoiseData,4)]); 
		end
		
		
		% Lipid Decon (Can this be simulated?)

		
		
		% Frequency Alignment
		if(Par.Flags.AlignFrequency_flag)
			InData.csi = CurNoiseData; InData.mask = mask_BefInterpol;
			CurNoiseData = FrequencyAlignment(InData,ShiftMap,4,2);
		end

		
		
		% Hamming Filtering
		if(Par.Flags.hamming_flag)
			if(size(CurNoiseData,3) > 1 && Par.CSI.ThreeD_flag)
				CurNoiseData = HammingFilter(CurNoiseData,[1 2 3],Par.Settings.hamming_factor,'OuterProduct',0);    
			else
				CurNoiseData = HammingFilter(CurNoiseData,[1 2],Par.Settings.hamming_factor,'OuterProduct',0);
			end 
		end
		
		
		% Interpolation Image Domain
		if(Par.Flags.InterpolateCSIResolution_flag)
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
        noise_sim_spectral = fftshift(fft(noise_sim_spectral,[],5),5) / sqrt(2*size(noise_sim_spectral,5)); % Why is this sqrt(2) necessary?  BUG: 5 --> 4
        noise_sim_time = noise_sim_PI * 0.1370;                                                             % But without it, the spectrum is differently scaled than the LCModel spectrum.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%    NON-PI NOISE    %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	else
        % In case of no PI: SNR calculated based on normal noise
        noise_sim_spectral = noise_sim;
        size(noise_sim_spectral)
        noise_sim_spectral = fftshift(fft(noise_sim_spectral,[],4),4) / sqrt(2*size(noise_sim_spectral,4)); % BUG: 5 --> 4
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




%% 15. Correct for 1storderphase
%     Rationale: Only do 1stPhCorr if: - FID sequence (--> TE < 10000 us = 10ms) and - The user asks for it or it is the water reference scan
    if( (Par.Flags.FirstOrderPhaseCorr_flag || (strcmp(Par.Settings.WaterReference_method, 'W1') && ~exist([Par.Paths.out_path '/WaterReference.mat' ],'file'))) && min(Par.CSI.TEs) < 10000 )
           size_csi=size(csi);

           omega_step = 2*pi/((size_csi(4)-1)*ReadInInfo.Par.Dwelltime*2/10^9);
           if(~mod(size(csi,4),2))
               fprintf('\n\nbu')
               FreqVec = repmat(-ceil(size_csi(4)/2) : (ceil(size_csi(4)/2)-1), [numel(Par.CSI.TEs) 1]);
           else
               fprintf('\n\nba')
               FreqVec = repmat(-ceil((size_csi(4)-1)/2) : (ceil((size_csi(4)-1)/2)), [numel(Par.CSI.TEs) 1]);
           end
            TEs = repmat(Par.CSI.TEs/1000, [1 size(FreqVec,2)]);
            FirstOrderCorrVec = exp(-1i * (TEs/1000)*omega_step .* FreqVec);
            clear TEs FreqVec
            FirstOrderCorrVec = reshape(FirstOrderCorrVec, [1 1 size(FirstOrderCorrVec,1) size(FirstOrderCorrVec,2)]);

            csi_spec = fftshift(fft(csi,[],4),4);
            csi_spec = csi_spec .* repmat(FirstOrderCorrVec, [size_csi(1) size_csi(2) size_csi(3) 1]);
            csi = ifft(ifftshift(csi_spec,4),[],4);

    end



%% 16. Manipulate 1storderphase for AD-"change" (GH 2019)
if(Par.Flags.FirstOrderPhaseModulation_flag)
	display([ char(10) char(10) 'Truncate FID to simulate different AD'])
	Cut_time=Par.Settings.FID_Truncation_in_ms   %should be in ms here
	ReadInInfo.Par.Dwelltime
	Cut_length=floor(Cut_time/(ReadInInfo.Par.Dwelltime*2/10^6)) %just calc it as ms %*2 bc displayed BW is half of the right value (verified vy manual calc)
	Cut_length=Cut_length+1; %account for starting at 1

    %csi_spec = fftshift(fft(csi,[],4),4);
    csi = csi(:,:,:,Cut_length:size(csi,4));
    %csi = ifft(ifftshift(csi_spec,4),[],4);
	size(csi,4)
end



%% 17. Handle Water Reference

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



%% 18. write the files necessary for LCM-processing

% Dont do for W1-WatRef
if( IsW2WatRef || (~IsWatRef) )  % i.e. either we have W2-Watref, or no Watref and last Avg

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
	Paths.basis_file = Par.Paths.basis_path;
	Paths.LCM_ProgramPath = Par.Paths.LCM_Path;
	Paths.batchdir = tmp_dir;
	MetaInfo.DimNames = {'x','y','z'};
	MetaInfo.pat_name = Par.CSI.PatName;
	MetaInfo.LarmorFreq =  ReadInInfo.Par.LarmorFreq;
	MetaInfo.dwelltime =  ReadInInfo.Par.Dwelltime*2;
    
    csi=circshift(csi,[1 -1 0 0]);
    csi=Zerofilling_Spectral(csi,[size(csi,1) size(csi,2) size(csi,3) size(csi,4)*2]);
% 	mask=zeros(size(mask));
%     mask(31,23,25)=1;
    
    Data.csi = csi;
	Write_LCM_files(Data,Paths,MetaInfo,ControlInfo,mask,8) 

end


%% 17. Save the processed MRSI data for SNR-Computation after LCModel processing

fprintf('\n\nThe MRSI Pre-Processing and LCModel preparations took %10.6f s.\n',toc(CoilCombtic))
clearvars -except Par csi weights mask g_FactorMap NoiseScalingMatrix_spectral* NoiseScalingMatrix_time IsWatRef NoiseCorrMat_post ReadInInfo

if(1)		% No need to save water reference
	save([Par.Paths.out_path '/CombinedCSI.mat'],'-v7.3')
end

