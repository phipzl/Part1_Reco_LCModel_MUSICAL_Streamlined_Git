function [Noise,PreProcessingInfo, kSpace] = PreProcessMRIData_Wrapper(kSpace,PreProcessingInfo,ReadInInfo)
%
% read_csi_dat Read in csi-data from Siemens raw file format
%
% This function was written by Bernhard Strasser, July 2012.
%
%
% The function can read in MRS(I) data in the Siemens raw file format ".DAT" and performs
% some easy Postprocessing steps like zerofilling, Hadamard decoding, Noise Decorrelation etc.
%
%
% [iSpace,Noise,PreProcessingInfo, kSpace] = PreProcessMRIData_Wrapper(kSpace,PreProcessingInfo,ReadInInfo)
%
% Input: 
% -         kSpace                      ...     Path of MRS(I) file.
% -         PreProcessingInfo           ...     Info about how the read in data should be pre-processed.
%                                               Sub-Fields: For each subdataset (EvalInfoMask entry) one, e.g. 'ONLINE', 'PATREFANDIMASCAN', and 'NOISEADJSCAN'
%                                               Sub-Sub-Fields: Each Sub-Field can have the following entries:
%                                               - NoFFT_flag: If true, do no fft is performed from kSpace to iSpace. Example: PreProcessingInfo.ONLINE.NoFFT_flag = true;
%                                               - fredir_shift: Corrects for a different phase caused by a shift in the frequency encoding direction
%                                               - FlipkSpaceAlong = 2;
%												- FlipkSpaceWhileAccessing = ':,:,:,:,:,:,2';
%                                               - (tbc)
% -         ReadInInfo                  ...     Factor with which the MRSI data should be zerofilled in k-space for interpolation (e.g. zerofill from 64x64 to 128x128)
% Output:
% -         iSpace                      ...     Output data in image domain. In case of Single Voxel Spectroscopy, this is the only output
% -         Noise                       ...     The Noise Correlation Matrix in order to check if it was properly computed from the csi data. Is 0 if no decorrelation performed.
% -         PreProcessingInfo           ...     The  updated PreProcessingInfo, as it is changed when reading in.
% -         kSpace                      ...     Output data in k-space. In case of SVS this is zero. size: channel x ROW x COL x SLC x vecSize x Averages
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: memused_linux,Analyze_csi_mdh, read_ascconv, hadamard_encoding.m

% Further remarks: This function uses FFTs to get from k- to image-space. This is mathematically wrong, but Siemens seems to do the same when
% creating DICOMS. The only difference is that the images are flipped left/right and up/down.



%% 0. Preparations


% Find out memory used by MATLAB
% memused_before = memused_linux(1); 


% % Assign standard values to variables if nothing is passed to function.

global NoiseCorrMat_post;


%% 1. Loop Over All Existing DataSets In PreProcessingInfo and Preprocess Data

DataSetNames = transpose(fields(kSpace));
for CurDataSet = DataSetNames

    
    CurDataSetString = CurDataSet{:};
	
    if(~isfield(kSpace,CurDataSet{:}) || strcmpi(CurDataSetString,'NOISEADJSCAN')  )
		continue;
	end
	fprintf('\n\nPreprocessing Data         \t...\t%s',CurDataSet{:})
	
	
	% 1.0 InLoopPreps
	% Make String out of Cell
	
	% Set fredirshift
% 	if(~isfield(PreProcessingInfo.(CurDataSetString),'fredir_shift') && size(kSpace.(CurDataSetString),6) == 1)	% Only for imaging data
% 		PreProcessingInfo.(CurDataSetString).fredir_shift = 2*ReadInInfo.General.Ascconv.Pos_Sag(1)/ (-ReadInInfo.General.Ascconv.FoV_Read(1)/ReadInInfo.(CurDataSetString).nReadEnc );
%     end
	
    

	% Do conj if necessary
	
	
	
	
	% 1.1 Compute Noise Correlation Matrix if necessary
	% Only read in data if NoiseCorrelationMatrix exists and is a 1
	if(isfield(PreProcessingInfo.(CurDataSetString),'NoiseCorrMat') && numel(PreProcessingInfo.(CurDataSetString).NoiseCorrMat) == 1 && PreProcessingInfo.(CurDataSetString).NoiseCorrMat > 0 )
		fprintf('\nCompute noise correlation matrix')


		%%%%%%% GET NOISE %%%%%%

		% Get data from Noise Prescan
		if(isfield(kSpace,'NOISEADJSCAN') && PreProcessingInfo.(CurDataSetString).NoiseCorrMat == 1)
			fprintf('\nGet noise from\t...\tNoise Prescan.')    
			Noise_mat = reshape(kSpace.NOISEADJSCAN,[size(kSpace.NOISEADJSCAN,1) numel(kSpace.NOISEADJSCAN)/size(kSpace.NOISEADJSCAN,1)]);
			Noise_mat = Noise_mat(:,20:end);
%             Noise_wsvd1=Noise_mat;
		% Or from end of FID of the outer kSpace
		elseif(isfield(kSpace,'ONLINE') && size(kSpace.ONLINE,6) > 512 && PreProcessingInfo.(CurDataSetString).NoiseCorrMat == 2)
			fprintf('\nGet noise from\t...\tONLINE Data Itself.')
			Noise_mat = GatherNoiseFromCSI(kSpace.ONLINE,ReadInInfo.General.Ascconv.Full_ElliptWeighted_Or_Weighted_Acq);
		end

		% Copy NoiseCorrMat
		if(exist('Noise_mat','var'))
			for CurDataSet2 = DataSetNames
				if(~strcmpi(CurDataSet2{:},'ONLINE'))
					continue;
				end

				% Rescale Noise to be the same like ONLINE noise
				if( numel(PreProcessingInfo.(CurDataSetString).NoiseCorrMat) == 1 && PreProcessingInfo.(CurDataSetString).NoiseCorrMat == 1 && ...
					isfield(ReadInInfo,'NOISEADJSCAN') && isfield(ReadInInfo.NOISEADJSCAN,'Dwelltime') && isfield(ReadInInfo,CurDataSet2{:}) && isfield(ReadInInfo.(CurDataSet2{:}),'Dwelltime') )
					NoiseScalingFactor = sqrt(ReadInInfo.NOISEADJSCAN.Dwelltime / (ReadInInfo.General.Ascconv.WipMemBlockInterpretation.Rollercoaster.ADCDwellTime/2)); % oversampling
				else
					NoiseScalingFactor = 2; %was 1 pmos
				end
				Noise_mat = Noise_mat * NoiseScalingFactor;
%                 Noise_wsvd2=Noise_mat;
				% Compute noise correlation matrix
				NoiseCorrMat = 1/(size(Noise_mat,2)) * (Noise_mat * Noise_mat');
				PreProcessingInfo.(CurDataSet2{:}).NoiseCorrMat = NoiseCorrMat;
				Noise.(CurDataSet2{:}) = Noise_mat;
                NoiseCorrMat_post=NoiseCorrMat; %for pst gridding noise-deco
			end

			% Delete Uneccessary stuff
			clear NoisePar Noise_mat  NoiseScalingFactor NoiseCorrMat
            Noise.PATREFSCAN = Noise.ONLINE;
        end
        
	end


	
% 	% 1.2. Perform Noise Decorrelation
% 	if(isfield(PreProcessingInfo.(CurDataSetString),'NoiseCorrMat') )
% 		fprintf('\nPerform Noise Decorrelation.')
% %         [MemUser,MemFree] = memused_linux(1);
% %         MemUser = whos('kSpace'); MemUser = MemUser.bytes / 2^20;
% %         if(MemFree > 1.1*3*MemUser)
% %             for echo = 1:numel(kSpace.(CurDataSetString))
% %                 kSpace.(CurDataSetString){echo} = PerformNoiseDecorrelation(kSpace.(CurDataSetString){echo},PreProcessingInfo.(CurDataSetString).NoiseCorrMat);
% %             end
% %         else
% %             for echo = 1:numel(kSpace.(CurDataSetString))
% %                 for kx = 1:size(kSpace.(CurDataSetString){echo},2)
% %                     kSpace.(CurDataSetString){echo}(:,kx,:,:,:,:,:) = PerformNoiseDecorrelation(kSpace.(CurDataSetString){echo}(:,kx,:,:,:,:,:),PreProcessingInfo.(CurDataSetString).NoiseCorrMat);                
% %                 end
% %             end
% %         end
%             
%     end

	
    
     
    % Sum averages
    if(isfield(PreProcessingInfo.(CurDataSetString),'SumAverages_flag') && PreProcessingInfo.(CurDataSetString).SumAverages_flag && size(kSpace.(CurDataSetString),7) > 1) 
        AvgTemp = zeros(size(kSpace.(CurDataSetString)(:,:,:,:,:,:,1)));
        for avg = 1:size(kSpace.(CurDataSetString),7)
            AvgTemp = AvgTemp + kSpace.(CurDataSetString)(:,:,:,:,:,:,avg);
        end
        kSpace.(CurDataSetString) = AvgTemp / size(kSpace.(CurDataSetString),7); clear AvgTemp;
    end

    
	
	% Remove spatial oversampling
	if(isfield(PreProcessingInfo.(CurDataSetString), 'RmOs') && PreProcessingInfo.(CurDataSetString).RmOs == 1)
		fprintf('\nRemove Spatial Oversampling.')
		for echo = 1:numel(kSpace.(CurDataSetString))
			image_center = floor(size(iSpace.(CurDataSetString){echo},2) / 2) + 1;
			left_border = image_center - ceil(size(iSpace.(CurDataSetString){echo},2) / (2*2));   
			right_border = image_center + ceil(size(iSpace.(CurDataSetString){echo},2) / (2*2)) - 1;   
			if(size(iSpace.(CurDataSetString){echo},2) > 1)
				iSpace.(CurDataSetString){echo} = iSpace.(CurDataSetString){echo}(:,left_border:right_border,:,:,:,:,:);
			else
				iSpace.(CurDataSetString){echo} = fftshift(fft(ifftshift(kSpace.(CurDataSetString){echo},2),[],2),2);
				iSpace.(CurDataSetString){echo} = iSpace.(CurDataSetString){echo}(:,left_border:right_border,:,:,:,:,:);
			end
			if(nargout > 3)
				kSpace.(CurDataSetString){echo} = FFTOfMRIData(iSpace.(CurDataSetString){echo},ConjFlag, [2 3 4],1);
				if(isfield(kSpace,[CurDataSetString '_Unfiltered']))
					kSpace.([CurDataSetString '_Unfiltered']){echo} = kSpace.([CurDataSetString '_Unfiltered']){echo}(:,1:2:end,:,:,:,:,:);
				end
			end
		end
	end
	

	

end



%% 7. Postparations

if(nargout > 1 && ~exist('Noise','var'))
	Noise = 0;
end
if(nargout > 2 && ~exist('kSpace','var'))
	kSpace = 0;
end
if(nargout > 3 && ~exist('PreProcessingInfo','var'))
	PreProcessingInfo = 0;
end

% memused_after = memused_linux(1); 
% display([char(10) 'The function used ' num2str(memused_after-memused_before) '% of the total memory.'])


