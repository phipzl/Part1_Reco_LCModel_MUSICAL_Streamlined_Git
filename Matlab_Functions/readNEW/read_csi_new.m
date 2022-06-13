function [Noise,PreProcessingInfo,ReadInInfo,kSpace] = read_csi_new(file,tmp_dir,PreProcessingInfo,ReadInDataSets)
%
% read_csi Read in csi-data
%
% This function was written by Bernhard Strasser, July 2012.
%
%
% The function sets default values for PreProcessingInfo and decides if input file is .dat file or DICOM according to its ending. 
% It then calls the read_csi_dat or read_csi_dicom functions. Refer to these for more info.
%
%
% [iSpace,Noise,PreProcessingInfo,ReadInInfo,kSpace] = read_csi(file,PreProcessingInfo,ReadInDataSets)
%
% Input: 
% -         file                        ...     Path of MR(S)(I) file.
% -         PreProcessingInfo           ...     Info about how the read in data should be pre-processed. See help PreProcessMRIData_Wrapper.
% -         ReadInDataSets              ...     Only read in those datasets, e.g. 'ONLINE', 'NOISEADJSCAN', 'PATREFANDIMASCAN'. These are compared to
%                                               the EvalInfoMask of the raw-data, and only if it matches, that data set is read in.
%
% Output:
% -         iSpace                      ...     Structure of data in image domain. Its fieldnames are according to the EvalInfoMask information (e.g. 'ONLINE', 'PATREFANDIMASCAN', etc.)
%                                               size of each field: channel x nFreqEnc x nPhasEnc x nPartEnc x nSlc x SpectroVecSize x Averages.
% -         Noise                       ...     Structure of Noise gathered from the Prescan or CSI data. Noise = 0, if not gathered. It is scaled individually for each measurement set.
% -         PreProcessingInfo           ...     The updated PreProcessingInfo, as it is changed when reading in.
% -         ReadInInfo                  ...     The info of the mdh.
% -         kSpace                      ...     Structure of data in k-space, see iSpace.
%                                               size of each field: channel x nFreqEnc x nPhasEnc x nPartEnc x nSlc x SpectroVecSize x Averages.
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: read_csi_dat, read_csi_dicom, PreProcessMRIData_Wrapper, and many more (???)
%
% See also read_csi_dat, read_csi_dicom, PreProcessMRIData_Wrapper.

global NoiseCorrMat_post;
global totCoili;

%% 0. PREPARATIONS

% % Assign standard values to variables if nothing is passed to function.

% If nothing is passed to function
if(nargin < 1)
    display('Please feed me with a file to read in.')
    iSpace = 0;
    kSpace = 0;
    return;
end
if(~exist('ReadInDataSets','var'))
	ReadInDataSets = 'All';
end

if(exist('PreProcessingInfo','var') && numel(PreProcessingInfo) == 1 && ~isstruct(PreProcessingInfo) && PreProcessingInfo == 0)
	clear PreProcessingInfo;
end
	
% Test if any kSpace Preprocessing should be done with ONLINE
if(exist('PreProcessingInfo','var') && isfield(PreProcessingInfo,'ONLINE') &&  isfield(PreProcessingInfo.ONLINE,'Hamming_flag') )
	ONLINEkSpaceNecessary = PreProcessingInfo.ONLINE.Hamming_flag;
else
	ONLINEkSpaceNecessary = false;
end



%% 1. Read In Data


if(numel(strfind(file, '.dat')) > 0)
    
    % Read Raw Data
    fprintf('read_csi_dat_new_v2...\n');
    [kSpace,ReadInInfo] = read_csi_dat_new_v2(file, tmp_dir);
   
end





%% 2. Define Standard PreProcessingInfo

PreProcessingInfo_Standard.ONLINE.NoiseCorrMat = 1;
PreProcessingInfo_Standard.ONLINE.Hadamard_flag = true;				% If this flag is set to true, hadamard decoding is performed, if several slices were measure. If set to false, 
PreProcessingInfo_Standard.ONLINE.NoFFT_flag = false;
%PreProcessingInfo_Standard.ONLINE.fredir_shift = 0;
PreProcessingInfo_Standard.ONLINE.SaveUnfilteredkSpace = false;
PreProcessingInfo_Standard.ONLINE.Hamming_flag = false;
PreProcessingInfo_Standard.ONLINE.SumAverages_flag = true;

% Remove Oversampling
if(isfield(kSpace,'ONLINE') && size(kSpace.ONLINE,6) > 1)
	PreProcessingInfo_Standard.ONLINE.RmOs = false;		% CSI
else
	PreProcessingInfo_Standard.ONLINE.RmOs = true;		% Imaging
end
PreProcessingInfo_Standard.PATREFANDIMASCAN.RmOs = true;
OversamplingFactor_ONLINE = 2;

if(isfield(kSpace,'ONLINE') && isfield(ReadInInfo,'ONLINE') && isfield(ReadInInfo.ONLINE, 'nReadEnc') && size(kSpace.ONLINE,2) == size(kSpace.ONLINE,3))  % This is really really bad...
	PreProcessingInfo_Standard.ONLINE.RmOs = false;
	OversamplingFactor_ONLINE = 1;
end	


if(isfield(kSpace,'PATREFSCAN'))            % Concept Prescan. Cheat it to be the same as ONLINE
    PreProcessingInfo_Standard.PATREFSCAN = PreProcessingInfo_Standard.ONLINE;          
end






%% 5. PreProcess Data

% PreProcess Data

[Noise,PreProcessingInfo,kSpace] = PreProcessMRIData_Wrapper_new(kSpace,PreProcessingInfo_Standard,ReadInInfo);


%% 6. Postparations

% if(isfield(kSpace,'Noise'))
%     if(PreProcessingInfo.Values.NoiseCorrMat == 1)
%         Noise.Mat = kSpace.Noise;
%     end
%     kSpace = rmfield(kSpace,'Noise');
% end
% if(isfield(kSpace,'NoiseCorrMat'))
%     if(PreProcessingInfo.Values.NoiseCorrMat == 1)
%         Noise.CorrMat = kSpace.NoiseCorrMat;
%     end
%     kSpace = rmfield(kSpace,'NoiseCorrMat');
% end
% 
% if(isfield(PreProcessingInfo.Values,'NoiseCorrMat'))
%     if(numel(PreProcessingInfo.Values.NoiseCorrMat) > 1)
%         Noise.CorrMat = PreProcessingInfo.Values.NoiseCorrMat;
%     end
% end





