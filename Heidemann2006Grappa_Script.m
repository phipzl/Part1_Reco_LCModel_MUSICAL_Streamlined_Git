%% Definitions
% load /ceph/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/Sourcecode/MRSI_Processing_GitRepos/Part1_Reco_LCModel_MUSICAL_GRD12/csi_k_PI.mat
% load /ceph/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/Sourcecode/MRSI_Processing_GitRepos/Part1_Reco_LCModel_MUSICAL_GRD12/image_k_PI.mat
% 

% clearvars -except csi_k image_k
% close all; clearvars

NoOfSegments = 55;
MinKernelSrcPts = 20;
AdditionalPhiPts = 4;       % These are the points along phi that we additionally 
                            % use on each side of the segment to calculate the weights
UndersamplingFactor = 2;

% %% Load Data
% 
% if(~exist('csi_k','var'))
%    load('./dataBernie/csi_k.mat')
% end
% if(~exist('image_k','var'))
%     load('./dataBernie/image_k.mat')
% end


%% Reconstruct data segment-wise using only first time point of image_k

csi_k_Reco = Heidemann2006GrappaFunc(csi_k_p,image_k_p,UndersamplingFactor,MinKernelSrcPts,NoOfSegments,AdditionalPhiPts);



% clear csi_k_Reco
% 
% %% Evaluate
% figure; subplot(2,1,1); imagesc(abs(squeeze(csi_k(2,:,:,1,1,1,1)))), title('Ground Truth')
% subplot(2,1,2); imagesc(abs(squeeze(csi_k_Reco(2,:,:,1,1,1,1)))), title('Reco')
% for i=1:16
% RMSE_mean = sqrt(mean(abs(squeeze(csi_k_Reco2(i,:))' - csi_k(:)).^2)) / sqrt(mean(abs(csi_k(:)).^2))
% %RMSE_median = sqrt(median(abs(squeeze(csi_k_Reco2(i,:)') - csi_k(:)).^2)) / sqrt(median(abs(csi_k(:)).^2))
% end


%% Reconstruct data segment-wise using all time points of image_k
% To-Be-Done



% 
% %% Regrid Data
% 
% load('./dataBernie/variablesForBernie.mat')
% 
% csi_k_Orig = csi_k;
% % clear csi_k_Reco
% 
% 
% RES = 64;
% nc = 32;
% NoOfTempInterleaves_vec = repmat(3,[1 nc]);
% totCoili = 32;
% nSlices = 1;
% slicenmb = 1;
% maxres = 2*nc;
% precision = 'double';
% post_grd_dc_flag = true;
% csi_k = ones(size(csi_k_Reco));
% run RollercoasterRecon_fast.m
% 
% 
% post_grd_dc_flag = false;
% csi_k = csi_k_Orig;
% run RollercoasterRecon_fast.m
% csi_Orig = csi; clear csi
% 
% post_grd_dc_flag = false;
% csi_k = csi_k_Reco;
% run RollercoasterRecon_fast.m
% csi_Reco = csi; clear csi
% 
% csi_k = csi_k_Orig;
% 
% 
% figure; imagesc(abs(squeeze(csi_Orig(1,:,:,1,10))))
% figure; imagesc(abs(squeeze(csi_Reco(1,:,:,1,10))))
