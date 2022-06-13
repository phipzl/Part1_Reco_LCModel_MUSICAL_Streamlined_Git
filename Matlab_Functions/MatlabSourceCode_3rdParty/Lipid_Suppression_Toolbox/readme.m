%%
% script_lipid_suppression.m : script file to run example lipid suppression
% reconstructions using 0.16cc, 20 average in vivo CS data. 
% Three methods are compared:
%
% i)   apply iterative lipid suppression to low resolution, 20 average data

% ii)  generate dual-density image by combining low res, 20 average data
%      with high-res 2 average lipid image, then apply iterative lipid suppression

% iii) apply FOCUSS compressed sensing reconstruction to 2 average,
%      undersampled lipid image, generate dual-density image by combining low res,
%      20 average data with FOCUSS reconstructed lipid image, then apply iterative lipid suppression

