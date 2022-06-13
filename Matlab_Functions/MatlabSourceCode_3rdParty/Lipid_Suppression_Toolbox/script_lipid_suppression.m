%% load data

addpath InVivo_Data\
addpath utils\
addpath FOCUSS_utils\


% 2 average CSI data at 0.16cc resolution, TE=50ms, no lipid suppression 
load Image_2aver_noOVS_p16cc_33min_TE50ms;

% 20 average CSI data at 0.16cc resolution, TE=50ms, no lipid suppression 
load Image_20aver_noOVS_p16cc_33min_TE50ms;


load lipid_mask_2012_02_13              % lipid layer mask
load meta_mask_2012_02_13               % metabolite mask


N = size(img);

% plot sum over all frequencies for 20 and 2 aver data
plot_projection([img, img_2avg])



%% Comparing three reconstructions:

% i)   apply iterative lipid suppression to low resolution, 20 average data

% ii)  generate dual-density image by combining low res, 20 average data
%      with high-res 2 average lipid image, then apply iterative lipid suppression

% iii) apply FOCUSS compressed sensing reconstruction to 2 average,
%      undersampled lipid image, generate dual-density image by combining low res,
%      20 average data with FOCUSS reconstructed lipid image, then apply iterative lipid suppression



%% i)   apply iterative lipid suppression to low resolution, 20 average data

csi_rad = 16;                     % size of low res image, has lipid ringing

mask_csi = repmat(circular_mask([64,64],csi_rad),[1,1,N(3)]);
FT_csi = FT_v2(mask_csi);         % samples the center radius = [csi_rad] region in kx-ky

img_zf = FT_csi'*(FT_csi*img);    % low res image with lipid ringing

param = init;
param.xfmWeight = 1e-3;           % Lipid basis penalty 

param.FT = FT_csi;
param.data =  param.FT*img; 
param.Itnlim = 10;

param.Lipid = get_LipidBasis(img_zf,lipid_mask);    % Lipid basis functions obtained from low-res image
param.Bmask = meta_mask;                            % Brain mask 
x = img_zf;

plot_projection([img_zf,img],1)

tic
for t = 1:10
    x = lipid_suppression(x,param);    
    figure(1), imagesc( [ meta_mask.*sum(abs(x),3) , meta_mask.*sum(abs(img_zf),3) ] ), axis image, colorbar, drawnow
end
toc


% lipid suppressed image within the brain mask
x_meta = x .* repmat(meta_mask,[1,1,N(3)]);


% overplot the initial spectra and the lipid-suppressed spectra
overplot_spectra( img_zf .* repmat(meta_mask,[1,1,N(3)]), x_meta, 30, 40, 8, 8, 35, 160)



%% ii)  generate dual-density image by combining low res, 20 average data
%%      with high-res 2 average lipid image, then apply iterative lipid suppression

mask_high = ones(N) - mask_csi; 
FT_high = FT_v2(mask_high);
FT_csi  = FT_v2(mask_csi); 
FT_full = FT_v2(ones(N));

data1 = FT_high * (img_2avg .* repmat(lipid_mask,[1,1,N(3)]));
data2 = FT_csi * img;  

combined_2avg = FT_full'*(data1+data2);         % dual-density image (low res CSI + high res lipid image)
plot_projection(combined_2avg)


param = init;
param.xfmWeight = 1e-3;       % Lipid basis penalty 

param.FT = FT_full;
param.data = param.FT*combined_2avg; 
param.Itnlim = 10;

param.Lipid = get_LipidBasis(combined_2avg,lipid_mask);    % Lipid basis functions obtained from dual-density image
param.Bmask = meta_mask;                                   % Brain mask 
res_2avg = combined_2avg;

tic
for t = 1:10
    res_2avg = lipid_suppression(res_2avg,param);
    figure(1), imagesc(  [ sum(abs(res_2avg),3).*meta_mask , sum(abs(combined_2avg),3).*meta_mask ] ), axis image, colorbar, drawnow
end
toc


% lipid suppressed image within the brain mask
res_2avg_meta = res_2avg .* repmat(meta_mask,[1,1,N(3)]);


% overplot low-res lipid suppression and dual-density + lipid suppression
overplot_spectra( x_meta, res_2avg_meta, 30, 40, 8, 8, 35, 160)



%% iii) apply FOCUSS compressed sensing reconstruction to 2 average,
%%      undersampled lipid image, generate dual-density image by combining low res,
%%     20 average data with FOCUSS reconstructed lipid image, then apply iterative lipid suppression
%
% Compressed sensing reconstruction: undersample 2-average high res data 10-fold in outer k-space

c = circular_mask([64,64],csi_rad);

p = 3;
R_high = 10;    % acceleration in the outer k-space 
R_total = 1 / (((64^2 - sum(c(:))) / R_high + sum(c(:))) / 64^2);
pdf = genPDF([64,64], p, 1/R_total, 2, csi_rad/45, 0);

figure(1), imagesc(pdf, [0,1]), axis image off, 

m3f = [];
for t = 1:N(3)
    m3d(:,:,t) = genSampling(pdf,20,10);
end
figure(2), imagesc(sum(m3d,3)), axis image off, title('Sampling mask')


DFT = FT_v2(m3d);
img_zf = DFT'*(DFT*img);
figure(1), imagesc( sum(abs([img_zf,img]),3) ), axis image, title(['Zero filling: ',num2str(100*norm(img_zf(:)-img(:)) / norm(img(:))),' % RMSE']), colorbar, drawnow

m2d = mean(m3d,3);
disp(['Acceleration factor of outer k-space: ', num2str(1 / mean( m2d(m2d<1) ))])    



%% Using FOCUSS algorithm to recon high res lipid image from undersampled k-space

avg2_lipid = img_2avg .* repmat(lipid_mask,[1,1,N(3)]);          % from 2 average data

FT = FT_v2(m3d);

data = FT * img_2avg;           % undersampled k-space
res0 = FT'*data;

lambda = 0;
cg_iter = 50;
focuss_iter = 10; 
p = 1;


tic
    % using reg. Least Squares formulation
    Fres = focuss_3D( res0, data, FT, lambda, cg_iter, focuss_iter, p, avg2_lipid, lipid_mask );
toc

Fres_lipid = Fres.*repmat(lipid_mask,[1,1,N(3)]);

plot_projection( [Fres_lipid, avg2_lipid], 2 ), title(['Lipid ring RMSE: ', num2str(100*norm(Fres_lipid(:)-avg2_lipid(:))/norm(avg2_lipid(:))), '%']), drawnow



%% use CS reconstructed lipid image to find combined image, then apply lipid suppression

mask_high = ones(N) - mask_csi; 
FT_high = FT_v2(mask_high);
FT_csi  = FT_v2(mask_csi); 

param = init;
param.xfmWeight = 1e-3;       % Lipid basis penalty 

param.FT = FT_v2(ones(N));
param.data = (FT_high * Fres_lipid) + (FT_csi * img); 
param.Itnlim = 10;

param.Lipid = get_LipidBasis(Fres, lipid_mask);    % Lipid basis functions from CS-reconstructed image
param.Bmask = meta_mask;                           % Brain mask 
res_combined = param.FT'*param.data;

tic
for t = 1:10
    res_combined = lipid_suppression(res_combined,param);
    figure(1), imagesc( [sum(abs(res_combined),3).*meta_mask , sum(abs(res_2avg_meta),3)] ), axis image, colorbar, drawnow
end
toc


% lipid suppressed image within the brain mask
res_combined_meta = res_combined .* repmat(meta_mask,[1,1,N(3)]);


% overplot dual-density + lipid suppression and CS-recon + lipid suppression
overplot_spectra( res_2avg_meta, res_combined_meta, 30, 40, 8, 8, 35, 160)




%% plot lipid images (projection over lipid frequencies)


% projection over lipid frequencies
figure(1), imagesc( (sum(abs([ imrotate(res_2avg(:,:,105:end),-90), imrotate(res_combined(:,:,105:end),-90); ...
    imrotate(x(:,:,105:end),-90), imrotate(img(:,:,105:end),-90) ]),3)), [0,2] ), axis image off, 
    title('2avg HighRes,    2avg R=10 HighRes,    LowRes,    No Suppression')

% projection over lipid frequencies in log scale
figure(2), imagesc( 20*log10(sum(abs([ imrotate(res_2avg(:,:,105:end),-90), imrotate(res_combined(:,:,105:end),-90); ...
    imrotate(x(:,:,105:end),-90), imrotate(img(:,:,105:end),-90) ]),3)), [-35,15] ), axis image off, colorbar,
    title('2avg HighRes,    2avg R=10 HighRes,    LowRes,    No Suppression  [Log Scale]')



%% plot NAA images (projection over NAA frequencies)

naa_range = 110:125;

r2m_naa = sum( abs(imrotate(res_2avg_meta(:,:,naa_range),-90)), 3);
rcm_naa = sum( abs(imrotate(res_combined_meta(:,:,naa_range),-90)), 3);
xm_naa = sum( abs(imrotate(x_meta(:,:,naa_range),-90)), 3);

naa_scale_factor = max(r2m_naa(:));

r2m_naa = r2m_naa / naa_scale_factor;
rcm_naa = rcm_naa / naa_scale_factor;
xm_naa = xm_naa / naa_scale_factor;

figure(3), imagesc( [r2m_naa, rcm_naa, xm_naa], [0,1.5] ), axis image off, colorbar, colormap gray, title('2avg HighRes,   2avg R=10 HighRes,   LowRes')



%% plot example spectra - I

overplot_spectra( imrotate(x_meta,-90), imrotate(res_2avg_meta,-90),  44, 36, 6, 6, 1, N(3)-25)

overplot_spectra( imrotate(res_combined_meta,-90), imrotate(res_2avg_meta,-90),  44, 36, 6, 6, 1, N(3)-25, 3,4)


%% plot example spectra - II

overplot_spectra( imrotate(x_meta,-90), imrotate(res_2avg_meta,-90),  38, 18, 6, 6, 1, N(3)-25)

overplot_spectra( imrotate(res_combined_meta,-90), imrotate(res_2avg_meta,-90),  38, 18, 6, 6, 1, N(3)-25, 3,4)


