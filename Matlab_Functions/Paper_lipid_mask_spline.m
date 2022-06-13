load /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/CONCEPT_INVIVO_MAPS/CONCEPT_PAPER/Vol05/Concept_HAM_SNRBugfix_GRDonHAM2_withPostGRD/CombinedCSI.mat
load /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/CONCEPT_INVIVO_MAPS/CONCEPT_PAPER/Vol05/Concept_HAM_SNRBugfix_GRDonHAM2_withPostGRD/maps/AllMaps.mat
csi_ham=csi;
meta_ham=AllMaps.Metabos.Normal(:,:,:,11);
load /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/CONCEPT_INVIVO_MAPS/CONCEPT_PAPER/Vol05/Concept_NAT_SNRBugfix_GRDonHAM2_withPostGRD/CombinedCSI.mat
load /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/CONCEPT_INVIVO_MAPS/CONCEPT_PAPER/Vol05/Concept_NAT_SNRBugfix_GRDonHAM2_withPostGRD/maps/AllMaps.mat
csi_nat=csi;
meta_nat=AllMaps.Metabos.Normal(:,:,:,11);
load /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/CONCEPT_INVIVO_MAPS/CONCEPT_PAPER/Vol05/ePE/CombinedCSI.mat
load /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/CONCEPT_INVIVO_MAPS/CONCEPT_PAPER/Vol05/ePE/maps/AllMaps.mat
csi_pe=csi;
meta_pe=AllMaps.Metabos.Normal(:,:,:,11);
load /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/CONCEPT_INVIVO_MAPS/CONCEPT_PAPER/Vol05/ePER5/CombinedCSI.mat
load /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/CONCEPT_INVIVO_MAPS/CONCEPT_PAPER/Vol05/ePER5/maps/AllMaps.mat
csi_r5=csi;
meta_r5=AllMaps.Metabos.Normal(:,:,:,11);
clear csi;

%load /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/CONCEPT_INVIVO_MAPS/CONCEPT_PAPER/phantom_lipid/Concept_HAM_lipidphanti_WS/CombinedCSI.mat
%load /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/CONCEPT_INVIVO_MAPS/water_phanti_R5/CombinedCSI.mat
%load /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/CONCEPT_INVIVO_MAPS/agar_ws_pe/CombinedCSI.mat
%csi_ham=csi;


size(csi_nat);
size(csi_pe);
size(csi_r5);

csi_ham=abs(fftshift(fft(csi_ham(:,:,1,:),[],4),4));
csi_nat=abs(fftshift(fft(csi_nat(:,:,1,:),[],4),4));
csi_pe=abs(fftshift(fft(csi_pe(:,:,1,:),[],4),4));
csi_r5=abs(fftshift(fft(csi_r5(:,:,1,:),[],4),4));


for i=1:64
   for j=1:64
      
      y=squeeze(csi_ham(i,j,1,:));
      y(370:500)=[];
      y=y';
         
      p=polyfit(1:numel(y),y,6);
      y1 = polyval(p,1:numel(y));
      
      csi_ham_new(i,j,1,1:numel(y))=y-y1;  
      
   end 
end


for i=1:64
   for j=1:64
      
      y=squeeze(csi_nat(i,j,1,:));
      y(370:500)=[];
      y=y';
         
      p=polyfit(1:numel(y),y,6);
      y1 = polyval(p,1:numel(y));
      
      csi_nat_new(i,j,1,1:numel(y))=y-y1;  
      
   end 
end


for i=1:64
   for j=1:64
      
      y=squeeze(csi_pe(i,j,1,:));
      y(800:1070)=[];
      y=y';
         
      p=polyfit(1:numel(y),y,6);
      y1 = polyval(p,1:numel(y));
      
      csi_pe_new(i,j,1,1:numel(y))=y-y1;  
      
   end 
end

for i=1:64
   for j=1:64
      
      y=squeeze(csi_r5(i,j,1,:));
      y(800:1070)=[];
      y=y';
         
      p=polyfit(1:numel(y),y,6);
      y1 = polyval(p,1:numel(y));
      
      csi_r5_new(i,j,1,1:numel(y))=y-y1;  
      
   end 
end

% 
% csi_ham_norm=max(csi_ham,[],4);
% csi_nat_norm=max(csi_nat,[],4);
% csi_pe_norm=max(csi_pe,[],4);
% csi_r5_norm=max(csi_r5,[],4);
% 
% csi_ham_norm=repmat(csi_ham_norm,[1 1 size(csi_ham,4)]);
% csi_nat_norm=repmat(csi_nat_norm,[1 1 size(csi_nat,4)]);
% csi_pe_norm=repmat(csi_pe_norm,[1 1 size(csi_pe,4)]);
% csi_r5_norm=repmat(csi_r5_norm,[1 1 size(csi_r5,4)]);
% 
% csi_ham=squeeze(csi_ham)./csi_ham_norm;
% csi_nat=squeeze(csi_nat)./csi_nat_norm;
% csi_pe=squeeze(csi_pe)./csi_pe_norm;
% csi_r5=squeeze(csi_r5)./csi_r5_norm;
% 
csi_ham_sum=sum(abs(csi_ham_new(:,:,:,600:700)),4);
csi_nat_sum=sum(abs(csi_nat_new(:,:,:,600:700)),4);

csi_pe_int=sum(abs(csi_pe_new(:,:,:,970:1100)),4);
csi_r5_int=sum(abs(csi_r5_new(:,:,:,970:1100)),4);

figure
subplot(1,4,1)
imagesc(squeeze(mask.*csi_ham_sum(:,:,1,:)./meta_ham),[0 200000])
subplot(1,4,2)
imagesc(squeeze(mask.*csi_nat_sum(:,:,1,:)./meta_nat),[0 200000])
subplot(1,4,3)
imagesc(squeeze(mask.*csi_pe_int(:,:,1,:)./meta_pe),[0 200000])
subplot(1,4,4)
imagesc(squeeze(mask.*csi_r5_int(:,:,1,:)./meta_r5),[0 200000])

% 
% subplot(3,4,9)
% imagesc(squeeze(mask.*csi_ham_int_norm(:,:,1,:)./meta_ham),[0 8])
% subplot(3,4,10)
% imagesc(squeeze(mask.*csi_nat_int_norm(:,:,1,:)./meta_nat),[0 8])
% subplot(3,4,11)
% imagesc(squeeze(mask.*csi_pe_int_norm(:,:,1,:)./meta_pe),[0 8])
% subplot(3,4,12)
% imagesc(squeeze(mask.*csi_r5_int_norm(:,:,1,:)./meta_r5),[0 8])

% figure; imagesc(squeeze(mask./csi_ham(:,:,1,:).*meta_ham),[0 1/10000])
% figure; imagesc(squeeze(mask./csi_nat(:,:,1,:).*meta_nat),[0 1/10000])
% figure; imagesc(squeeze(mask./csi_pe(:,:,1,:).*meta_pe),[0 1/10000])
% figure; imagesc(squeeze(mask./csi_r5(:,:,1,:).*meta_r5),[0 1/10000])

