load /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/CONCEPT_INVIVO_MAPS/CONCEPT_PAPER/Vol01/Concept_HAM_SNRBugfix_GRDonHAM2_withPostGRD/CombinedCSI.mat
load /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/CONCEPT_INVIVO_MAPS/CONCEPT_PAPER/Vol01/Concept_HAM_SNRBugfix_GRDonHAM2_withPostGRD/maps/AllMaps.mat
csi_ham=csi;
meta_ham=AllMaps.Metabos.Normal(:,:,:,11);
load /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/CONCEPT_INVIVO_MAPS/CONCEPT_PAPER/Vol01/Concept_NAT_SNRBugfix_GRDonHAM2_withPostGRD/CombinedCSI.mat
load /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/CONCEPT_INVIVO_MAPS/CONCEPT_PAPER/Vol01/Concept_NAT_SNRBugfix_GRDonHAM2_withPostGRD/maps/AllMaps.mat
csi_nat=csi;
meta_nat=AllMaps.Metabos.Normal(:,:,:,11);
load /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/CONCEPT_INVIVO_MAPS/CONCEPT_PAPER/Vol01/ePE/CombinedCSI.mat
load /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/CONCEPT_INVIVO_MAPS/CONCEPT_PAPER/Vol01/ePE/maps/AllMaps.mat
csi_pe=csi;
meta_pe=AllMaps.Metabos.Normal(:,:,:,11);
load /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/CONCEPT_INVIVO_MAPS/CONCEPT_PAPER/Vol01/ePER5/CombinedCSI.mat
load /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/CONCEPT_INVIVO_MAPS/CONCEPT_PAPER/Vol01/ePER5/maps/AllMaps.mat
csi_r5=csi;
meta_r5=AllMaps.Metabos.Normal(:,:,:,11);
clear csi;

%load /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/CONCEPT_INVIVO_MAPS/CONCEPT_PAPER/phantom_lipid/Concept_HAM_lipidphanti_WS/CombinedCSI.mat
%load /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/CONCEPT_INVIVO_MAPS/water_phanti_R5/CombinedCSI.mat
%load /net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/CONCEPT_INVIVO_MAPS/agar_ws_pe/CombinedCSI.mat
%csi_ham=csi;

size(csi_ham)
size(csi_nat)
size(csi_pe)
size(csi_r5)

csi_ham=abs(fftshift(fft(csi_ham(:,:,1,:),[],4),4));
csi_nat=abs(fftshift(fft(csi_nat(:,:,1,:),[],4),4));
csi_pe=abs(fftshift(fft(csi_pe(:,:,1,:),[],4),4));
csi_r5=abs(fftshift(fft(csi_r5(:,:,1,:),[],4),4));

csi_ham=csi_ham(:,:,1,720:850);
csi_nat=csi_nat(:,:,1,720:850);
csi_pe=csi_pe(:,:,1,1240:1350);
csi_r5=csi_r5(:,:,1,1240:1350);
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
% csi_ham_int=sum(csi_ham(:,:,:),3);
% csi_nat_int=sum(csi_nat(:,:,:),3);
% csi_pe_int=sum(csi_pe(:,:,:),3);
% csi_r5_int=sum(csi_r5(:,:,:),3);


subplot(3,4,1)
imagesc(squeeze(mask.*csi_ham_norm(:,:,1,:)./meta_ham),[0 8000])
subplot(3,4,2)
imagesc(squeeze(mask.*csi_nat_norm(:,:,1,:)./meta_nat),[0 8000])
subplot(3,4,3)
imagesc(squeeze(mask.*csi_pe_norm(:,:,1,:)./meta_pe),[0 8000])
subplot(3,4,4)
imagesc(squeeze(mask.*csi_r5_norm(:,:,1,:)./meta_r5),[0 8000])

subplot(3,4,5)
imagesc(squeeze(mask.*csi_ham_int(:,:,1,:)./meta_ham),[0 10])
subplot(3,4,6)
imagesc(squeeze(mask.*csi_nat_int(:,:,1,:)./meta_nat),[0 10])
subplot(3,4,7)
imagesc(squeeze(mask.*csi_pe_int(:,:,1,:)./meta_pe),[0 10])
subplot(3,4,8)
imagesc(squeeze(mask.*csi_r5_int(:,:,1,:)./meta_r5),[0 10])
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

