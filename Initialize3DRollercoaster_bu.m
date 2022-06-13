%Initialize3DRollercoaster.m

nSlices=ReadInInfo.ONLINE.nPartEnc;

csi_k=repmat(csi_k,[1 1 1 nSlices*size(csi_k,4) 1]);
image_k=repmat(image_k,[1 1 1 nSlices*size(image_k,4) 1]);
noise_sim=repmat(noise_sim,[1 1 1 nSlices*size(noise_sim,4) 1]);

%%%the fct below comes from the c++ sourcecode%%%

for i=-ReadInInfo.ONLINE.nPhasEnc/2:ReadInInfo.ONLINE.nPhasEnc/2

    NumberOfLoopsPerSlice(1+i+ReadInInfo.ONLINE.nPhasEnc/2)=ceil(ReadInInfo.ONLINE.nPhasEnc/2*sqrt(1-i^2/(ReadInInfo.ONLINE.nPhasEnc/2)^2));
    if NumberOfLoopsPerSlice(1+i+ReadInInfo.ONLINE.nPhasEnc/2)==0
        NumberOfLoopsPerSlice(1+i+ReadInInfo.ONLINE.nPhasEnc/2)=1;
    end
    %NumberOfLoopsPerSlice(1+i+ReadInInfo.ONLINE.nPhasEnc/2)=32; %%%FOR ZYL CONCEPT
end

NumberOfLoopsPerSlice(1:(numel(NumberOfLoopsPerSlice)-ReadInInfo.ONLINE.nPartEnc)/2)=[];
NumberOfLoopsPerSlice(ReadInInfo.ONLINE.nPartEnc+1: end)=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cumisumi_zwero=[0 cumsum(NumberOfLoopsPerSlice)];
cumisumi=cumsum(NumberOfLoopsPerSlice);

csi_k_fill=zeros([size(csi_k,1) max(NumberOfLoopsPerSlice) size(csi_k,3) size(csi_k,4) size(csi_k,5)]);
image_k_fill=zeros([size(csi_k,1) max(NumberOfLoopsPerSlice) size(csi_k,3) size(csi_k,4) size(image_k,5)]);
noise_sim_fill=zeros([size(csi_k,1) max(NumberOfLoopsPerSlice) size(csi_k,3) size(csi_k,4) size(noise_sim,5)]);


for i=1:numel(NumberOfLoopsPerSlice)
   
   csi_k_fill(:,:,:,i,:)=Zerofilling_Spectral(csi_k(:,cumisumi_zwero(i)+1:cumisumi(i),:,i,:),[size(csi_k,1) max(NumberOfLoopsPerSlice) size(csi_k,3) 1 size(csi_k,5)]);
   image_k_fill(:,:,:,i,:)=Zerofilling_Spectral(image_k(:,cumisumi_zwero(i)+1:cumisumi(i),:,i,:),[size(image_k,1) max(NumberOfLoopsPerSlice) size(image_k,3) 1 size(image_k,5)]);
   noise_sim_fill(:,:,:,i,:)=Zerofilling_Spectral(noise_sim(:,cumisumi_zwero(i)+1:cumisumi(i),:,i,:),[size(noise_sim,1) max(NumberOfLoopsPerSlice) size(noise_sim,3) 1 size(noise_sim,5)]);
   
end

csi_k=csi_k_fill;
image_k=image_k_fill;
noise_sim=noise_sim_fill;

clear csi_k_fill noise_sim_fill image_k_fill
% 
% csi_k_scaled=csi_k;
% image_k_scaled=image_k;
% noise_sim_scaled=noise_sim;