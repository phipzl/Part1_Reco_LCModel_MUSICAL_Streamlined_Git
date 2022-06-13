function [csi_k_p,image_k_p,noise_sim_k_p,NumberOfLoopsPerSlice,NoOfTempInterleaves_vec] = Initialize3DRollercoaster(slicenmb,csi_k,image_k,noise_sim,nSlices,ReadInInfo,isStackofRings)

    if(nSlices~=1 )
        for i=-(nSlices-1)/2:(nSlices-1)/2

            NumberOfLoopsPerSlice(1+i+(nSlices-1)/2)=ceil(sqrt(ReadInInfo.Par.ReadLines/2*ReadInInfo.Par.ReadLines/2-i^2*(ReadInInfo.Par.ReadLines/2*ReadInInfo.Par.ReadLines/2)/((nSlices-1)/2)/((nSlices-1)/2)));
            if (NumberOfLoopsPerSlice(1+i+(nSlices-1)/2)==0 ||NumberOfLoopsPerSlice(1+i+(nSlices-1)/2)==1) 
                NumberOfLoopsPerSlice(1+i+(nSlices-1)/2)=2;
            end
            if(isStackofRings)
                NumberOfLoopsPerSlice(1+i+(nSlices-1)/2)=ReadInInfo.Par.ReadLines/2;
            end
            
            
        end
        %inplane radius
       
        
    else

        % only if hammin
        if(size(csi_k,2) <= ReadInInfo.Par.ReadLines/2 )
            NumberOfLoopsPerSlice=ReadInInfo.Par.ReadLines/2;
        else
            NumberOfLoopsPerSlice=size(csi_k,2);
        end
    end

    if(nSlices~=1 && max(NumberOfLoopsPerSlice) <= ReadInInfo.Par.ReadLines/2 )

       cumisumi_zwero=[0 cumsum(NumberOfLoopsPerSlice)];
       cumisumi=cumsum(NumberOfLoopsPerSlice);

       csi_k_p=zeros([size(csi_k,1) max(NumberOfLoopsPerSlice) size(csi_k,3) 1 size(csi_k,5)],'single');
       image_k_p=zeros([size(csi_k,1) max(NumberOfLoopsPerSlice) size(csi_k,3) 1 size(image_k,5)],'single');
       noise_sim_k_p=zeros([size(csi_k,1) max(NumberOfLoopsPerSlice) size(csi_k,3) 1 size(noise_sim,5)],'single');


       csi_k_p(:,:,:,1,:)=Zerofilling_Spectral(csi_k(:,cumisumi_zwero(slicenmb)+1:cumisumi(slicenmb),:,1,:),[size(csi_k,1) max(NumberOfLoopsPerSlice) size(csi_k,3) 1 size(csi_k,5)]);
       image_k_p(:,:,:,1,:)=Zerofilling_Spectral(image_k(:,cumisumi_zwero(slicenmb)+1:cumisumi(slicenmb),:,1,:),[size(image_k,1) max(NumberOfLoopsPerSlice) size(image_k,3) 1 size(image_k,5)]);
       noise_sim_k_p(:,:,:,1,:)=Zerofilling_Spectral(noise_sim(:,cumisumi_zwero(slicenmb)+1:cumisumi(slicenmb),:,1,:),[size(noise_sim,1) max(NumberOfLoopsPerSlice) size(noise_sim,3) 1 size(noise_sim,5)]);
    
    end
    
    if(nSlices==1)
       csi_k_p=csi_k;
       image_k_p=image_k;   
       noise_sim_k_p=noise_sim;
    end
    
    if(nSlices==1)
        bla=0;
    else
        cmsm=cumsum(NumberOfLoopsPerSlice);
        bla=cmsm((nSlices-1)/2);
    end;
    for i=1:NumberOfLoopsPerSlice(slicenmb);
        NoOfTempInterleaves_vec(i)=ReadInInfo.Data(5,i+bla);
    end;

end