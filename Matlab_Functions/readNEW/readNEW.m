kSpace.ONLINE=0;
kSpace.PATREFSCAN=0;
kSpace.NOISEADJSCAN=0;

% inData=mapVBVD('/ceph/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/pmoser/mosaic_move/raw/meas_MID657_bs_Rollercoaster_hack_v4_noMOVE_VC_FID75030.dat');
% inData=mapVBVD('/ceph/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/pmoser/PI_paper/vol01/raw/meas_MID333_bs_Ro11_FID17674.dat');
% inData = mapVBVD('/ceph/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/pmoser/PI_paper/3D_PI/raw/meas_MID291_pmos_Rollercoaster_3DPI_FID21344.dat');
% inData=mapVBVD('/ceph/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Measurement_Data/ComparisonRollerVsSpiral/AbstractMeasNew/tumor4/meas_MID71_3D_MRSI_RVAR9_64x64x39_TR280_FID6934.dat');
inData=mapVBVD('/ceph/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/pmoser/prisma/newnew/meas_MID00094_FID01448_body_60er_varTI.dat');
% inData=mapVBVD('/ceph/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/pmoser/prisma/newnew/meas_MID00095_FID01449_body_60er_constTI.dat');


if(size(inData,2)==2)
    inData=inData{2};
end
     
for CurrentMeasSet2 = transpose(fields(kSpace))
    CurrentMeasSet = CurrentMeasSet2{1};  
    
    fprintf('\nRead\t%s data.', CurrentMeasSet)
    if(strcmp(CurrentMeasSet2,'ONLINE'))
        tic;inImage=feval('single',inData.image());
    elseif(strcmp(CurrentMeasSet2,'PATREFSCAN'))
        tic;inImage=feval('single',inData.refscan());
    elseif(strcmp(CurrentMeasSet2,'NOISEADJSCAN'))
         kSpace.(CurrentMeasSet)=permute(feval('single',inData.noise()),[2 1 3]);
        continue;
    end
    fprintf('\n\t\t\t\t...took\t%10.2f seconds',toc) 
    
    inImage=squeeze_single_dim(inImage,4);
    inImage=squeeze_single_dim(inImage,4);
    inImage=squeeze_single_dim(inImage,4);
    inImage=squeeze_single_dim(inImage,4);
    inImage=squeeze_single_dim(inImage,4);
    inImage=squeeze_single_dim(inImage,4);
    inImage=squeeze_single_dim(inImage,5);


    nCha = size(inImage,2);
    nAve = size(inImage,4);
    nRings = size(inImage,3);
    nPointsPerRing = inData.image.NPointsPerCircle;
    nADCmax=size(inImage,5);
    nPointsPerADC = size(inImage,1);
    nVec= round(nPointsPerADC*nADCmax/nPointsPerRing-0.5);
    nTImax=size(inImage,6);


    temp=inImage(1,1,:,:,:,:);
    isvarTI=0;
    isconstTI=0;
    if(sum(temp(:)==0)>0)
        isvarTI=1;
    else
        isconstTI=1;
    end
    clear temp;

    tic
    fprintf('\nReshape\t%s data.', CurrentMeasSet)
    if(strcmp(CurrentMeasSet2,'ONLINE'))
        counter = 2;
        counter1 = 1;
        temp=inData.image.iceParam;
        Info=NaN(6,nRings);
        Info(1:4,1)=temp(1:4,1);
        Info(6,1)=1;
        for n=2:size(temp,2)
            if(temp(1,n)~=temp(1,n-1))
                Info(1:4,counter)=temp(1:4,n);
                if(temp(1,n)<temp(1,n-1))
                    counter1 = counter1 + 1;
                end
                Info(6,counter)=counter1;
                counter = counter + 1;
            end
        end
        clear temp counter counter 1
        Info(4,:)=Info(4,:)-65536;   
    end
    
    if(isconstTI)
        fprintf('\nThis is constant TI data\n');
        nTI=size(inImage,6);
        inImage=permute(inImage,[2 3 4 6 1 5]);
        inImage=inImage(:,:,:,:,1:nVec*nPointsPerRing);

        inImage=reshape(inImage,[nCha nRings nAve nTI nPointsPerRing nVec]);
        kSpace.(CurrentMeasSet)=zeros(nCha,nRings,nPointsPerRing, 1,1,nVec*nTI,nAve,'single');
        for iTI=1:nTI
            kSpace.(CurrentMeasSet)(:,:,:,1,1,iTI:nTI:end,:)=permute(inImage(:,:,:,iTI,:,:),[1 2 5 4 6 3]);
        end
       clear inImage
       Info(5,:)=nTI;
    end 

    if(isvarTI)
        fprintf('\nThis is variable TI data\n');
        kSpace.(CurrentMeasSet)=zeros(nCha,nRings,1,1,nAve,nPointsPerRing*nVec,'single');
        for iRing=1:nRings
            temp=inImage(:,:,iRing,:,:,:);
            nTI=sum(temp(1,1,1,1,1,:)~=0);
            Info(5,iRing)=nTI;
            for iTI=1:nTI
                temp1=temp(:,:,:,:,:,iTI);
                temp1=permute(temp1,[2 3 4 6 1 5]);
%                 temp1=temp1(:,:,:,:,:);temp1=temp1(:,:,:,:,1:size(temp1,5)/nTI);
%                 kSpace.(CurrentMeasSet)(:,iRing,1,1,:,iTI:nTI:end)=temp1;
                temp1=temp1(:,:,:,:,:);
                kSpace.(CurrentMeasSet)(:,iRing,1,1,:,iTI:nTI:end)=temp1(:,:,:,:,1:nPointsPerRing*nVec/nTI);
            end
        end
        clear temp1 temp inImage
        kSpace.(CurrentMeasSet)=reshape(kSpace.(CurrentMeasSet),[nCha nRings 1 1 nAve nPointsPerRing nVec]);
        kSpace.(CurrentMeasSet)=permute(kSpace.(CurrentMeasSet),[1 2 6 3 4 7 5]);
    end
    
    fprintf('\n\t\t\t\t...took\t%10.2f seconds',toc) 
    

end





