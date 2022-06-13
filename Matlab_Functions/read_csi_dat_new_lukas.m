function [kSpace, Info] = read_csi_dat_new(file)

    global launchtrackpoints;
    kSpace.ONLINE=0;
    kSpace.PATREFSCAN=0;
    kSpace.NOISEADJSCAN=0;
    

    inData=mapVBVD(file);
    
    
    
    if(size(inData,2)==2)
        inData=inData{2};
    end

    if( isfield(inData,'RTfeedback'))
       kSpace.SPINECHO=0; 
    end
    
    Info.Par.PhaseLines =  inData.hdr.Config.PhaseEncodingLines;
    Info.Par.PhaseFoV =  inData.hdr.Config.PhaseFoV;
    Info.Par.PhaseVOI =  inData.hdr.Config.PhaseFoV;
    Info.Par.ReadLines =  inData.hdr.Config.PhaseEncodingLines;
    Info.Par.ReadFOV =  inData.hdr.Config.PhaseFoV;
    Info.Par.ReadVOI =  inData.hdr.Config.PhaseFoV;
    Info.Par.SliceLines =  inData.hdr.Spice.Partitions;
    Info.Par.SliceFOV =  inData.hdr.MeasYaps.sSliceArray.asSlice{1}.dThickness;
    Info.Par.SliceVOI =  inData.hdr.Phoenix.sSpecPara.sVoI.dThickness;

    Info.Par.VoiPositionCor =  inData.image.slicePos(2,1);
    Info.Par.VoiPositionSag =  inData.image.slicePos(1,1);
    Info.Par.VoiPositionTra =  inData.image.slicePos(3,1);
    Info.Par.Dwelltime = inData.hdr.MeasYaps.sRXSPEC.alDwellTime{1};
    Info.Par.nAve =  inData.hdr.Config.NAve;
    Info.Par.PatientName = inData.hdr.Config.PatientName;
    Info.Par.LarmorFreq = inData.hdr.Dicom.lFrequency;
    %  = numel(inData.hdr.MeasYaps.sCoilSelectMeas.aRxCoilSelectData{1}.asList);
    Info.Par.grad_corr_flag=0;
    %  run Calc_Traj_and_Eddy.m
    
    if(strcmpi(inData.image.softwareVersion,'vd'))
         Info.Par.nCha = numel(inData.hdr.MeasYaps.sCoilSelectMeas.aRxCoilSelectData{1}.asList);
    end

    if(strcmpi(inData.image.softwareVersion,'vb'))
        Info.Par.nCha = numel(inData.hdr.MeasYaps.asCoilSelectMeas{1}.asList);
    end


    for CurrentMeasSet2 = transpose(fields(kSpace))
        CurrentMeasSet = CurrentMeasSet2{1};  
       
        fprintf('\nRead\t%s data.', CurrentMeasSet)
        if(strcmp(CurrentMeasSet2,'ONLINE'))
            tic;inImage=feval('single',inData.image());
            
        elseif(strcmp(CurrentMeasSet2,'PATREFSCAN'))
            tic;inImage=feval('single',inData.refscan());
        elseif(strcmp(CurrentMeasSet2,'SPINECHO')) %=RTFEEDBACK
            tic;inImage=feval('single',inData.RTfeedback());   
%             inImage=inImage(:,:,:,:,:,:,:,:,:,:,2,:,:);
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
        nPointsPerRing = 2*inData.image.NPointsPerCircle;
        nADCmax=size(inImage,5);
        nPointsPerADC = size(inImage,1);
        nVec= round(nPointsPerADC*nADCmax/nPointsPerRing-0.5);
        nTImax=size(inImage,6);
%         
%         Info.Data(7,:)=unique(inData.image.Lin,'stable'); % we also need this later before gridding!!!
%         inImage=inImage(:,:,  Info.Data(7,:)  ,:,:,:,:,:);
%         
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
        if(strcmp(CurrentMeasSet2,'ONLINE') || strcmp(CurrentMeasSet2,'SPINECHO'))
            counter = 2;
            counter1 = 1;
            temp=inData.image.iceParam;
            Info.Data=NaN(6,nRings);
            Info.Data(1:4,1)=temp(1:4,1);
            Info.Data(6,1)=1;

            for n=2:size(temp,2)
                if(inData.image.Lin(n)~=inData.image.Lin(n-1) )
                    Info.Data(1:4,counter)=temp(1:4,n);
                    if(temp(1,n)<temp(1,n-1))
                        counter1 = counter1 + 1;
                    end
                    Info.Data(6,counter)=counter1;
                    counter = counter + 1;
                end
            end
            clear temp counter counter 1

            
            [brauchmaned,brauchma]=sort(unique(inData.image.Lin,'stable'));
            Info.Data=Info.Data(:,brauchma);
            
            
            for n=1:nRings
                if(Info.Data(3,n)>60000)
                    Info.Data(3,n) = Info.Data(3,n) -65536;
                end
                if(Info.Data(4,n)>20000)
                    Info.Data(4,n) = Info.Data(4,n) -65536;
                end
            end
        end

    %     if(0)
    %         fprintf('\nThis is constant TI data\n');
    %         nTI=size(inImage,6);
    %         inImage=permute(inImage,[2 3 4 6 1 5]);
    %         inImage=inImage(:,:,:,:,1:nVec*nPointsPerRing);
    % 
    %         inImage=reshape(inImage,[nCha nRings nAve nTI nPointsPerRing nVec]);
    %         kSpace.(CurrentMeasSet)=zeros(nCha,nRings,nPointsPerRing, 1,1,nVec*nTI,nAve,'single');
    %         for iTI=1:nTI
    %             kSpace.(CurrentMeasSet)(:,:,:,1,1,iTI:nTI:end,:)=permute(inImage(:,:,:,iTI,:,:),[1 2 5 4 6 3]);
    %         end
    % %        clear inImage
    %        Info.Data(5,:)=nTI;
    %     end 

    %     if(isvarTI)
    %         fprintf('\nThis is variable TI data\n');
    %         
    %         inImage=permute(inImage,[2 3 4 6 1 5]);
    %         inImage=reshape(inImage,[size(inImage,1) size(inImage,2) size(inImage,3) size(inImage,4) size(inImage,5)*size(inImage,6)]);
    %         nVec = floor(nVec  / factorial(nTImax)) * factorial(nTImax);
    %         kSpace.(CurrentMeasSet)=zeros(nCha,nRings,nPointsPerRing,1,1,nVec,nAve,'single');
    %         
    %         for iRing=1:nRings
    %                        
    %             temp=inImage(:,iRing,:,:,:);
    %             nTI=sum(temp(1,1,1,:,1,1)~=0);
    %             Info.Data(5,iRing)=nTI;
    %             temp=temp(:,:,:,1:nTI,:);
    %             TakeOnlyPoints=nVec/nTI*nPointsPerRing;
    %             temp(temp==0)=[];
    % 
    %             temp=reshape(temp,[size(inImage,1) 1 1 nTI size(inImage,5)*size(inImage,6)/nTI]);
    % 
    %             temp=temp(:,:,:,:,1:TakeOnlyPoints);
    % 
    %             temp=reshape(temp,[size(inImage,1) 1 nTI nPointsPerRing nVec/nTI]);
    %             temp=permute(temp,[1 2 4 3 5]);
    %             temp=reshape(temp,[size(inImage,1) 1 nPointsPerRing nVec]);
    % 
    %             kSpace.(CurrentMeasSet)(:,iRing,:,1,1,:,:) = temp(:,1,:,:);
    %             
    %         end
    %     end
        
          
        if(isvarTI)
            fprintf('\nThis is variable TI data\n');
            if(nVec>factorial(nTImax))
                nVec = floor(nVec  / factorial(nTImax)) * factorial(nTImax); 
            else
                nVec = floor(nVec);
            end
        else 
            fprintf('\nThis is constant TI data\n');
        end


        kSpace.(CurrentMeasSet)=zeros(nCha,nRings,nPointsPerRing,1,1,nVec,nAve,'single');

        for iCha = 1:nCha
            fprintf('Channel %d of %d\n',iCha,nCha);
            inImage_temp = inImage(:,iCha,:,:,:,:); 
            inImage_temp=permute(inImage_temp,[2 3 4 6 1 5]);
            inImage_temp=reshape(inImage_temp,[size(inImage_temp,1) size(inImage_temp,2) size(inImage_temp,3) size(inImage_temp,4) size(inImage_temp,5)*size(inImage_temp,6)]);
            
%             lukas
%             [brauchmaned,brauchma]=sort(Info.Data(7,:));
%             inImage_temp=inImage_temp(:,brauchma,:,:,:,:);
       
%         prevlen = 0;
            for iRing=1:nRings
        %             msg = sprintf('iRing %d/%d of Cha %d\n',iRing,nRings,iCha);
        %             fprintf([repmat('\b',1,prevlen) '%s'],msg);
        %             prevlen = numel(msg);

                temp=inImage_temp(:,iRing,:,:,:);
                nTI=sum(temp(1,1,1,:,1,1)~=0);
                Info.Data(5,iRing)=nTI;
                         
                temp=temp(:,:,:,1:nTI,:);
                TakeOnlyPoints=nVec/nTI*nPointsPerRing;
                temp=temp(:);
                if(  nTI==2 || nTI==3  || nTI ==4 )
                   temp(numel(temp)/nTI+1:end)=[];
                end

                temp=reshape(temp,[size(inImage_temp,1) 1 size(inImage_temp,3) nTI size(inImage_temp,5)*size(inImage_temp,6)/nTI]);

                temp=temp(:,:,:,:,1:TakeOnlyPoints);

                temp=reshape(temp,[size(inImage_temp,1) size(inImage_temp,3) nTI nPointsPerRing nVec/nTI]);
                temp=permute(temp,[1 2 4 3 5]);
                temp=reshape(temp,[size(inImage_temp,1) size(inImage_temp,3) nPointsPerRing nVec]);
                temp=permute(temp,[1 3 4 2]);

                kSpace.(CurrentMeasSet)(iCha,iRing,:,1,1,:,:) = temp(:,:,:,:);

            end
            
        end
        
        clear inImage_temp inImage
        fprintf('\n\t\t\t\t...took\t%10.2f seconds',toc) 
        
    end
    
    clear inData

end
