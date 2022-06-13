function OutData=Heidemann2006GrappaFunc(OutData, ACS, UndersamplingFactor,MinKernelSrcPts,NoOfSegments,AdditionalPhiPts) 



%% Definitions

% This code works by dividing the CONCEPT data in segments, and then performs a 1D-Cartesian-GRAPPA in (r,phi)-space along r
% (in principle could be also along phi, but for CONCEPT that's not useful)
% So we divide the 360° e.g. in 20 segments, calculate the weights for each segment, and apply the segments. However, this causes a
% problem: The standard GRAPPA algorithm expands zeros on the edges of the undersampled data to be able to reconstruct all missing
% data points. But since we subdivide our k-space, we assume that at all the edges of all the segments there are zeros, which is
% wrong. Solution: We give more points along phi to the algorithm, and remove all the additional, unneccessary points after the
% reco. With that, we throw away the points where we edge-zero-filled the data, and thus don't care about them!
% This is done with AdditionalPhiPts. Use values larger than half the kernel size along phi (you can check the used 
% kernelsize by the third output of "opencaipirinha_MRSI" (it has 4 entries: 
% left-to-TargetPt, right-to-TargetPt, above-TargetPt, below-TargetPt))

TotPhiPts = size(ACS,3);
if(~exist('NoOfSegments','var'))
    PossibleNoOfSegments = 30:70;
    SegSz = TotPhiPts./PossibleNoOfSegments;
    PossibleNoOfSegments(mod(SegSz,1) ~= 0) = Inf;
    [dum, ind] = min(abs(PossibleNoOfSegments - 56)); %#ok    % Heidemann used 56
    NoOfSegments = PossibleNoOfSegments(ind);
    clear PossibleNoOfSegments SegSz dum ind
end
if(~exist('MinKernelSrcPts','var'))
    MinKernelSrcPts = 20;
end
if(~exist('AdditionalPhiPts','var'))        % These are the points along phi that we additionally 
    AdditionalPhiPts = 4;                   % use on each side of the segment to calculate the weights
end                            

ComplicatedOrSimpleReco = 'Complicated';

if(NoOfSegments == 1)
    ComplicatedOrSimpleReco = 'Simple';    
end


%% More Definitions

if(mod(TotPhiPts,NoOfSegments) ~= 0)
   error(['Can''t divide TotPtsAlongPhi of ' num2str(TotPhiPts) ' into requested segments of ' num2str(NoOfSegments) '.'])
end

%% Undersample OutData along radius

UndersamplingCell = zeros([UndersamplingFactor 1]); UndersamplingCell(1) = 1;
OutData = Skip_kPoints(OutData,UndersamplingCell);


%% Reconstruct data segment-wise using only first time point of ACS

if(isempty(regexpi(ComplicatedOrSimpleReco,'Compl')))
    csi_k_Reco = zeros(size(OutData));
    for CurSeg = 1:NoOfSegments

        %CurSeg
        PhiPtsToReco = [(CurSeg-1)*TotPhiPts/NoOfSegments+1 (CurSeg)*TotPhiPts/NoOfSegments]; % These Phi pts we want to reconstruct now
        PhiPtsForACS = [(CurSeg-1)*TotPhiPts/NoOfSegments+1 (CurSeg)*TotPhiPts/NoOfSegments]; % These Phi pts we use for the ACS calib


        PhiPtsForACS(1) = PhiPtsForACS(1) - AdditionalPhiPts;
        PhiPtsForACS(2) = PhiPtsForACS(2) + AdditionalPhiPts;

        % If we ask for too many points in beginning of segment, append the end of the data to it, because phi is an angle and the
        % start point at phi = 0 is artificial: The data at phi = -Deltaphi should be also very similar to phi = 0
        % same if we ask for too many points at end of segment
        if(PhiPtsForACS(1) <= 0)
            csi_k_skipped_Temp = cat(3,OutData(:,:,end:-1:end+PhiPtsForACS(1),:,:),OutData(:,:,1:PhiPtsForACS(2),:,:));
            ACSTemp = cat(3,ACS(:,:,end:-1:end+PhiPtsForACS(1),:,1),ACS(:,:,1:PhiPtsForACS(2),:,1));
        elseif(PhiPtsForACS(2) > TotPhiPts)
            csi_k_skipped_Temp=cat(3,OutData(:,:,PhiPtsForACS(1):end,:,:),OutData(:,:,1:(PhiPtsForACS(2)-TotPhiPts),:,:));
            ACSTemp=cat(3,ACS(:,:,PhiPtsForACS(1):end,:,1),ACS(:,:,1:(PhiPtsForACS(2)-TotPhiPts),:,1));
        else
            csi_k_skipped_Temp = OutData(:,:,PhiPtsForACS(1):PhiPtsForACS(2),:,:);
            ACSTemp = ACS(:,:,PhiPtsForACS(1):PhiPtsForACS(2),:,1);
        end
        
%         csi_k_skipped_Temp(:,41:43,1:end,:,:)=conj(csi_k_skipped_Temp(:,4:-1:2,end:-1:1,:,:));
%         csi_k_skipped_Temp=circshift(csi_k_skipped_Temp,[0 3 0 0 0]);
%         ACSTemp(:,41:43,1:end,:,:)=conj(ACSTemp(:,4:-1:2,end:-1:1,:,:));
%         ACSTemp=circshift(ACSTemp,[0 3 0 0 0]);
    
        [csi_k_skipped_Temp,weights2DCaip,kernelsize2DCaip,SrcRelativeTarg2DCaip] = opencaipirinha_MRSI(...
        csi_k_skipped_Temp, ACSTemp,UndersamplingCell,true,'double',MinKernelSrcPts); 
        weights{CurSeg} = weights2DCaip;

        csi_k_Reco(:,:,PhiPtsToReco(1):PhiPtsToReco(2),:,:) = csi_k_skipped_Temp(:,:,AdditionalPhiPts+1:end-AdditionalPhiPts,:,:);

    end

end



%% Complicated Reco

if(~isempty(regexpi(ComplicatedOrSimpleReco,'Compl')))
    csi_k_Reco = zeros(size(OutData));
    for CurSeg = 1:NoOfSegments

       % CurSeg
        PhiPtsToReco = [(CurSeg-1)*TotPhiPts/NoOfSegments+1 (CurSeg)*TotPhiPts/NoOfSegments]; % These Phi pts we want to reconstruct now
        PhiPtsForACS = [(CurSeg-1)*TotPhiPts/NoOfSegments+1 (CurSeg)*TotPhiPts/NoOfSegments]; % These Phi pts we use for the ACS calib


        PhiPtsForACS(1) = PhiPtsForACS(1) - AdditionalPhiPts;
        PhiPtsForACS(2) = PhiPtsForACS(2) + AdditionalPhiPts;

        % If we ask for too many points in beginning of segment, append the end of the data to it, because phi is an angle and the
        % start point at phi = 0 is artificial: The data at phi = -Deltaphi should be also very similar to phi = 0
        % same if we ask for too many points at end of segment
        if(PhiPtsForACS(1) <= 0)
            csi_k_skipped_Temp = cat(3,OutData(:,:,end:-1:end+PhiPtsForACS(1),:,1),OutData(:,:,1:PhiPtsForACS(2),:,1));
            ACSTemp = cat(3,ACS(:,:,end:-1:end+PhiPtsForACS(1),:,1),ACS(:,:,1:PhiPtsForACS(2),:,1));
        elseif(PhiPtsForACS(2) > TotPhiPts)
            csi_k_skipped_Temp=cat(3,OutData(:,:,PhiPtsForACS(1):end,:,1),OutData(:,:,1:(PhiPtsForACS(2)-TotPhiPts),:,1));
            ACSTemp=cat(3,ACS(:,:,PhiPtsForACS(1):end,:,1),ACS(:,:,1:(PhiPtsForACS(2)-TotPhiPts),:,1));
        else
            csi_k_skipped_Temp = OutData(:,:,PhiPtsForACS(1):PhiPtsForACS(2),:,1);
            ACSTemp = ACS(:,:,PhiPtsForACS(1):PhiPtsForACS(2),:,1);
        end


        [csi_k_skipped_Temp,weights2DCaip,kernelsize2DCaip,SrcRelativeTarg2DCaip] = opencaipirinha_MRSI(...
        csi_k_skipped_Temp, ACSTemp,UndersamplingCell,true,'double',MinKernelSrcPts); 
        weights{CurSeg} = weights2DCaip;

%         csi_k_Reco(:,:,PhiPtsToReco(1):PhiPtsToReco(2),:,1) = csi_k_skipped_Temp(:,:,AdditionalPhiPts+1:end-AdditionalPhiPts,:,:);

    end

    csi_k_Reco = zeros(size(OutData));
    weights = cat(2,weights(end), weights, weights(1));
    for CurPhiPt = 1:TotPhiPts

       % CurPhiPt
        SegSz = TotPhiPts/NoOfSegments;
        CurSeg = ceil(CurPhiPt/SegSz);

        PhiPtsForACS = [(CurSeg-1)*TotPhiPts/NoOfSegments+1 (CurSeg)*TotPhiPts/NoOfSegments]; % These Phi pts we use for the ACS calib


        PhiPtsForACS(1) = PhiPtsForACS(1) - AdditionalPhiPts;
        PhiPtsForACS(2) = PhiPtsForACS(2) + AdditionalPhiPts;

        % If we ask for too many points in beginning of segment, append the end of the data to it, because phi is an angle and the
        % start point at phi = 0 is artificial: The data at phi = -Deltaphi should be also very similar to phi = 0
        % same if we ask for too many points at end of segment
        if(PhiPtsForACS(1) <= 0)
            csi_k_skipped_Temp = cat(3,OutData(:,:,end:-1:end+PhiPtsForACS(1),:,:),OutData(:,:,1:PhiPtsForACS(2),:,:));
            ACSTemp = cat(3,ACS(:,:,end:-1:end+PhiPtsForACS(1),:,:),ACS(:,:,1:PhiPtsForACS(2),:,:));
        elseif(PhiPtsForACS(2) > TotPhiPts)
            csi_k_skipped_Temp=cat(3,OutData(:,:,PhiPtsForACS(1):end,:,:),OutData(:,:,1:(PhiPtsForACS(2)-TotPhiPts),:,:));
            ACSTemp=cat(3,ACS(:,:,PhiPtsForACS(1):end,:,:),ACS(:,:,1:(PhiPtsForACS(2)-TotPhiPts),:,:));
        else
            csi_k_skipped_Temp = OutData(:,:,PhiPtsForACS(1):PhiPtsForACS(2),:,:);
            ACSTemp = ACS(:,:,PhiPtsForACS(1):PhiPtsForACS(2),:,1);
        end

        SegCenters = -SegSz : SegSz : TotPhiPts;
        SegCenters = SegCenters + (SegSz+1)/2;
        CurDist = abs(CurPhiPt - SegCenters);

        % Take smallest and 2nd smallest distance
        [CurDist,idx]=sort(CurDist);
        if(CurDist(1) > 0)
            NormalizeFac = 1/(1/CurDist(1) + 1/CurDist(2));
            CurWeight = {(weights{idx(1)}{1}/CurDist(1) + weights{idx(2)}{1}/CurDist(2)) * NormalizeFac};
        else
            CurWeight = weights{idx(1)};
        end

        csi_k_skipped_Temp = opencaipirinha_MRSI(...
        csi_k_skipped_Temp, ACSTemp,UndersamplingCell,true,'double',MinKernelSrcPts,CurWeight,kernelsize2DCaip,SrcRelativeTarg2DCaip); 

        csi_k_Reco(:,:,CurPhiPt,:,:) = csi_k_skipped_Temp(:,:,AdditionalPhiPts+mod(CurPhiPt-1,SegSz)+1,:,:);
    end

    weights(end) = []; weights(1) = [];

end



%% Result

OutData = csi_k_Reco;

