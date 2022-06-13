CurPoint.(CurrentMeasSet) = 0; 
	for i = 1 : Info.(CurrentMeasSet).total_channel_no : size(Info.(CurrentMeasSet).mdhInfo,2)
		% {echo}[cha x kx x ky x kz x slc x samples x avg x rep x ADCNo x TempIntNo]
		CurEco = Info.(CurrentMeasSet).mdhInfo(6,i); Curkx = Info.(CurrentMeasSet).mdhInfo(2,i); Curky = Info.(CurrentMeasSet).mdhInfo(3,i); Curkz = Info.(CurrentMeasSet).mdhInfo(4,i);
		CurSlc = Info.(CurrentMeasSet).mdhInfo(5,i); CurAvg = Info.(CurrentMeasSet).mdhInfo(7,i); CurRep = Info.(CurrentMeasSet).mdhInfo(8,i);
		CurADC = Info.(CurrentMeasSet).mdhInfo(11,i); CurTempIntNo = Info.(CurrentMeasSet).mdhInfo(12,i);
		
		Temp{CurEco}(:,Curkx, Curky, Curkz, CurSlc, :, CurAvg, CurRep, CurADC,CurTempIntNo) = ...
		kSpace.(CurrentMeasSet)(:,CurPoint.(CurrentMeasSet)+1:CurPoint.(CurrentMeasSet)+Info.(CurrentMeasSet).mdhInfo(9,i)-Info.(CurrentMeasSet).mdhInfo(10,i));
	
		mdhInfo_reshaped{CurEco}(1,Curkx, Curky, Curkz, CurSlc, 1, CurAvg, CurRep, CurADC,CurTempIntNo,:) = Info.(CurrentMeasSet).mdhInfo(:,i);
	
		CurPoint.(CurrentMeasSet) = CurPoint.(CurrentMeasSet) + (Info.(CurrentMeasSet).mdhInfo(9,i)-Info.(CurrentMeasSet).mdhInfo(10,i));
	end
	Info.(CurrentMeasSet).mdhInfo = mdhInfo_reshaped; clear mdhInfo_reshaped; 
	
	for echo = minecho:maxi(6)
		% CONCEPT RESIZING. ACTUALLY, SHOULD WE DO THIS LATER, SO THAT WE FIRST GET THE REAL RAW DATA AS MEASURED, AND THEN JUST RESHAPE IT ACC. TO OUR NEEDS?
		if(isfieldRecursive(Info,'General','Ascconv','WipMemBlockInterpretation','Rollercoaster','sNoADCPointsPerCircle') && ~strcmpi(CurrentMeasSet,'NoiseAdjScan'))
			Temp{echo} = permute(Temp{echo},[1 2 3 4 5 10 6 9 7 8]);
			Temp{echo} = reshape(Temp{echo},[size_MultiDims(Temp{echo},1:6) size(Temp{echo},7)*size(Temp{echo},8) size_MultiDims(Temp{echo},9:10)]);
			
            if(strcmpi(CurrentMeasSet,'ONLINE'))
                vecSize = Info.General.Ascconv.vecSize/2;							% The system thinks we do oversampling in spectral dimension, which we dont...
            else
                vecSize = floor(Info.PATREFSCAN.Samples / Info.General.Ascconv.WipMemBlockInterpretation.Rollercoaster.sNoADCPointsPerCircle);		% For the pre-scan. This is guess-work! Make better in future!      
            end
			if( isfieldRecursive(Info,'General','Ascconv','WipMemBlockInterpretation','Rollercoaster','sNoADCPointsPerCircle') && Info.General.Ascconv.WipMemBlockInterpretation.Rollercoaster.sNoADCPointsPerCircle > 0 )
				NoOfPtsPerLoop = Info.General.Ascconv.WipMemBlockInterpretation.Rollercoaster.sNoADCPointsPerCircle*2*1;				% *last entry ... Spatial oversampling. always *2 internal
			else	% otherwise make a guess 
				NoOfPtsPerLoop = size(Temp{echo},7) / vecSize;
				NoOfPtsPerLoop = 143*2;
			end
			NoOfTempInterleaves = maxi(12);										% Same here
			TakeOnlyPoints = (vecSize/NoOfTempInterleaves)*NoOfPtsPerLoop;		% vecSize/NoOfTempInterleaves ... Spectral Points, NoOfPtsPerLoop ... PointsPerLoop
			Temp{echo} = Temp{echo}(:,:,:,:,:,:,1:TakeOnlyPoints,:,:);
			Temp{echo} = reshape(Temp{echo},[size_MultiDims(Temp{echo},1:5) NoOfTempInterleaves NoOfPtsPerLoop vecSize/NoOfTempInterleaves size_MultiDims(Temp{echo},8:9)]);
			Temp{echo} = permute(Temp{echo},[1 2 3 4 5 7 6 8 9 10]);
			Temp{echo} = reshape(Temp{echo},[size_MultiDims(Temp{echo},1:5) NoOfPtsPerLoop vecSize size_MultiDims(Temp{echo},9:10)]);
			% Permute dimensions to make size {echo}(cha x LoopNo x PointsPerLoops x kz x slice x samples x avg x rep)
			Temp{echo} = permute(Temp{echo},[1 2 6 4 5 7 8 9 3]);
			
		else
			if(size(Temp{echo},3) == 1 && (size(Temp{echo},2) > 1 || size(Temp{echo},4) > 1) )			% This should basically mean: If data is imaging data. However it is a little cheated
				Temp{echo} = permute(Temp{echo},[1 6 2 4 5 3 7 8 9 10]);
			end
		end
	end
	
	kSpace.(CurrentMeasSet) = Temp;		