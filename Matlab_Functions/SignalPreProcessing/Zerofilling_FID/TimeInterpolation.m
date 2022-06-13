function csiI=TimeInterpolation(csi, truncatefactor, zerofillfactor, FillToOrig)
%
% This function performs ZIP (zero interpolation, zero filling) of CSI data in the time-domain.
% Created by Philipp Lazen, October 2018.
%
% csi 			input: csi data
% csiI 			output: interpolated csi data
% truncatefactor	How much of the signal should be left after truncation.	
% zerofillfactor	How far the zero interpolation should be performed, relative to the length after truncation.
%
%
% As of 2018/10/30 this function uses the flag "TimeInterpolation_flag" and the option "TimeInterpolationFactor" in various bash scripts:
%
% Part1_ProcessMRSI.sh
% 	Lines 149 ("FLAGS", "optional"), 230 ("INITIALIZING", "optional"), 271 (documentation). Also:
%	while getopts 'x:y:z: ... :T' OPTION
%
% write_InitialParameters.sh
%	Lines 28 ("Flags"), 138 ("Additional User Input")
%
% The function itself is so far only called by MRSI_Reconstruction.m.
%
% MRSI_Reconstruction.m
% 	Lines 1050ff (See: "17. Perform Interpolation", "In time domain"), 1155ff ("Interpolation Time Domain)
%
% 


%	display([ char(10) 'function TimeInterpolation()...' char(10)])


	if(FillToOrig)
		OrigLength = size(csi, 4);
	end
		
 	if(truncatefactor ~= 1)
		Size = size(csi);	
		csi  = csi(:,:,:,1:floor(Size(4)*truncatefactor));
 	end

	if(FillToOrig == 1)
        Size = size(csi);
		csiI = zeros(Size(1),Size(2),Size(3),OrigLength);
		csiI(:,:,:,1:size(csi,4)) = csi;
	elseif(zerofillfactor == 1)
		csiI = csi;
    else
        Size = size(csi);
		csiI = zeros(Size(1),Size(2),Size(3),ceil(Size(4)*zerofillfactor));  
		csiI(:,:,:,1:Size(4)) = csi;
	end
end
