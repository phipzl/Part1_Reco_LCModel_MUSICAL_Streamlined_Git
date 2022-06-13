pattern = [1:100];sl_n_csi_encodes=100;
pseudoRadius= rand(1,100);
x=pseudoRadius;

for ii = 0:sl_n_csi_encodes-1
	
		for jj = ii+1 :sl_n_csi_encodes
		
			
			if (pseudoRadius(jj) < pseudoRadius(ii+1))
			
				iTemp = pseudoRadius(ii+1);
				pseudoRadius(ii+1) = pseudoRadius(jj);
				pseudoRadius(jj) = iTemp;
				
				jTemp = pattern(ii+1);
				pattern(ii+1) = pattern(jj);
				pattern(jj) = jTemp;
				
				
              end
         end


end

3