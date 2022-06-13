%%
clear all
RootPath = '/net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/CONCEPT_INVIVO_MAPS/CONCEPT_PAPER';
Vols = {'Vol01','Vol02','Vol03','Vol04','Vol05'};
Methods = {'Concept_NAT_SNRBugfix','Concept_HAM_SNRBugfix','ePE','ePER5'};
MethodsShort = {'NAT','HAM','ePE','ePER5'};

%%

AllSNRMedians = zeros([numel(Vols) numel(Methods)]);
for VolNo = 1:numel(Vols)
    for MethodNo = 1:numel(Methods)
        
        %VolNo,MethodNo
        
        load([RootPath '/' Vols{VolNo} '/' Methods{MethodNo} '/SNR_Computations/AllSNRMaps.mat' ])
        load([RootPath '/' Vols{VolNo} '/' Methods{MethodNo} '/maps/AllMaps.mat' ])
        
        MaskShrunk = MaskShrinkOrGrow(AllMaps.Masks.Normal,1,0,1);
        Temp = SNR_spectral_NAA.* MaskShrunk;
        Temp(Temp == 0) = NaN;
        AllSNRMaps.(Vols{VolNo}).(MethodsShort{MethodNo}) = Temp;
        
%         AllSNRMedians(VolNo,MethodNo) = nanmedian_own(Temp(:));
      
        AllSNRMedians(VolNo,MethodNo) = nanmean_own(Temp(:));
    end
end
AllSNRMedians(:,3)=AllSNRMedians(:,3)/sqrt((60*30+7)/(5*60+52));
AllSNRMedians(:,4)=AllSNRMedians(:,4)/sqrt((6*60+10)/(5*60+52));

AllSNRMedians

AllSNRMedians_Mean = mean(AllSNRMedians(1:end,:))
AllSNRMedians_Std = std(AllSNRMedians(1:end,:))

AllSNRMedians_Mean_Rel_1 = AllSNRMedians_Mean / AllSNRMedians_Mean(1)
AllSNRMedians_Mean_Rel_2 = AllSNRMedians_Mean / AllSNRMedians_Mean(2)
AllSNRMedians_Mean_Rel_3 = AllSNRMedians_Mean / AllSNRMedians_Mean(3)
AllSNRMedians_Mean_Rel_4 = AllSNRMedians_Mean / AllSNRMedians_Mean(4)

% Is that correct, or do I need a Gaußian error propagation here?
AllSNRMedians_Std_Rel_1 = AllSNRMedians_Std / AllSNRMedians_Mean(1)
AllSNRMedians_Std_Rel_2 = AllSNRMedians_Std / AllSNRMedians_Mean(2)
AllSNRMedians_Std_Rel_3 = AllSNRMedians_Std / AllSNRMedians_Mean(3)
AllSNRMedians_Std_Rel_4 = AllSNRMedians_Std / AllSNRMedians_Mean(4)
