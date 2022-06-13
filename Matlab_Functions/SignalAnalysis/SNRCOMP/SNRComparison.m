%%
clear all
RootPath = '/net/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/lhingerl/CONCEPT_INVIVO_MAPS/CONCEPT_PAPER';
Vols = {'Vol01','Vol02','Vol03','Vol04','Vol05'};
%Methods = {'Concept_NAT_SNRBugfix','Concept_HAM_SNRBugfix','ePE','ePER5'};
Methods = {'Concept_NAT_SNRBugfix_GRDonHAM2_withPostGRD','Concept_HAM_SNRBugfix_GRDonHAM2_withPostGRD','ePE','ePER5'};
MethodsShort = {'NAT','HAM','ePE','ePER5'};

%%

AllSNRMedians = zeros([numel(Vols) numel(Methods) 64*64]);
for VolNo = 1:numel(Vols)
    for MethodNo = 1:numel(Methods)
        
        %VolNo,MethodNo
        
        load([RootPath '/' Vols{VolNo} '/' Methods{MethodNo} '/SNR_Computations/AllSNRMaps.mat' ])
        load([RootPath '/' Vols{VolNo} '/' Methods{MethodNo} '/maps/AllMaps.mat' ])
        
        MaskShrunk = MaskShrinkOrGrow(AllMaps.Masks.Normal,1,0,1);
        MaskShrunk2 = MaskShrinkOrGrow(AllMaps.Masks.Normal,18,0,1);
        Temp = SNR_spectral_NAA.* (MaskShrunk-0*MaskShrunk2);
        Temp(Temp == 0) = NaN;
%             AllSNRMaps.(Vols{VolNo}).(MethodsShort{MethodNo}) = Temp;
        
     %   AllSNRMedians(VolNo,MethodNo) = nanmedian_own(Temp(:));
      
         AllSNRMedians(VolNo,MethodNo,:) = (Temp(:));
    
    
    end
end


AllSNRMedians(:,3,:)=AllSNRMedians(:,3,:)/sqrt((60*30+7)/(5*60+52));
AllSNRMedians(:,4,:)=AllSNRMedians(:,4,:)/sqrt((6*60+10)/(5*60+52));

% AllSNRMedians

% AllSNRMedians_Mean = mean(AllSNRMedians(1:end,:))
% AllSNRMedians_Std = std(AllSNRMedians(1:end,:))


nanmean(AllSNRMedians,3)'
AllSNRMedians_Std = nanstd(AllSNRMedians,0,3)'

AllSNRMedians_meanofratio_eC_R1 = nanmean( AllSNRMedians(:,1,:) ./  	AllSNRMedians(:,3,:),3)';
AllSNRMedians_meanofratio_dwC_R1 = nanmean( AllSNRMedians(:,2,:) ./  AllSNRMedians(:,3,:),3)';
AllSNRMedians_meanofratio_eC_R5 = nanmean( AllSNRMedians(:,1,:) ./  AllSNRMedians(:,4,:),3)';
AllSNRMedians_meanofratio_dwC_R5 = nanmean( AllSNRMedians(:,2,:) ./  AllSNRMedians(:,4,:),3)';

ec_dw=nanmean( AllSNRMedians(:,2,:) ./  AllSNRMedians(:,1,:),3);
ec_pe=nanmean( AllSNRMedians(:,1,:) ./  AllSNRMedians(:,3,:),3);
ec_r5=nanmean( AllSNRMedians(:,1,:) ./  AllSNRMedians(:,4,:),3);
dw_pe=nanmean( AllSNRMedians(:,2,:) ./  AllSNRMedians(:,3,:),3);
dw_r5=nanmean( AllSNRMedians(:,2,:) ./  AllSNRMedians(:,4,:),3);
pe_r5=nanmean( AllSNRMedians(:,3,:) ./  AllSNRMedians(:,4,:),3);

[ec_pe dw_pe ec_r5 dw_r5 ec_dw pe_r5]'

ec_dw_s=nanstd( AllSNRMedians(:,2,:) ./  AllSNRMedians(:,1,:),0,3);
ec_pe_s=nanstd( AllSNRMedians(:,1,:) ./  AllSNRMedians(:,3,:),0,3);
ec_r5_s=nanstd( AllSNRMedians(:,1,:) ./  AllSNRMedians(:,4,:),0,3);
dw_pe_s=nanstd( AllSNRMedians(:,2,:) ./  AllSNRMedians(:,3,:),0,3);
dw_r5_s=nanstd( AllSNRMedians(:,2,:) ./  AllSNRMedians(:,4,:),0,3);
pe_r5_s=nanstd( AllSNRMedians(:,3,:) ./  AllSNRMedians(:,4,:),0,3);

[ec_pe_s dw_pe_s ec_r5_s dw_r5_s ec_dw_s pe_r5_s]'

ec_dw_MEAN=nanmean(nanmean( AllSNRMedians(:,2,:) ./  AllSNRMedians(:,1,:),3),1);
ec_pe_MEAN=nanmean(nanmean( AllSNRMedians(:,1,:) ./  AllSNRMedians(:,3,:),3),1);
ec_r5_MEAN=nanmean(nanmean( AllSNRMedians(:,1,:) ./  AllSNRMedians(:,4,:),3),1);
dw_pe_MEAN=nanmean(nanmean( AllSNRMedians(:,2,:) ./  AllSNRMedians(:,3,:),3),1);
dw_r5_MEAN=nanmean(nanmean( AllSNRMedians(:,2,:) ./  AllSNRMedians(:,4,:),3),1);
pe_r5_MEAN=nanmean(nanmean( AllSNRMedians(:,3,:) ./  AllSNRMedians(:,4,:),3),1);

[ec_pe_MEAN dw_pe_MEAN ec_r5_MEAN dw_r5_MEAN ec_dw_MEAN pe_r5_MEAN]'

ec_dw_std=std(nanmean( AllSNRMedians(:,2,:) ./  AllSNRMedians(:,1,:),3),1);
ec_pe_std=std(nanmean( AllSNRMedians(:,1,:) ./  AllSNRMedians(:,3,:),3),1);
ec_r5_std=std(nanmean( AllSNRMedians(:,1,:) ./  AllSNRMedians(:,4,:),3),1);
dw_pe_std=std(nanmean( AllSNRMedians(:,2,:) ./  AllSNRMedians(:,3,:),3),1);
dw_r5_std=std(nanmean( AllSNRMedians(:,2,:) ./  AllSNRMedians(:,4,:),3),1);
pe_r5_std=std(nanmean( AllSNRMedians(:,3,:) ./  AllSNRMedians(:,4,:),3),1);

[ec_pe_std dw_pe_std ec_r5_std dw_r5_std ec_dw_std pe_r5_std]'

%[median(AllSNRMedians_meanofratio_eC_R1) median(AllSNRMedians_meanofratio_dwC_R1) median(AllSNRMedians_meanofratio_eC_R5) median(AllSNRMedians_meanofratio_dwC_R5)]

% p_eC_R5=signrank([nanmean(AllSNRMedians(1,1,:),3) nanmean(AllSNRMedians(2,1,:),3) nanmean(AllSNRMedians(3,1,:),3) nanmean(AllSNRMedians(4,1,:),3) nanmean(AllSNRMedians(5,1,:),3)],[nanmean(AllSNRMedians(1,4,:),3) nanmean(AllSNRMedians(2,4,:),3) nanmean(AllSNRMedians(3,4,:),3) nanmean(AllSNRMedians(4,4,:),3) nanmean(AllSNRMedians(5,4,:),3)])
% 
% p_eC_R1=signrank([nanmean(AllSNRMedians(1,1,:),3) nanmean(AllSNRMedians(2,1,:),3) nanmean(AllSNRMedians(3,1,:),3) nanmean(AllSNRMedians(4,1,:),3) nanmean(AllSNRMedians(5,1,:),3)],[nanmean(AllSNRMedians(1,3,:),3) nanmean(AllSNRMedians(2,3,:),3) nanmean(AllSNRMedians(3,3,:),3) nanmean(AllSNRMedians(4,3,:),3) nanmean(AllSNRMedians(5,3,:),3)])
% 
% p_eC_DW=signrank([nanmean(AllSNRMedians(1,1,:),3) nanmean(AllSNRMedians(2,1,:),3) nanmean(AllSNRMedians(3,1,:),3) nanmean(AllSNRMedians(4,1,:),3) nanmean(AllSNRMedians(5,1,:),3)],[nanmean(AllSNRMedians(1,2,:),3) nanmean(AllSNRMedians(2,2,:),3) nanmean(AllSNRMedians(3,2,:),3) nanmean(AllSNRMedians(4,2,:),3) nanmean(AllSNRMedians(5,2,:),3)])
% 
% p_DW_R5=signrank([nanmean(AllSNRMedians(1,2,:),3) nanmean(AllSNRMedians(2,2,:),3) nanmean(AllSNRMedians(3,2,:),3) nanmean(AllSNRMedians(4,2,:),3) nanmean(AllSNRMedians(5,2,:),3)],[nanmean(AllSNRMedians(1,4,:),3) nanmean(AllSNRMedians(2,4,:),3) nanmean(AllSNRMedians(3,4,:),3) nanmean(AllSNRMedians(4,4,:),3) nanmean(AllSNRMedians(5,4,:),3)])
% 
% p_DW_R1=signrank([nanmean(AllSNRMedians(1,2,:),3) nanmean(AllSNRMedians(2,2,:),3) nanmean(AllSNRMedians(3,2,:),3) nanmean(AllSNRMedians(4,2,:),3) nanmean(AllSNRMedians(5,2,:),3)],[nanmean(AllSNRMedians(1,3,:),3) nanmean(AllSNRMedians(2,3,:),3) nanmean(AllSNRMedians(3,3,:),3) nanmean(AllSNRMedians(4,3,:),3) nanmean(AllSNRMedians(5,3,:),3)])
% % 
% AllSNRMedians_Mean_Rel_1 = AllSNRMedians_Mean / AllSNRMedians_Mean(1)
% AllSNRMedians_Mean_Rel_2 = AllSNRMedians_Mean / AllSNRMedians_Mean(2)
% AllSNRMedians_Mean_Rel_3 = AllSNRMedians_Mean / AllSNRMedians_Mean(3)
% AllSNRMedians_Mean_Rel_4 = AllSNRMedians_Mean / AllSNRMedians_Mean(4)

% % Is that correct, or do I need a Gaußian error propagation here?
% AllSNRMedians_Std_Rel_1 = AllSNRMedians_Std / AllSNRMedians_Mean(1)
% AllSNRMedians_Std_Rel_2 = AllSNRMedians_Std / AllSNRMedians_Mean(2)
% AllSNRMedians_Std_Rel_3 = AllSNRMedians_Std / AllSNRMedians_Mean(3)
% AllSNRMedians_Std_Rel_4 = AllSNRMedians_Std / AllSNRMedians_Mean(4)
