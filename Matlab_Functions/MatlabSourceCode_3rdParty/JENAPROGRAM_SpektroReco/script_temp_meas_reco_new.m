function [MetaboFID_corr,MetaboFID_uncorr,WaterFID_corr,WaterFID_uncorr] = script_temp_meas_reco_new(MetaboliteFile,WaterFile,dt)
%
% 

dt = dt/10^9;



% *************************************************************************
% define the file names for the meas and rda data

% [fn,fp]     = uigetfile('*.dat','Select the meas file with metabolite data!','MultiSelect','off');      % meas metabolites
% if ~isstr(fn), return; end
% file_meas_m = strcat(fp,fn); clear fp fn
% 
% [fn,fp]     = uigetfile('*.dat','Select the meas file with water data!','MultiSelect','off');           % meas water
% if ~isstr(fn), return; end
% WaterFile = strcat(fp,fn); clear fp fn
% 
% [fn,fp]     = uigetfile('*.rda','Select the corresponding metabolite rda file!','MultiSelect','off');   % RDA
% if ~isstr(fn), return; end
% RDAFile    = strcat(fp,fn); clear fp fn


%MetaboliteFile = '/media/Daten/Arbeit/PROCESS/Measurement_data/Riechkolben_Studie/N01/stim/meas_133988.dat';
%WaterFile = '/media/Daten/Arbeit/PROCESS/Measurement_data/Riechkolben_Studie/N01/ws/meas_133986.dat';
%RDAFile = '/media/Daten/Arbeit/PROCESS/Measurement_data/Riechkolben_Studie/N01/stim/N01_stimulation_MrSpec.20120118.202239.rda';


% *************************************************************************
% load the meas data (metabolite and water) and header of rda-file 
% Array size: Repetitions x Averages x Channels x vector_Size


data_met        = read_VB15MeasDat(MetaboliteFile);
data_wat        = read_VB15MeasDat(WaterFile);
%[d,h]           = rda_reader_2(0,RDAFile);             clear d RDAFile       % What for????????????????????????????

MetaboFID_uncorr           = squeeze(data_met(1,:,:,:));
N_nex_wat       = size(data_wat,2);
WaterFID_uncorr           = squeeze(mean(data_wat(1,:,:,:),2));   % average over all averages of the water signal.

% *************************************************************************
% determine the number of averages (N_nex), number of coil channels (N_chan) and number of
% sampled FID data points (N_samp)

N_nex           = size(MetaboFID_uncorr,1);
N_chan          = size(MetaboFID_uncorr,2);
N_samp          = size(MetaboFID_uncorr,3);

% *************************************************************************
% calculate the channel weighting factors by analysing the water meas data with SVD method
[dum1, weights]    = chan_combine_3(WaterFID_uncorr);
clear dum1 dum2

% % *************************************************************************
% % modify the meta information in rda header with respect to the size of meas data
% 
% h_wat           = h;
% h_met           = h_wat;
% 
% for k_rda = 1:max(size(h))
%     if      strcmp(h{k_rda,1},'VectorSize')
%         N_samp_temp     = str2num(h{k_rda,2});
%         if N_samp_temp ~= N_samp
%             h_met{k_rda,2}  = num2str(N_samp);
%             h_wat{k_rda,2}  = num2str(N_samp);
%         end
%     elseif  strcmp(h{k_rda,1},'DwellTime')
%         dt              = 0.000001*str2num(h{k_rda,2})
%         ind_dt          = k_rda;
%     elseif  strcmp(h{k_rda,1},'NumberOfAverages')
%         h_wat{k_rda,2}  = num2str(N_nex_wat);
%     elseif  strcmp(h{k_rda,1},'PatientName')
%         name_pat        = h{k_rda,2};
%     elseif  strcmp(h{k_rda,1},'StudyDate')
%         study_date      = h{k_rda,2};
%     elseif  strcmp(h{k_rda,1},'SeriesNumber')
%         series_ID       = h{k_rda,2};
%     end
% end, clear k_rda h
% 
% if N_samp_temp ~= N_samp
%     dt              = 0.5*dt
%     h_met{ind_dt,2} = num2str(round(1000000*dt));
%     h_wat{ind_dt,2} = num2str(round(1000000*dt));
%     clear ind_dt
% end


t               = (0:N_samp-1).*dt;             % absolute time scale

% *************************************************************************
% define the data variables

WaterFID_corr   = zeros(1,N_samp);              % channel combined, phased water data      
MetaboFID_corr = zeros(1,N_samp);              % channel combined, phased and frequency corrected water data            
d_met_ph_comb   = zeros(N_nex,N_samp);          % channel combined, phased water data
d_met_ph        = zeros(N_nex,N_chan,N_samp);   % phased water data from single coil channels

% *************************************************************************
% phase correction and weighted channel combination of water data

for k_chan  = 1:N_chan
    WaterFID_corr   = WaterFID_corr + weights(k_chan).*phase_cor_autom_1(WaterFID_uncorr(k_chan,:)); 
end, clear k_chan

% *************************************************************************
% phase and eddy current correction (according to the phase distortions in water data from corresponding coil channels) and weighted channel combination of metabolite data

for k_nex   = 1:N_nex
    for k_chan  = 1:N_chan
        d_met_ph(k_nex,k_chan,:)    = phase_cor_autom_1(ecc_cor(permute(squeeze(MetaboFID_uncorr(k_nex,k_chan,:)),[2 1]),WaterFID_uncorr(k_chan,:)));    % phase and eddy current correction

        d_met_ph_comb(k_nex,:)      = d_met_ph_comb(k_nex,:) + weights(k_chan).*permute(squeeze(d_met_ph(k_nex,k_chan,:)),[2 1]);   % weighted channel combination
    end, clear k_chan    	
end, clear k_nex 

% *************************************************************************
% frequency correction of single acquisitions of metabolite data


w_dif       = f_shift_H2O_1_2(d_met_ph_comb,70,dt);      % determine the current frequency shift of the water peak relative to the excpected water frequency (at 0 Hz)
% N = size(d_met_ph_comb,2);
% freqvec = (-(N-0.5)/2:1:(N-0.5)/2) .* (1/dt)/N;
% bla = figure; movegui(bla,'northwest')
% bla2 = figure; movegui(bla2,'southwest')
% FreqDown = find(min(abs(freqvec + 30)) == abs(freqvec + 30));
% FreqUp = find(min(abs(freqvec - 30)) == abs(freqvec - 30));


for k_nex = 1:N_nex
	
% 	tempspec = fftshift(fft( d_met_ph_comb(k_nex,:) ,[],2),2);
% 	freqvec(tempspec == max(tempspec))
% 	set(0,'currentfigure',bla), scatter(freqvec(FreqDown:FreqUp), real( tempspec(FreqDown:FreqUp) )), title(['Uncorr, Avg = ' num2str(k_nex)])
	
    MetaboFID_corr(k_nex,:)        = d_met_ph_comb(k_nex,:).* exp(1i*2*pi*w_dif(k_nex)*t);

% 	tempspec = fftshift(fft( MetaboFID_corr(k_nex,:) ,[],2),2);
% 	freqvec(tempspec == max(tempspec))
% 	set(0,'currentfigure',bla2), scatter(freqvec(FreqDown:FreqUp), real( tempspec(FreqDown:FreqUp) )), title(['Corr, Avg = ' num2str(k_nex)])
% 	waitforbuttonpress
	
	
end, clear k_nex

% % *************************************************************************
% % save the reconstructed mean spectra in new rda files
%         
% dir_save        = uigetdir(pwd,'Select directory to save genrated RDA files!');
% 
% if isstr(dir_save)
%     file_rda_m      = strcat(dir_save,47,name_pat,95,'Ser',series_ID,95,study_date,95,'MET.rda');
%     file_rda_w      = strcat(dir_save,47,name_pat,95,'Ser',series_ID,95,study_date,95,'WAT.rda');
% 
%     rda_saver(conj(mean(MetaboFID_corr,1)),convert_header(h_met),file_rda_m);
%     rda_saver(conj(WaterFID_corr),convert_header(h_wat),file_rda_w);
% end

% *************************************************************************
% clear the working variables

clear N_nex N_chan N_samp_temp N_nex_wat weights t ind_dt h_wat dir_save file_rda_m file_rda_w name_pat study_date series_ID data_met data_wat

clear d_met_ph_comb d_met_ph w_dif dt N_samp h_met 

