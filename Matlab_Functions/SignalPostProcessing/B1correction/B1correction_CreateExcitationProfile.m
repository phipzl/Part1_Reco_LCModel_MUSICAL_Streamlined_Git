
function profile=createExcitationProfile()
%% Pulse preparation (leave in in case this has to be redone in the future)
if false
    pulse=load('PW_Pulse_600us_600pts');        % Load
    filter=hamming(600);                        % Hamming filter
    waveform_filtered=pulse(:,1).*filter;
    pulse_filtered=pulse;
    pulse_filtered(:,1)=waveform_filtered;      
    dlmwrite('PW_Pulse_600us_600pts_HammF', ... % Write to file for pulse wizard
        pulse_filtered,'delimiter','\t')     
                                                % Plots for debugging plot(pulse); plot(pulse_filtered); plot(filter);
 end                                            % At this point, open PulseWizard, load the pulse and run the simulation

%% Convert to flip angle distribution
% 5 cols: freq, Mx, My, Mz, Mxy
excitation=load('PW_Excit_39deg_600us_600pts_HammF_FreqDom');

% calculate flip angle distribution
alpha = asin(excitation(:,5))*180/pi;
% Corresponding frequency scale
Fsc = excitation(:,1);

% plot(Fsc,alpha);

%% Find FWHM
% Beware: FWHM is different for alpha-distribution than for Mxy!
HM = max(alpha)/2;
index1 = find(alpha >= HM, 1, 'first'); % Find where the data rises above half the max.
index2 = find(alpha >= HM, 1, 'last');
FWHMi  = index2-index1 + 1;             % FWHM in indexes.
FWHM   = Fsc(index2) - Fsc(index1);     % FWHM in kHz.

clearvars index1 index2 HM FWHMi

%% Calculate necessary gradient
% Numbers
width = 133;       % mm
gamma = 42.58*1E6; % Hz/T

gradient_Hz_per_cm = 1E3*FWHM/(width/10);        % Hz/cm
gradient_mT_per_m  = (gradient_Hz_per_cm/gamma)*1E5;       % mT/m - Note: 1E2 from Hz/cm to Hz/m, 1E3 from T to mT 

fprintf('\n');
fprintf('FWHM: \t\t\t%3.4f kHz.\n', FWHM);
fprintf('Desired slab width: \t%6i mm.\n', width);
fprintf('Required gradient: \t%.4f mT/m.\n', gradient_mT_per_m);
fprintf('Equivalent to \t\t%.4f kHz/cm.\n', gradient_Hz_per_cm/1000);

profile = [Fsc alpha];


end





































