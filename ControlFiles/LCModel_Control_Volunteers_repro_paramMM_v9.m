% Referencing, see LCModel Manual page 105
ControlWrite.DOREFS = {'DOREFS(1) = T','DOREFS(2) = T'};           % T: Use standard water referencing, F: No Other metabolites used for referencing
% ControlWrite.NREFPK = {'NREFPK(2) = 1'};                         % NREFPK(JCCF): number of delta functions used for Reference JCCF, e.g. 1 delta function for NAA at 2.01; Each delta function is numbered by JREF
 ControlWrite.PPMREF = {'PPMREF(2,1) = 3.22'};                       % PPMREF(JREF,JCCF): What chemical shift should the delta function with number JREF for the Reference with number JCCF have?
% ControlWrite.SHIFMNMX = {'SHIFMN(2) = -0.1','SHIFMX(2) = 0.1'};  % SHIFMN(JCCF): Minimum chemshift to search for peak. LCModel searches from PPMREF(2,1) + SHIFMN(2) to PPMREF(2,1) + SHIFMX(2)
                                                                  % SHIFMN should be negative. 

                                                                  
% Water Scaling, Absolute Quantification
ControlWrite.WSMET = 'WSMET = ''DSS''';                          % This tells LCModel what to use for scaling the absolute fitting concentrations. 
ControlWrite.WSPPM = 'WSPPM = 0.0';                              % The chemical shift of WSMET   
ControlWrite.N1HMET = 'N1HMET = 9';                              % The number of protons contributing to the signal



% Plotting Parameters
ControlWrite.SUBBAS =  'SUBBAS = T';                             % Subtracts the baseline from the spectra
ControlWrite.NEACH =  'NEACH = 99';                              % "the number of metabolites for which individual plots are to be made." (LCM Manual p. 118)
ControlWrite.WDLINE =  {'WDLINE(6) = 0.0'};                        % Set the fine grid lines to thickness = 0. 

% Number of Spline Knots used
ControlWrite.DKNTMN = 'DKNTMN = 10.0';  			 % This forces a very flat baseline  %GH 2HG exp: try one small and one high

% Analysis Window
ControlWrite.PPMST = 'PPMST = 4.2';                              % Fit data in chemical shift region [PPMEND, PPMST], PPMST > PPMEND
ControlWrite.PPMEND = 'PPMEND = 0.2';
ControlWrite.PPMGAP = {'PPMGAP(1,1)=1.75', 'PPMGAP(2,1)=1.15 '};

% Zero Order Phase
ControlWrite.DEGZER = 'DEGZER = 0';                              % zero order phase prior knowledge, set to zero for no prior knowledge
ControlWrite.SDDEGZ = 'SDDEGZ = 20';                             % standard deviation of DEGZER, set to 999 for no prior knowledge.

% First order phase
ControlWrite.DEGPPM = 'DEGPPM = 0';                              % 1st order phase prior knowledge, set to zero for no prior knowledge
ControlWrite.SDDEGP = 'SDDEGP = 20';                             % standard deviation of DEGPPM; note that LCM varies the phase a lot: E.g. for sddegp=1 a total 1.order_phase > 15 is not rare!
%%%%% FIRST ORDER PHASE:
%%%%% if 2 metabolites have circle freq w1, w2 and a phase 0. order of ph0 and they are measured with acquisition delay t
%%%%% then their phases in rad are:
%%%%% phw1(t) = ph0 + w1*t
%%%%% phw2(t) = ph0 + w2*t
%%%%% so their phase difference is delta_ph[rad] = (w2-w1)*t
%%%%% the acq delay is then t = delta_ph[rad]/(w2-w1) in seconds, converted to degree/ppm it is 
%%%%% t [s] = delta_ph[rad]/(w2-w1) = delta_ph[rad]/((f2-f1)*2pi) = delta_ph[rad]/(ppm_difference*297.223*10^6*2pi);  (297.223 @ 7T !)
%%%%% so assuming ppm_difference = 2ppm and a delta_ph[rad]=0.05236rad  equal to 3deg --> 1stOrderPhase = 0.0262 rad/ppm equal to 1.5 deg/ppm
%%%%% this leads to an acq delay of: t ~ 0.014 ms; 
%%%%% so such an deviation (because of wrong basis file or technical inaccuracy) can be compensated



% Concentration ratio priors
ControlWrite.NRATIO = 'NRATIO = 10';                          
ControlWrite.CHRATO = { 'CHRATO(1) = ''MM_1_2/MM_0_9 = 0.12 +- 0.10''', 'CHRATO(2) = ''MM_1_4/MM_0_9 = 0.235 +- 0.15''', 'CHRATO(3) = ''MM_1_6/MM_0_9 = 0.24 +- 0.3''', 'CHRATO(4) = ''MM_2_0/MM_0_9 = 0.74 +- 0.25''', 'CHRATO(5) = ''MM_2_2/MM_0_9 = 0.6 +- 0.3''', 'CHRATO(6) = ''MM_2_9/MM_0_9 = 0.30 +- 0.2''', 'CHRATO(7) = ''MM_3_2/MM_0_9 = 0.32 +- 0.2''', 'CHRATO(8) = ''MM_3_7/MM_0_9 = 1.01 +- 0.3''', 'CHRATO(9) = ''MM_2_9/MM_3_2 = 0.99 +- 0.1''', 'CHRATO(10) = ''MM_2_2/MM_2_0 = 0.81 +- 0.1''' }

% Basis Set Parameters
ControlWrite.NSIMUL = 'NSIMUL = 0';                              % Simulate Lipid signals to account for lipid contamination in 1.3ppm (provided by Marjanska,Jan2016)

%ControlWrite.CHSIMU = { 'CHSIMU(1) = ''Lip13a @ 1.28 +- .01 FWHM= .15 < .2 +- .035 AMP= 2.''', 'CHSIMU(2) = ''Lip13b @ 1.28 +- .01 FWHM= .089 < .09 +- .035 AMP= 2.''', 'CHSIMU(3) = ''Lip13c @ 1.30 +- .01 FWHM= .089 < .09 +- .035 AMP= 2.''', 'CHSIMU(4) = ''Lip13d @ 1.26 +- .01 FWHM= .089 < .09 +- .035 AMP= 2.''', 'CHSIMU(5) = ''Lip09 @ .89 +- .02 FWHM= .14 < .19 +- .035 AMP= 3.''', 'CHSIMU(6) = ''Lip12n @1.12 +-.05 FWHM= .15 < .2 +- .35 AMP= -1.''', 'CHSIMU(7) = ''Lip14 @1.44 +-.03 FWHM= .17 < .2 +- .2 AMP= 2.''' }


ControlWrite.NOMIT =  'NOMIT = 8';                               % Number of Metabolites within the Basis Set that should be omitted from the analysis
ControlWrite.CHOMIT =  {'CHOMIT(1) = ''Cho''','CHOMIT(2) = ''Act''','CHOMIT(3) = ''mm3''','CHOMIT(4) = ''mm4''','CHOMIT(5) = ''Glc_B''', 'CHOMIT(6) = ''Lip_c''', 'CHOMIT(7) = ''Glc''', 'CHOMIT(8) = ''TwoHG'''}; % Names of omitted metabolites



% ControlWrite.NUSE =  'NUSE = 2'; 
% ControlWrite.CHUSE1 = {'CHUSE(1)=''Act''', 'CHUSE(2)=''Lac'''};% Only Use the following metabolites in the Preliminary Analysis.

%ControlWrite.NKEEP =  'NKEEP = 1';                               % Keep These Metabolites in the Analysis, even if only minor peaks are in the Analysis Window. (LCM Manual ver. 6.3-1, p. 118)
%ControlWrite.CHKEEP =  {'CHKEEP(1) = ''Lac'''}; 


% This is a little bit of a hack, since it has in principle nothing to do in this file
% But it works!
CPU_cores = 16;



% You can also set other parameters by assigning ControlWrite.Others[n], where [n] is a natural number.
% These things are written literally to the LCModel Control File, so they always have to be strings.
% Examples:
% ControlWrite.Others1 = 'NAMREL = ''NAA''';		% Relative concentrations to NAA instead of Cr.
% ControlWrite.Others2 = 'ROOMT = T';			% Measurement was performed at Room Temperature.


% Controls for creating different files
ControlWrite.LTABLE =  'LTABLE = 7';         % Create a .table file
ControlWrite.LCSV =  'LCSV = 0';             % Don't create a .CSV file 
ControlWrite.LCOORD =  'LCOORD = 9';         % Create a Coord file
