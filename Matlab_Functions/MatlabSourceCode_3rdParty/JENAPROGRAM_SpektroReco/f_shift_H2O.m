function w_dif = f_shift_H2O(d,threshold,h)

% INPUT:
%       d               = FID-data array
%       threshold       = limit for f-shift value (in Hz) for exclusion of artifact contamined data
%       h               = rda-header

% OUTPUT:
%       w_dif           = vector with f-shift, calculated for each FID
%       index_include   = included FID-data indexes


%**************************************************************************
% get the rda-header
if      nargin == 2
    [~,h]       = rda_reader_2(0);
elseif  nargin == 1
    threshold   = 5;
    [~,h]       = rda_reader_2(0);
end

% *************************************************************************
% FID-data array dimensions

[N_fid,dim_fid]     = min(size(d));     % number of FID's
[N ,dim_samp]       = max(size(d));     % number of sampled FID data points

d                   = permute(d,[dim_fid,dim_samp]);
N_store             = N;

% *************************************************************************
% header information (actual file)
for k = 1:max(size(h))
    if strcmp(h{k,1},'MRFrequency'),    w0      = str2double(h{k,2});               end
    if strcmp(h{k,1},'DwellTime'),      dt      = 0.000001*str2double(h{k,2});      end
    if strcmp(h{k,1},'VectorSize'),     N_temp  = str2double(h{k,2});               end        
end

if N_temp == N/2,   dt = dt/2; end

% ****************************t*********************************************
% calculate the time-, frequency- and ppm-scales

t           = (0:N-1).*dt;                                      % time scale
t           = reshape(t,[1 N]);
dw_store    = (1/dt)/N;
w_store     = (-(N-0.5)/2:1:(N-0.5)/2).*dw_store;               % why -0.5 ?

% *************************************************************************
% apodizate data (exp. damping with 20Hz)
d_apod              = zeros(N_fid,N);
for k = 1:N_fid,    d_apod(k,:) = exp(-20.*t).*d(k,:); end

% d_apod              = d;

% *************************************************************************
% zero-filling 

f_zero_filling      = 1;

% N_stored            = N;
zero_vector         = complex(zeros(1,N*(2^f_zero_filling)-N),zeros(1,N*(2^f_zero_filling)-N));
N                   = N*(2^f_zero_filling);
d_zer_fil           = zeros(N_fid,N);
for k = 1:N_fid,    d_zer_fil(k,:) = cat(2,d_apod(k,:),zero_vector);   end

% *************************************************************************
% calculate the frequency shifts relative to the first FID

w_dif               = zeros(1,N_fid);
D                   = zeros(N_fid,N);
for k = 1:N_fid,    D(k,1:N) = fftshift(fft(d_zer_fil(k,:))); end

% define the H2O position (within the intervall 'threshold/2:threshold/2', defined around the nominal water peak frequency (@ 0Hz))
dw                  = (1/dt)/N;                         % dw after the zero-filling
w                   = (-(N-0.5)/2:1:(N-0.5)/2).*dw;     % frequency scale (Hz)

i_start_H2O         = fix(median(find(w>-8 & w<8)) - 0.5*threshold);
i_fin_H2O           = fix(median(find(w>-8 & w<8)) + 0.5*threshold);

w_H2O_ref           = 0;

% extract the actual H2O-frequencies and calculate the
% f-differences relative to the nominal H2O-frequencies
for k = 1:N_fid
    temp            = abs(D(k,:));    
    H2O             = max(temp(1,i_start_H2O:i_fin_H2O));
    w_dif_temp      = w_H2O_ref - w(temp == H2O); clear H2O
    
    if      abs((mod(w_dif_temp,dw_store)/dw_store)) > 0.5,     w_dif(k) = round(w_dif_temp/dw_store)*dw_store;
    else                                                        w_dif(k) = fix(w_dif_temp/dw_store)*dw_store;
    end
    clear temp
    
    temp_start      = fix(median(find(w_store>-8 & w_store<8)) - 0.5*threshold);
    temp_fin        = fix(median(find(w_store>-8 & w_store<8)) + 0.5*threshold);
    
    temp            = abs(fftshift(fft(d_apod(k,:).*exp(1i*2*pi*w_dif(k)*t))));         % f-corrected
    H2O             = max(temp(1,temp_start:temp_fin));
    temp_dif1       = abs(w_H2O_ref - w_store(temp == H2O));
    
    temp            = abs(fftshift(fft(d_apod(k,:))));                                  % f-uncorrected
    H2O             = max(temp(1,temp_start:temp_fin));
    temp_dif2       = abs(w_H2O_ref - w_store(temp == H2O));
    
    if temp_dif1 > temp_dif2, w_dif(k) = 0; end
    
    clear temp temp_start temp_fin temp_dif w_dif_temp
    
end