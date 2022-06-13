function [d1_comb, w, SNR] = chan_combine_3(d)

% SVD based chanel combination method
% d         = channel data matrix
% d1        = weighted chanel data matrix

[N_fid,dim_fid]     = max(size(d));     % number of fid data points     [Y,I] = max(x): Y ... maximum of x; I ... Index of Y within x
[N_chan,dim_chan]   = min(size(d));     % number of chanels

[U,dums,V]             = svd(d);           % SVD
clear dums

if dim_chan < dim_fid
    w       = abs(U(:,1))./sum(abs(U(:,1)));        % global weigth factors
    for k = 1:N_chan,   d1(k,1:N_fid) = w(k).*d(k,:); end
    d1_comb = sum(d1,dim_chan);
    SNR     = abs(d1_comb(1,1))/2*std(abs(d1_comb(1,N_fid-100:N_fid)));
else
    w       = (abs(V(:,1))./abs(V(:,1)))';
    SNR     = abs(d1_comb(1,1))/sqrt(abs(var(d1_comb(N_fid-100:N_fid,1))));
end