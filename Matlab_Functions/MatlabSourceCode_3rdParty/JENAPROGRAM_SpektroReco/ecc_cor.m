function d1 = ecc_cor(d,d_ref)

if iscell(d),       d_temp = d{1}; clear d; d = d_temp; clear d_temp;               end
if iscell(d_ref),   d_temp = d_ref{1}; clear d_ref; d_ref = d_temp; clear d_temp;   end

d       = d(:);
d_ref   = d_ref(:);

cor     = exp(1i*angle(d_ref));

d1      = d./cor;