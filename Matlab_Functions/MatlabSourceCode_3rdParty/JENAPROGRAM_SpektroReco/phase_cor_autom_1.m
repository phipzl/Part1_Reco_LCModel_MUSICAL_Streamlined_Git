function [d1,phi0_opt] = phase_cor_autom_1(d)

d           = squeeze(d);

step        = 10;
phi0_min    = -180;
phi0_max    = 180;

while step > 0.1
    phi0    = phi0_min:step:phi0_max;
    N_phi0  = length(phi0);
    
    for k = 1:N_phi0
        temp        = d.*exp(1i*phi0(k)*(pi/180));
        temp1       = fftshift(fft(temp));
    
        sum_re(k)   = sum(real(temp1));
        sum_im(k)   = sum(imag(temp1));        
    end
    
    delta           = sum_re - abs(sum_im);
    [ma,k_opt]    = max(delta); clear ma
    
    if          k_opt == 1
        phi0_min        = phi0(k_opt);
        phi0_max        = phi0(k_opt+1);
    elseif      k_opt == length(phi0)
        phi0_min        = phi0(k_opt-1);
        phi0_max        = phi0(k_opt);
    else
        phi0_min        = phi0(k_opt-1);
        phi0_max        = phi0(k_opt+1);
    end
    
    step            = step/2;
end

phi0_opt        = phi0(k_opt);
d1              = d.*exp(1i*phi0(k_opt)*(pi/180));
