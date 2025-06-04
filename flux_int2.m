function [flux] = flux_int2(n,L,f,f_n,ksi,d_t,d)
    % Computes interior flux

    kdn = ksi'*n;

    f_up = (max(kdn, 0)~=0).*f + (min(kdn, 0)~=0).*f_n;
    f_down = (min(kdn, 0)~=0).*f + (max(kdn, 0)~=0).*f_n;
    
    erm = (f_up + (f_down - f_up).*(0.5 - abs(kdn)*d_t./(2*d)));
    
    flux = sum(erm.*kdn.*L, 2);

end

