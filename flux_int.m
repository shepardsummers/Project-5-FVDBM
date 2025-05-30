function [flux] = flux_int(n,L,f,f_n,ksi)
    % Computes interior flux

    kdn = ksi'*n;

    f_t = (max(kdn, 0)~=0).*f + (min(kdn, 0)~=0).*f_n;

    flux = sum(f_t.*(ksi'*n).*L, 2);
end

