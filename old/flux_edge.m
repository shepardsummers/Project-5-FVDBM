function [flux] = flux_edge(n,L,f,f_n,ksi,f_e,x)
    % Computes edge flux
    

    kdn = ksi'*n;

    f_t = (max(kdn, 0)~=0).*f + (min(kdn, 0)~=0).*f_n;

    f_t(:,x) = f_e;

    flux = sum(f_t.*(ksi'*n).*L, 2);
end