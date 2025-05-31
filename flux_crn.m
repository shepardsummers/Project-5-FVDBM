function [flux] = flux_crn(n,L,f,f_n,ksi,f_e,x,f_e2,x2)
    % Computes corner flux
    

    kdn = ksi'*n;

    f_t = (max(kdn, 0)~=0).*f + (min(kdn, 0)~=0).*f_n;

    f_t(:,x) = f_e;
    f_t(:,x2) = f_e2;

    flux = sum(f_t.*(ksi'*n).*L, 2);
end