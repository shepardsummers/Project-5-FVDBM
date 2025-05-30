clear, clc;

order = 1;

ksi = [0  1  0 -1  0  1 -1 -1  1; ...
       0  0  1  0 -1  1  1 -1 -1
      ];

L = [0.8, 1.2, 0.8, 1.2];
d = [1.2, 0.8, 1.2, 0.8];
n = [1  0 -1  0; ...
     0  1  0 -1
    ];

dt = 0.1;

f = [1.67; 0.43; 0.42; 0.4; 0.42; 0.11; 0.1; 0.1; 0.11];
f_n = [1.63 1.67 1.55 1.66; ...
       0.61 0.42 0.77 0.50; ...
       0.41 0.42 0.39 0.42; ...
       0.27 0.42 0.20 0.35; ...
       0.41 0.42 0.39 0.42; ...
       0.15 0.10 0.19 0.12; ...
       0.07 0.11 0.05 0.09; ...
       0.07 0.10 0.05 0.08; ...
       0.16 0.11 0.20 0.13
      ];

tic;
guh = zeros(9,4);
guh2 = guh;
erm = zeros(9,1);
upp = zeros(9,1);
dwn = zeros(9,1);

for a = 1:4
    erm = ksi'*n(:,a);
    switch order
        case 1
            for i = 1:9
                if erm(i) > 0
                    erm(i) = f(i);
                else
                    erm(i) = f_n(i,a);
                end
            end
        case 2
            for i = 1:9
                if erm(i) > 0
                    upp(i) = f(i);
                    dwn(i) = f_n(i,a);
                else
                    upp(i) = f_n(i,a);
                    dwn(i) = f(i);
                end
            end
            erm = upp + (dwn - upp).*(1/2 - abs(erm)*dt/(2*d(a)));
    end
    guh(:,a) = erm.*ksi'*n(:,a)*L(a);
end

flux1 = sum(guh, 2);
time1 = toc;

tic;
    kdn = ksi'*n;

    f_t = (max(kdn, 0)~=0).*f + (min(kdn, 0)~=0).*f_n;

    flux2 = sum(f_t.*(ksi'*n).*L, 2);
time2 = toc;

tic;
    flux3 = flux_int(n,L,f,f_n,ksi);
time3 = toc;