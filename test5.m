clear, clc;

order = 2;

ksi = [0  1  0 -1  0  1 -1 -1  1; ...
       0  0  1  0 -1  1  1 -1 -1
      ];

L = [0.8, 1.2, 0.8, 1.2];
d = [1.2, 0.8, 1.2, 0.8];
n = [1  0 -1  0; ...
     0  1  0 -1
    ];

w = [4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36];

dt = 0.1;

f = [1.63 1.67 1.66; ...
     0.61 0.42 0.50; ...
     0.41 0.42 0.42; ...
     0.27 0.42 0.35; ...
     0.41 0.42 0.42; ...
     0.15 0.10 0.12; ...
     0.07 0.11 0.09; ...
     0.07 0.10 0.08; ...
     0.16 0.11 0.13
    ];

[rho, U] = rhoNu(f, ksi);
f_eq = eqm_d2q9(rho, U, ksi, w);

n_cell = 3;
rho_n = zeros(1,n_cell-1);
U_n = zeros(2,2);

f_n_neq = zeros(9,2);

for i = 1:n_cell-1
    rho_n(i) = (rho(i) + rho(i+1))/2;
    f_n_neq(:,i) = (f(:,i) - f_eq(:,i) + f(:,i+1) - f_eq(:,i+1))/2;
end

f_n_eq = eqm_d2q9(rho_n, U_n, ksi, w);
%f_n_neq = eqm_d2q9(rho_n, U_n, ksi, w);

f_n = f_n_eq + f_n_neq;

f_g = (sum(f_n,2)) - f(:,2);

% guh = zeros(9,4);
% guh2 = guh;
% erm = zeros(9,1);
% upp = zeros(9,1);
% dwn = zeros(9,1);

erm = ksi'*n(:,4);
for i = 1:9
    if erm(i) > 0
        erm(i) = f(i);
    else
        erm(i) = f_g(i);
    end
end

flux_b = erm.*ksi'*n(:,4)*L(4);