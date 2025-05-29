%% Definition of Parameters
clear, clc;

% Domain Related
N_nd_x = 100; % num of x nodes
N_nd_y = N_x; % num of y nodes
L = 1; % total length
H = 1; % total height
N_cv_x = N_x - 1; % number of x cells
N_cv_y = N_y - 1; % number of y cells
d_x = L/N_cv_x; % dist between x nodes
d_y = H/N_cv_y; % dist between y nodes
A=dx*dy;

n = [1  0 -1  0; ...
     0  1  0 -1]; % normals (CCW starting from East)

% DBM Related
% Lattice velocity
ksi = [0 1 0 -1 0 1 -1 -1 1; ...
       0 0 1 0 -1 1 1 -1 -1 ];

w = [4/9 1/9 1/9 1/9 1/9 1/36 1/36 1/36 1/36]; % weights for D2Q9

c_s = 1/sqrt(3); % speed of sound (D2Q9)

Re = 100; % Reanolds number
Tau = 0.4; % relaxation time
vis = Tau*c_s^2; % kinematic viscosity
U_top = Re*vis/L; % Top velocity
min_error = 0.000001; % error when sim ends

%% Initialization
Rho_init=2;
Rho_nd = ones(1, N_nd_y, N_nd_x)*Rho_init; % Density nodes
Rho_cv = ones(1, N_cv_y, N_cv_x)*Rho_init; % Density cells
U_nd = zeros(2, N_nd_y, N_nd_x); % Velocity nodes
U_cv = zeros(2, N_cv_y, N_cv_x); % Velocity cells

%f_nd = zeros(9, N_nd_y, N_nd_x); % PDF for all 9 directions at all node locations
%f_cv = zeros(9, N_cv_y, N_cv_x); % PDF for all 9 directions at all cell locations

f_nd = eqm_d2q9(Rho_nd, U_nd, ksi, w); % PDF for all 9 directions at all node locations
f_cv = eqm_d2q9(Rho_cv, U_cv, ksi, w); % PDF for all 9 directions at all cell locations

f_nd_new = f_nd; % Update variable nodes
f_cv_new = f_nd; % Update variable cells
f_nd_eq = f_nd_new; % Equilibrium nodes
f_cv_eq = f_cv_new; % Equilibrium cells

%% Timer
timer = 10000;

%% Solving
tic
res_list = zeros(1, max_timer);

% % Interior nodes
% int_x = 2:(N_x-1);
% int_y = 2:(N_x-1);

for t = 1:timer
    % Compute boundary nodal PDF

    % Flux calculation
    
    % Collision
    Rho_old = Rho_cv;
    % Rho, U calculation
    [Rho_cv, U_cv] = rhoNu(f_cv_new, ksi);
    % f_eq calculation
    f_cv_eq = eqm_d2q9(Rho_cv, U_cv, ksi, w);

    % BGK Collision and Update
    f = f_new - (f_new-f_eq)/Tau; % needs to be changed

    [guh, res_list(t)] = res(Rho_old, Rho, min_error);
    
    %addpoints(r, t, max_error)
    %drawnow

    %progress(timer, t);
    fprintf("Itt: %i     ||     Res: %.4e\n", t, res_list(t))
end
total_time = toc;
%% Post-Processing / Visualization
clc;
load Ghia_Re100.mat
figure
quiver(flipud(squeeze(U(1,:,:))),flipud(squeeze(U(2,:,:))),10)
axis equal tight

figure
contourf(flipud(squeeze(Rho)),30)
axis equal tight

mid = N_x/2;
Vertical_Sample = U(1, :, mid)/U_top;
Horizontal_Sample = U(2, mid, :)/U_top;

%u2_Ghia = [1 0.48223 0.46120 0.45992 0.46036 0.33556 0.20087 0.08183 -0.03039 -0.07404 -0.22855 -0.33050 -0.40435 -0.43643 -0.42901 -0.41165 0.00000 ];
%v2_Ghia = [0.00000 -0.49774 -0.55069 -0.55408 -0.52876 -0.41442 -0.36214 -0.30018 0.00945 0.27280 0.28066 0.35368 0.42951 0.43648 0.43329 0.42447 0.00000];

figure
plot(squeeze(Vertical_Sample), flip((1:L)/L), "black", flip(u_Ghia), flip(y_Ghia))
title("Vertical Sample (U)")
xlabel("u")
ylabel("y")
figure
plot((1:L)/L, squeeze(Horizontal_Sample), "black", flip(x_Ghia), flip(v_Ghia));
title("Horizontal Sample (V)")
xlabel("x")
ylabel("v")

figure
u = flip(squeeze(U(1, :, :)));
v = squeeze(U(2, :, :));
[startX, startY] = meshgrid(1:50:N_x, 1:50:N_y);
verts = stream2(1:N_x,1:N_y,u,v,startX,startY);
streamline(verts)