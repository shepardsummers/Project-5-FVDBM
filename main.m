%% Definition of Parameters
clear, clc;

% Domain Related
N_nd_x = 100; % num of x nodes
N_nd_y = N_nd_x; % num of y nodes
L = 100; % total length
H = 100; % total height
N_cv_x = N_nd_x - 1; % number of x cells
N_cv_y = N_nd_y - 1; % number of y cells
d_x = L/N_cv_x; % dist between x nodes
d_y = H/N_cv_y; % dist between y nodes
A=d_x*d_y;
d_t = 0.1;
Len = [d_y,d_x,d_y,d_x];

f_ept = zeros(9,1);

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
f_cv_new = f_cv; % Update variable cells
f_nd_eq = f_nd_new; % Equilibrium nodes
f_cv_eq = f_cv_new; % Equilibrium cells
f_nd_neq = f_nd_new; % Non equilibrium nodes
f_cv_neq = f_cv_new; % Non equilibrium cells

flux = zeros(9,N_cv_y,N_cv_x);
%% Timer
timer = 3000;

%% Solving
tic
res_list = zeros(1,timer);

% % Interior nodes
% int_x = 2:(N_x-1);
% int_y = 2:(N_x-1);

for t = 1:timer
    % Compute boundary nodal PDF
    for j = 1:N_nd_y
        for i = 1:N_nd_x
            if j == 1 % Top boundary
                if i == 1 % Top left corner
                    % f_eq
                        Rho_nd(1,j,i) = Rho_cv(1,1,1);
                        f_nd_eq(:,j,i) = eqm_d2q9(Rho_nd(1,j,i),[U_top;0],ksi,w);
                        % Try (U_lid + 0) / 2 for x velocity
                    % f_neq
                        f_cv_1 = f_cv(:,1,1);
                        f_cv_eq_1 = f_cv_eq(:,1,1);
                        
                        f_nd_neq(:,j,i) = (f_cv_1 - f_cv_eq_1);
                    % f
                        f_nd(:,j,i) = f_nd_eq(:,j,i) + f_nd_neq(:,j,i);
                elseif i == N_nd_x % Top right corner
                    % f_eq
                        Rho_nd(1,j,i) = Rho_cv(1,1,end);
                        f_nd_eq(:,j,i) = eqm_d2q9(Rho_nd(1,j,i),[U_top;0],ksi,w);
                    % f_neq
                        f_cv_1 = f_cv(:,1,end);
                        f_cv_eq_1 = f_cv_eq(:,1,end);
                        
                        f_nd_neq(:,j,i) = (f_cv_1 - f_cv_eq_1);
                    % f
                        f_nd(:,j,i) = f_nd_eq(:,j,i) + f_nd_neq(:,j,i);
                else % Top nodes
                    % f_eq
                        Rho_nd(1,j,i) = (Rho_cv(1,1,i) + Rho_cv(1,1,i-1))/2;
                        f_nd_eq(:,j,i) = eqm_d2q9(Rho_nd(1,j,i),[U_top;0],ksi,w);
                    % f_neq
                        f_cv_1 = f_cv(:,1,i);
                        f_cv_2 = f_cv(:,1,i-1);
                        f_cv_eq_1 = f_cv_eq(:,1,i);
                        f_cv_eq_2 = f_cv_eq(:,1,i-1);
                        
                        f_nd_neq(:,j,i) = ((f_cv_1 - f_cv_eq_1) + (f_cv_2 - f_cv_eq_2))/2;
                    % f
                        f_nd(:,j,i) = f_nd_eq(:,j,i) + f_nd_neq(:,j,i);
                end
            elseif j == N_nd_y % Bottom boundary
                if i == 1 % Bottom left corner
                    % f_eq
                        Rho_nd(1,j,i) = Rho_cv(1,end,1);
                        f_nd_eq(:,j,i) = eqm_d2q9(Rho_nd(1,j,i),[0;0],ksi,w);
                    % f_neq
                        f_cv_1 = f_cv(:,end,1);
                        f_cv_eq_1 = f_cv_eq(:,end,1);
                        
                        f_nd_neq(:,j,i) = (f_cv_1 - f_cv_eq_1);
                    % f
                        f_nd(:,j,i) = f_nd_eq(:,j,i) + f_nd_neq(:,j,i);
                elseif i == N_nd_x % Bottom right corner
                    % f_eq
                        Rho_nd(1,j,i) = Rho_cv(1,end,end);
                        f_nd_eq(:,j,i) = eqm_d2q9(Rho_nd(1,j,i),[0;0],ksi,w);
                    % f_neq
                        f_cv_1 = f_cv(:,end,end);
                        f_cv_eq_1 = f_cv_eq(:,end,end);
                        
                        f_nd_neq(:,j,i) = (f_cv_1 - f_cv_eq_1);
                    % f
                        f_nd(:,j,i) = f_nd_eq(:,j,i) + f_nd_neq(:,j,i);
                else % Bottom nodes
                    % f_eq
                        Rho_nd(1,j,i) = (Rho_cv(1,end,i) + Rho_cv(1,end,i-1))/2;
                        f_nd_eq(:,j,i) = eqm_d2q9(Rho_nd(1,j,i),[0;0],ksi,w);
                    % f_neq
                        f_cv_1 = f_cv(:,end,i);
                        f_cv_2 = f_cv(:,end,i-1);
                        f_cv_eq_1 = f_cv_eq(:,end,i);
                        f_cv_eq_2 = f_cv_eq(:,end,i-1);
                        
                        f_nd_neq(:,j,i) = ((f_cv_1 - f_cv_eq_1) + (f_cv_2 - f_cv_eq_2))/2;
                    % f
                        f_nd(:,j,i) = f_nd_eq(:,j,i) + f_nd_neq(:,j,i);
                end
            elseif i == 0 % Left boundary
                % f_eq
                    Rho_nd(1,j,i) = (Rho_cv(1,j,1) + Rho_cv(1,j-1,1))/2;
                    f_nd_eq(:,j,i) = eqm_d2q9(Rho_nd(1,j,i),[0;0],ksi,w);
                % f_neq
                    f_cv_1 = f_cv(:,j,1);
                    f_cv_2 = f_cv(:,j-1,1);
                    f_cv_eq_1 = f_cv_eq(:,j,1);
                    f_cv_eq_2 = f_cv_eq(:,j-1,1);
                    
                    f_nd_neq(:,j,i) = ((f_cv_1 - f_cv_eq_1) + (f_cv_2 - f_cv_eq_2))/2;
                % f
                    f_nd(:,j,i) = f_nd_eq(:,j,i) + f_nd_neq(:,j,i);
            elseif i == N_nd_x % Right Boundary
                % f_eq
                    Rho_nd(1,j,i) = (Rho_cv(1,j,end) + Rho_cv(1,j-1,end))/2;
                    f_nd_eq(:,j,i) = eqm_d2q9(Rho_nd(1,j,i),[0;0],ksi,w);
                % f_neq
                    f_cv_1 = f_cv(:,j,end);
                    f_cv_2 = f_cv(:,j-1,end);
                    f_cv_eq_1 = f_cv_eq(:,j,end);
                    f_cv_eq_2 = f_cv_eq(:,j-1,end);
                    
                    f_nd_neq(:,j,i) = ((f_cv_1 - f_cv_eq_1) + (f_cv_2 - f_cv_eq_2))/2;
                % f
                    f_nd(:,j,i) = f_nd_eq(:,j,i) + f_nd_neq(:,j,i);
            else % Interior
                % Do nothing
            end
        end
    end
    % Flux calculation
    for j = 1:N_cv_y
        for i = 1:N_cv_x
            if j == 1 % Top
                if i == 1 % Top left
                    % f_e = f_nd(:,1,1) - f_cv(:,j,i);
                    % f_e2 = f_nd(:,1,1) - f_cv(:,j,i);
                    % f_n = [f_cv(:,j,i+1), f_ept, f_ept, f_cv(:,j+1,i)];
                    % flux(:,j,i) = flux_crn(n,Len,f_cv(:,j,i),f_n,ksi,f_e,2,f_e2,3);

                    f_g1 = (f_nd(:,1,1) + f_nd(:,1+1,1)) - f_cv(:,j,i); % left
                    f_g2 = (f_nd(:,1,1) + f_nd(:,1,1+1)) - f_cv(:,j,i); % top
                    f_n = [f_cv(:,j,i+1), f_g2, f_g1, f_cv(:,j+1,i)];
                    flux(:,j,i) = flux_int(n,Len,f_cv(:,j,i),f_n,ksi);
                elseif i == N_cv_x % Top right
                    % f_e = f_nd(:,1,end) - f_cv(:,j,i);
                    % f_e2 = f_nd(:,1,end) - f_cv(:,j,i);
                    % f_n = [f_ept, f_ept, f_cv(:,j,i-1), f_cv(:,j+1,i)];
                    % flux(:,j,i) = flux_crn(n,Len,f_cv(:,j,i),f_n,ksi,f_e,2,f_e2,1);

                    f_g1 = (f_nd(:,1,end) + f_nd(:,1+1,end)) - f_cv(:,j,i); % right
                    f_g2 = (f_nd(:,1,end) + f_nd(:,1,end-1)) - f_cv(:,j,i); % top
                    f_n = [f_g1, f_g2, f_cv(:,j,i-1), f_cv(:,j+1,i)];
                    flux(:,j,i) = flux_int(n,Len,f_cv(:,j,i),f_n,ksi);
                else % Other top
                    % f_e = f_nd(:,1,i-1) - f_cv(:,j,i);
                    % f_n = [f_cv(:,j,i+1), f_ept, f_cv(:,j,i-1), f_cv(:,j+1,i)];
                    % flux(:,j,i) = flux_edge(n,L,f_cv(:,j,i),f_n,ksi,f_e,2);

                    f_g = (f_nd(:,1,i) + f_nd(:,1,i+1)) - f_cv(:,j,i);
                    f_n = [f_cv(:,j,i+1), f_g, f_cv(:,j,i-1), f_cv(:,j+1,i)];
                    flux(:,j,i) = flux_int(n,Len,f_cv(:,j,i),f_n,ksi);
                end
            elseif j == N_cv_y % Bottom
                if i == 1 % Bottom left
                    % f_e = f_nd(:,end,1) - f_cv(:,j,i);
                    % f_e2 = f_nd(:,end,1) - f_cv(:,j,i);
                    % f_n = [f_cv(:,j,i+1), f_cv(:,j-1,i), f_ept, f_ept];
                    % flux(:,j,i) = flux_crn(n,Len,f_cv(:,j,i),f_n,ksi,f_e,4,f_e2,3);

                    f_g1 = (f_nd(:,end,1) + f_nd(:,end-1,1)) - f_cv(:,j,i); % left
                    f_g2 = (f_nd(:,end,1) + f_nd(:,end,1+1)) - f_cv(:,j,i); % bottom
                    f_n = [f_cv(:,j,i+1), f_cv(:,j-1,i), f_g1, f_g2];
                    flux(:,j,i) = flux_int(n,Len,f_cv(:,j,i),f_n,ksi);
                elseif i == N_cv_x % Bottom right
                    % f_e = f_nd(:,end,end) - f_cv(:,j,i);
                    % f_e2 = f_nd(:,end,end) - f_cv(:,j,i);
                    % f_n = [f_ept, f_cv(:,j-1,i), f_cv(:,j,i-1), f_ept];
                    % flux(:,j,i) = flux_crn(n,Len,f_cv(:,j,i),f_n,ksi,f_e,4,f_e2,1);

                    f_g1 = (f_nd(:,end,end) + f_nd(:,end-1,end)) - f_cv(:,j,i); % right
                    f_g2 = (f_nd(:,end,end) + f_nd(:,end,end-1)) - f_cv(:,j,i); % bottom
                    f_n = [f_g1, f_cv(:,j-1,i), f_cv(:,j,i-1), f_g2];
                    flux(:,j,i) = flux_int(n,Len,f_cv(:,j,i),f_n,ksi);
                else % Other bottom
                    % f_e = f_nd(:,end,i-1) - f_cv(:,j,i);
                    % f_n = [f_cv(:,j,i+1), f_cv(:,j-1,i), f_cv(:,j,i-1), f_ept];
                    % flux(:,j,i) = flux_edge(n,L,f_cv(:,j,i),f_n,ksi,f_e,4);

                    f_g = (f_nd(:,end,i) + f_nd(:,end,i+1)) - f_cv(:,j,i);
                    f_n = [f_cv(:,j,i+1), f_cv(:,j-1,i), f_cv(:,j,i-1), f_g];
                    flux(:,j,i) = flux_int(n,Len,f_cv(:,j,i),f_n,ksi);
                end
            elseif i == 1 % Left
                % f_e = f_nd(:,j-1,1) - f_cv(:,j,i);
                % f_n = [f_cv(:,j,i+1), f_cv(:,j-1,i), f_ept, f_cv(:,j+1,i)];
                % flux(:,j,i) = flux_edge(n,L,f_cv(:,j,i),f_n,ksi,f_e,3);

                f_g = (f_nd(:,j,1) + f_nd(:,j+1,1)) - f_cv(:,j,i);
                f_n = [f_cv(:,j,i+1), f_cv(:,j-1,i), f_g, f_cv(:,j+1,i)];
                flux(:,j,i) = flux_int(n,Len,f_cv(:,j,i),f_n,ksi);
            elseif i == N_cv_x % Right 
                % f_e = f_nd(:,j-1,end) - f_cv(:,j,i);
                % f_n = [f_ept, f_cv(:,j-1,i), f_cv(:,j,i-1), f_cv(:,j+1,i)];
                % flux(:,j,i) = flux_edge(n,L,f_cv(:,j,i),f_n,ksi,f_e,1);

                f_g = (f_nd(:,j,end) + f_nd(:,j+1,end)) - f_cv(:,j,i);
                f_n = [f_g, f_cv(:,j-1,i), f_cv(:,j,i-1), f_cv(:,j+1,i)];
                flux(:,j,i) = flux_int(n,Len,f_cv(:,j,i),f_n,ksi);
                
            else % Interior
                f_n = [f_cv(:,j,i+1), f_cv(:,j-1,i), f_cv(:,j,i-1), f_cv(:,j+1,i)];
                flux(:,j,i) = flux_int(n,Len,f_cv(:,j,i),f_n,ksi);
            end
        end
    end
    % Collision
    Rho_old = Rho_cv;
    % Rho, U calculation
    [Rho_cv, U_cv] = rhoNu(f_cv, ksi);
    % f_eq calculation
    %f_cv_eq_old = f_cv_eq;
    f_cv_eq = eqm_d2q9(Rho_cv, U_cv, ksi, w);

    % 
    %f_cv_eq = 2*f_cv_eq - f_cv_eq_old;
    %f_cv = Tau/(Tau + d_t)*f_cv + d_t/(Tau + d_t)*f_cv_eq - (Tau * d_t)/(Tau + d_t)*flux/A;
    f_cv = ((Tau - d_t)/Tau)*f_cv + (d_t/Tau)*f_cv_eq - d_t*flux/A;

    [guh, res_list(t)] = res(Rho_old, Rho_cv, min_error);
    
    %addpoints(r, t, max_error)
    %drawnow

    %progress(timer, t);
    fprintf("Itt: %i     ||     Res: %.4e\n", t, res_list(t))
end
total_time = toc;
%% Post-Processing / Visualization
load Ghia_Re100.mat
figure
quiver(flipud(squeeze(U_cv(1,:,:))),flipud(squeeze(U_cv(2,:,:))),10)
axis equal tight
title("Quiver")

figure
contourf(flipud(squeeze(Rho_cv)),30)
axis equal tight
title("Density")

figure
contourf(flipud(squeeze(sum(U_cv.^2))),30)
axis equal tight
title("Velocity Mag")

% 
mid = round(N_cv_x/2);
Vertical_Sample = U_cv(1, :, mid)/U_top;
Horizontal_Sample = U_cv(2, mid, :)/U_top;
% 
% %u2_Ghia = [1 0.48223 0.46120 0.45992 0.46036 0.33556 0.20087 0.08183 -0.03039 -0.07404 -0.22855 -0.33050 -0.40435 -0.43643 -0.42901 -0.41165 0.00000 ];
% %v2_Ghia = [0.00000 -0.49774 -0.55069 -0.55408 -0.52876 -0.41442 -0.36214 -0.30018 0.00945 0.27280 0.28066 0.35368 0.42951 0.43648 0.43329 0.42447 0.00000];
% 

L = N_cv_x;

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
% 
% figure
% u = flip(squeeze(U(1, :, :)));
% v = squeeze(U(2, :, :));
% [startX, startY] = meshgrid(1:50:N_x, 1:50:N_y);
% verts = stream2(1:N_x,1:N_y,u,v,startX,startY);
% streamline(verts)