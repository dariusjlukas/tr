% Based on Manuel Hermann's
% "Torsional Stiffness and Natural Frequency analysis of a Formula SAE
% Vehicle Carbon Fiber Reinforced Polymer Chassis Using Finite Element
% Analysis" 
% Masters Thesis, California Polytechnic State University
% San Luis Obispo
% Also Robert D Story's "Design of Composite Sandwich Panels for a Formula
% SAE Monocoque Chassis"
% Masters Thesis, Oregon State University

% Transcribed by Thomas Lai, US Naval Academy/University of Maryland

% Classical Laminate Theory - 3-point-bend test

clear
clc
close all

format long

[num, txt, raw] = xlsread('PanelDataClean_metric_base_GFRMethod.xlsx');
ind = find(strcmp(raw(:,1),'CalPoly15_45F_0F_45U_-45U_0F_C_ditto'));
Mdata = raw(ind,:);

%% Materials Database
% All units in US customary (kg, m, s) 
   %  1,         2,         3,         4,         5,         6,         7,         8,         9,         10,        11,      
   %  E1_f,      E2_f,      v12_f,     G12_f,     F1t_f,     F1c_f,     F2t_f,     F2c_f,     F12_f,     t_f,       Aw_f   
M = [ Mdata(17), Mdata(18), Mdata(19), Mdata(20), Mdata(21), Mdata(22), Mdata(23), Mdata(24), Mdata(25), Mdata(15), Mdata(16); % Fabric
   %  E1_u,      E2_u,      v12_u,     G12_u,     F1t_u,     F1c_u,     F2t_u,     F2c_u,     F12_u,     t_u,       Aw_u   
      Mdata(33), Mdata(34), Mdata(35), Mdata(36), Mdata(37), Mdata(38), Mdata(39), Mdata(40), Mdata(41), Mdata(31), Mdata(32); % Uni
   %  E1_c,      E2_c,      v12_c,     G12_c,     GL,        GW,        FL,        FW,        Fc,        t_c,       Aw_c
      Mdata(46), Mdata(47), Mdata(48), Mdata(49), Mdata(50), Mdata(51), Mdata(52), Mdata(53), Mdata(54), Mdata(42), Mdata(44)]; % Core

M = cell2mat(M);
M(3,11) = M(3,11)*M(3,10); %calculate areal weight of honeycomb
M_cell_f = num2cell(M(1,:));
M_cell_u = num2cell(M(2,:));
M_cell_c = num2cell(M(3,:));

[E1_f,      E2_f,      v12_f,     G12_f,     F1t_f,     F1c_f,     F2t_f,     F2c_f,     F12_f,     t_f,       Aw_f] = deal(M_cell_f{:});
[E1_u,      E2_u,      v12_u,     G12_u,     F1t_u,     F1c_u,     F2t_u,     F2c_u,     F12_u,     t_u,       Aw_u] = deal(M_cell_u{:});
[E1_c,      E2_c,      v12_c,     G12_c,     GL,        GW,        FL,        FW,        Fc,        t_c,       Aw_c] = deal(M_cell_c{:});

width  = cell2mat(Mdata(55));  % Width of the panel in m
span = cell2mat(Mdata(56));  % Span of the test fixture in m


%% Laminate Definition
% Define angle, thickness, and material of each ply in the laminate [angle, thickness, matl#]
F = [t_f, 1];
U = [t_u, 2];
C = [t_c, 3];

    L = [45.0, F;
         0.0, F;
         45.0, U;
         -45.0, U;
         0.0, F;
         0.0, C;
         0.0, F;
         -45.0, U;
         45.0, U;
         0.0, F;
         45.0, F];

    % get number of plies
        n = size(L,1);
    % find the total thickness 
        thick = sum(L(:,2));
    % Calculate Areal Weight of Panel
        Aw_panel = (Aw_f*length(find(L(:,3)==1)) + Aw_u*length(find(L(:,3)==2)) + Aw_c*length(find(L(:,3)==3)));
        
    n_core = find(L(:,3)==3);    % locate core
    t_outer = sum(L(1:(n_core-1),2));
    t_inner = sum(L(((n_core+1):n),2));
    d = thick - t_outer/2 - t_inner/2;
    ts_ave = (t_outer + t_inner)/2;

        %% Setting up the ABD Matrix
        % Initializing sub-matricies

        h = zeros(n+1,1);
        A = zeros(3);
        B = zeros(3);
        D = zeros(3);

        % Form R matrix which relates engineering to tensor strain
            R = [1, 0, 0; 0, 1, 0; 0, 0, 2];

        % Locate the bottom of the first ply
            h(1) = -thick/2;
            imax = n + 1;
        
        % loop for rest of the ply distances from midsurface
            for i = 2:imax
                h(i) = h(i-1) + L(i-1,2);
            end
            
        % loop over each ply in the top half to integrate the ABD matricies
            for i = 1:(n_core-1)

                % ply material ID
                    mi = L(i,3);        % grab material number from the laminate definition
                    v21 = M(mi,2)*M(mi,3)/M(mi,1);
                    D_s = 1 - M(mi,3)*v21;
                    

                % Q12 matrix
                    Q = [M(mi,1)/D_s,           M(mi,3)*M(mi,2)/D_s,    0;
                         M(mi,3)*M(mi,2)/D_s,   M(mi,2)/D_s,            0;
                         0,                     0,                      M(mi,4)];
                     
                % ply angle in radians
                    theta = L(i,1)*pi/180;

                % form transformation matricies T1 for ply
                    T1 = [(cos(theta))^2,           (sin(theta))^2,         2*sin(theta)*cos(theta); 
                          (sin(theta))^2,           (cos(theta))^2,         -2*sin(theta)*cos(theta);
                          -sin(theta)*cos(theta),   sin(theta)*cos(theta),  (cos(theta))^2 - (sin(theta))^2 ];

                % form T2
                    T2 = R*T1*R^(-1); 

                % form Qxy
                    Qxy = T1^(-1)*Q*R*T1*R^(-1);
                % alt method?    
%                     U1 = 1/8*(3*Q(1,1) + 3*Q(2,2) + 2*Q(1,2) + 4*Q(3,3));
%                     U2 = 1/2*(Q(1,1) - Q(2,2));
%                     U3 = 1/8*(Q(1,1) + Q(2,2) - 2*Q(1,2) - 4*Q(3,3));
%                     U4 = 1/8*(Q(1,1) + Q(2,2) + 6*Q(1,2) - 4*Q(3,3));
% 
%                     Qg_xx = U1 + U2*cos(2*theta) + U3*cos(4*theta);
%                     Qg_xy = U4 - U3*cos(4*theta);
%                     Qg_yy = U1 - U2*cos(2*theta) + U3*cos(4*theta);
%                     Qg_xs = 1/2*U2*sin(2*theta) + U3*sin(4*theta);
%                     Qg_ys = 1/2*U2*sin(2*theta) - U3*sin(4*theta);
%                     Qg_ss = 1/2*(U1 - U4) - U3*sin(4*theta);
% 
%                     Qxy2 = [Qg_xx, Qg_xy, Qg_xs;
%                             Qg_xy, Qg_xx, Qg_ys;
%                             Qg_xs, Qg_ys, Qg_ss];

                % build up the laminate stiffness matrices
                    A = A + Qxy*(h(i+1) - h(i));
                    B = B + Qxy*(h(i+1)^2 - h(i)^2);    % 1/2 applied later
                    D = D + Qxy*(h(i+1)^3 - h(i)^3);    % 1/3 applied later 
                
            end
                            
        %% Build Stiffness and Compliance Matrix
        % Change diplay format
        format short e
        
        B = 1/2*B;
        D = 1/3*D;
        
        % Assemble the stiffness Matrix
            K = [A, B;
                 B, D];
         
        % Calculate the compliance matrix
            C = K^(-1);

        E_fs = 1/t_outer*(A(1,1) - (A(1,2)^2)/A(2,2));
        kpanel = (span^3/(24*E_fs*width*ts_ave*d^2) + span/(4*GL*width*thick))^(-1);
        kpanelnmm = kpanel/1000


    
    
    