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
ind = find(contains(raw(:,1),'NR18_0F_45F_0U_45F_0F_C_0F_0U_0F'));
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
         
        % Inner
    L = [0.0, F;
         0.0, U;
         0.0, F;
         0.0, C;
         0.0, F;
         45.0, F;
         0.0, U;
         45.0, F;
         0.0, F];
        % Outer

    % get number of plies
        n = size(L,1);
    % find the total thickness 
        thick = sum(L(:,2));
    % Calculate Areal Weight of Panel
        Aw_panel = (Aw_f*length(find(L(:,3)==1)) + Aw_u*length(find(L(:,3)==2)) + Aw_c*length(find(L(:,3)==3)));
        
    n_core = find(L(:,3)==3);    % locate core
    t_inner = sum(L(1:(n_core-1),2));
    t_outer = sum(L(((n_core+1):n),2));
    d = thick - t_outer/2 - t_inner/2;
    ts_ave = (t_outer + t_inner)/2;


        %% Setting up the ABD Matrix
        % Initializing sub-matricies

        h = zeros(n+1,1);
        A = zeros(3,3,n);
        B = zeros(3,3,n);
        D = zeros(3,3,n);
        Ex_layer = zeros(n,1);
        h_bar_layer = zeros(n,1);
        
        % Form R matrix which relates engineering to tensor strain
            R = [1, 0, 0; 0, 1, 0; 0, 0, 2];

        % Locate the bottom of the inner ply
            h(1) = -thick/2;
            imax = n + 1;
        
        % loop for rest of the ply distances from midsurface
            for i = 2:imax
                h(i) = h(i-1) + L(i-1,2);
            end
            
        %% loop over each ply in the top half to integrate the ABD matricies
            for i = 1:n

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

                % build up the laminate stiffness matrices
                    A(:,:,i) = Qxy*(h(i+1) - h(i));
                    B(:,:,i) = Qxy*(h(i+1)^2 - h(i)^2);    % 1/2 applied later
                    D(:,:,i) = Qxy*(h(i+1)^3 - h(i)^3);    % 1/3 applied later 
                % Calculate E of the layer;
                    Ex_layer(i) = cos(L(i,1))*M(L(i,3),1);
                    h_bar_layer(i) = h(i) + L(i,2)/2;
            end
                            
        % Build Stiffness and Compliance Matrix        
            A_in = sum(A(:,:,(1:(n_core-1))),3);
            B_in = 1/2*sum(B(:,:,(1:(n_core-1))),3);
            D_in = 1/3*sum(D(:,:,(1:(n_core-1))),3);
            A_out = sum(A(:,:,((n_core+1):n)),3);
            B_out = 1/2*sum(B(:,:,((n_core+1):n)),3);
            D_out = 1/3*sum(D(:,:,((n_core+1):n)),3);
        
        % Assemble the Stiffness Matrix
            K_in = [A_in, B_in;
                    B_in, D_in];  
            K_out = [A_out, B_out;
                     B_out, D_out];
        
        % Calculate E of facesheets
            E_fs_in = 1/t_inner*(A_in(1,1) - (A_in(1,2)^2)/A_in(2,2));
            E_fs_out = 1/t_outer*(A_out(1,1) - (A_out(1,2)^2)/A_out(2,2));
        
        % Locate local centroids and global centroid
            h_bar_in = sum(Ex_layer((1:(n_core-1)),1).*L((1:(n_core-1)),2).*h_bar_layer((1:(n_core-1)),1))./(sum(Ex_layer((1:(n_core-1)),1).*L((1:(n_core-1)),2)));
            h_bar_out = sum(Ex_layer(((n_core+1):n),1).*L(((n_core+1):n),2).*h_bar_layer(((n_core+1):n),1))./(sum(Ex_layer(((n_core+1):n),1).*L(((n_core+1):n),2)));
            h_bar_tot = sum(Ex_layer((1:n),1).*L((1:n),2).*h_bar_layer((1:n),1))./(sum(Ex_layer((1:n),1).*L((1:n),2)));

        % Calculate I of each facesheet
            I_in = t_inner*width*abs(h_bar_in - h_bar_tot)^2;
            I_out = t_outer*width*abs(h_bar_out - h_bar_tot)^2; 
            kpanel_new = (span^3/(48*(E_fs_in*I_in + E_fs_out*I_out)) + span/(4*GL*width*thick))^(-1);
            kpanel_newnmm = kpanel_new/1000
            
            kpanel = (span^3/(48*E_fs_out*width*t_outer*d^2) + span^3/(48*E_fs_in*width*t_inner*d^2) + span/(4*GL*width*thick))^(-1); 
            kpanelnmm = kpanel/1000


    
    
    