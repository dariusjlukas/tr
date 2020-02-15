% Credit to Manuel Hermann, California Polytechnic State University
% San Luis Obispo
% From "Torsional Stiffness and Natural Frequency analysis of a Formula SAE
% Vehicle Carbon Fiber Reinforced Polymer Chassis Using Finite Element
% Analysis" 
% Masters Thesis

% Transcribed by Thomas Lai, US Naval Academy/University of Maryland

% Classical Laminate Theory - 3-point-bend test

clear
clc
close all

format long
%% Toggle Switch

toggle = 0;

% toggle = 0: evaluate single value only
% toggle = 1: plot core thickness graph
% toggle = 2: plot skin thickness graph

%% Materials Database
% All units in US customary (lb, in, s) 

   % Expired Cycom 5320-1 Fabric from boeing 
   % 1,    2,                          3,                          4,                         5,     6,     7,     8,    9,    10,   11,   12,   13,   14,   15,   16,   17 
   % rho,  E1,                         E2,                         G12,                       v12,   a_11,  a_22,  F1t,  F2t,  F12,  F1c,  F2c,  e1t,  e2t,  e12,  e1c,  e2c
M = [0,    0.7*70.81*10^9*0.000145038, 0.7*69.96*10^9*0.000145038, 0.7*5.58*10^9*0.000145038, 0.048, 0,     0,     0,    0,    0,    0,    0,    0,    0,    0,    0,    0; 
   
   % Expired Cycom 5320-1 uni from boeing 
   % 1,    2,                           3,                          4,                          5,    6,     7,     8,    9,    10,   11,   12,   13,   14,   15,   16,   17 
   % rho,  E1,                          E2,                         G12,                        v12,  a_11,  a_22,  F1t,  F2t,  F12,  F1c,  F2c,  e1t,  e2t,  e12,  e1c,  e2c
     0,    0.7*157.13*10^9*0.000145038, 0.7*9.373*10^9*0.000145038, 0.7*5.513*10^9*0.000145038, 0.3,  0,     0,     0,    0,    0,    0,    0,    0,    0,    0,    0,    0; 
   
   % PAMG-XR1-3.4-1/4-N-5056 core
   % 1,      2,     3,     4,    5,    6,       7,       8,       9,       10,  11,  12,  13,  14,  15,  16,  17   
   % rho,    E1,    E2,    G12,  v12,  G23(W),  G23(L),  F23(W),  F23(L),  F3,  F1c, F2c, e1t, e2t, e12, e1c, e2c
     0,      89.62, 46.80, 4.46, 1,    22000,   52000,   0,       0,       0,   0,   0,   0,   0,   0,   0,   0]; 
   
   % rho = density
   % E1 = Young's Modulus in 1 direction
   % E2 = Young's Modulus in 2 direction
   % G12 = Shear Modulus
   % v12 = Major Poisson's Ratio
   % G23 (W) = G23 = Shear Modulus transverse to the ribbon (core) % not used in code
   % G23 (L) = G13 = Shear Modulus in ribbon direction (core) % not used in code
   % a_11 = thermal expansion coefficient in the 1 direction
   % a_22 = thermal expansion coefficient in the 2 direction
   % F1t = ? % not used in code
   % F2t = ? % not used in code
   % F12 = ? % not used in code
   % F23 (W) = ?% not used in code
   % F23 (L) = ?% not used in code
   % F3 = ?% not used in code
   % F1c = ?% not used in code
   % F2c = ?% not used in code
   % e1t = max tensile strain in 1 direction
   % e2t = max tensile strain in 2 direction
   % e12 = max shear strain
   % e1c = max compressive strain in 1 direction
   % e2c = max compressive strain in 2 direction
   
   % rho only used for mass of panel
   % max strains only used for failure theory
   % To get stiffness with no thermal stress, only E1, E2, G12, v12 are needed
   
%% Laminate Definition
% Define angle, thickness, and material of each ply in the laminate [angle, thickness, matl#]                           
if toggle == 0
    cthick = 5/8;               % core thickneess
    pthick = [0.476/25.4;       % fabric ply thickness 
              0.15875/25.4];    % uni ply thickness
end
if toggle == 1
    cthick = 5/8*linspace(0.8,1.2,10);     % core thickness
    pthick = [0.476/25.4;                  % ply thickness
              0.15875/25.4];                                                               
end
if toggle == 2        
    cthick = 5/8;                                   % core thickness                                 
    pthick = [0.476/25.4;                           % ply thickness
              0.15875/25.4].*linspace(0.8,1.2,10);   
end

% Initialize Vectors
if toggle == 0
    mass = NaN(1,length(cthick));
    dtotaln = NaN(1,length(cthick));
    dtotalw = NaN(1,length(cthick));
    kpaneln = NaN(1,length(cthick));
    kpanelw = NaN(1,length(cthick));
    cshear = NaN(1,length(cthick));
    cshear_thesis = NaN(1,length(cthick));
    sigma_c_el = NaN(1,length(cthick));
end
if toggle == 1
    mass = NaN(1,length(cthick));
    dtotaln = NaN(1,length(cthick));
    dtotalw = NaN(1,length(cthick));
    kpaneln = NaN(1,length(cthick));
    kpanelw = NaN(1,length(cthick));
    cshear = NaN(1,length(cthick));
    cshear_thesis = NaN(1,length(cthick));
    sigma_c_el = NaN(1,length(cthick));
end
if toggle == 2
    mass = NaN(1,length(pthick(1,:)));
    dtotaln = NaN(1,length(pthick(1,:)));
    dtotalw = NaN(1,length(pthick(1,:)));
    kpaneln = NaN(1,length(pthick(1,:)));
    kpanelw = NaN(1,length(pthick(1,:)));
    cshear = NaN(1,length(cthick));
    cshear_thesis = NaN(1,length(pthick(1,:)));
end

% Sweep through possible core thicknesses
for k = 1:length(cthick)
    % Sweep through possible ply thicknesses
    for kk = 1:length(pthick(1,:))
    
        L = [ 0, pthick(1,kk), 1;
              45, pthick(1,kk), 1;
              0, pthick(1,kk), 1;
              0, cthick(k), 3; %% Core---
              0, pthick(1,kk), 1;
              0, pthick(2,kk), 2;
              0, pthick(2,kk), 2;
              0, pthick(2,kk), 2;
              0, pthick(1,kk), 1 ];
              
        % material #1: Hexcel IM7 Fabric
        % material #2: Expired Cycom 5320-1 uni from boeing 
        % material #3: PAMG-XR1-e.4-1/4-N-5056 core
       
        % get number of plies
            n = size(L,1);
        % find the total thickness 
            thick = sum(L(:,2));
        
        %% Load Input
        % mechanical line loads specifically for 3-point-bend test w/ central load & simple support at both ends

        Force = 8705/4.44822;    % applied force in lb
        span = 400/25.4;            % Span of the test fixture in inches
        width  = 275/25.4;          % Width of the panel in inches
        Mx = Force*span/(4*width);
        
        % Nx; Ny; Nxy; Mx; My; Mxy
            load_m = [0; 0; 0; Mx; 0; 0];
        % thermal loads
            delT = 0;   % temperature change in Farenheit

        %% Setting up the ABD Matrix
        % Initializing sub-matricies

        h = zeros(n+1,1);
        A = zeros(3);
        B = zeros(3);
        D = zeros(3);
        N_t = zeros(3,1);
        M_t = zeros(3,1);
        FIF = zeros(1,3);

        % Form R matrix which relates engineering to tensor strain
            R = [1, 0, 0; 0, 1, 0; 0, 0, 2];

        % Locate the bottom of the first ply
            h(1) = -thick/2;
            imax = n + 1;
        
        % loop for rest of the ply distances from midsurface
            for i = 2:imax
                h(i) = h(i-1) + L(i-1,2);
            end

        % initialize ply mass
            ply_mass = NaN(1,n);
            
        % loop over each ply to integrate the ABD matricies
            for i = 1:n

                % ply material ID
                    mi = L(i,3);        % grab material number from the laminate definition
                    v21 = M(mi,3)*M(mi,5)/M(mi,2);
                    D_s = 1 - M(mi,5)*v21;

                % Q12 matrix
                    Q = [M(mi,2)/D_s,           v21*M(mi,2)/D_s,    0;
                         M(mi,5)*M(mi,3)/D_s,   M(mi,3)/D_s,        0;
                         0,                     0,                  M(mi,4)];

                % ply angle in radians
                    theta = L(i,1)*pi/180;

                % form transformation matricies T1 for ply
                    T1 = [(cos(theta))^2,           (sin(theta))^2,         2*sin(theta)*cos(theta); 
                          (sin(theta))^2,           (cos(theta))^2,         -2*sin(theta)*cos(theta);
                          -sin(theta)*cos(theta),   sin(theta)*cos(theta),  (cos(theta))^2 - (sin(theta))^2 ];

                % form T2 for thermal
                    T2 = R*T1*R^(-1); 

                % form Qxy
                    Qxy = T1^(-1)*Q*R*T1*R^(-1);

                % build up the laminate stiffness matrices
                    A = A + Qxy*(h(i+1) - h(i));
                    B = B + Qxy*(h(i+1)^2 - h(i)^2);    % 1/2 applied later
                    D = D + Qxy*(h(i+1)^3 - h(i)^3);    % 1/3 applied later 

                % load alphas into an array
                    a = [M(mi,6); M(mi,7); 0];

                % Calculate thermal loads
                    N_t = N_t + delT*Qxy*T2^(-1)*a*(h(i+1) - h(i));
                    M_t = M_t + delT*Qxy*T2^(-1)*a*(h(i+1)^2 - h(i)^2);

                % Calculate mass of the panel
                    ply_mass(i) = M(mi,1)*L(i,2)*width*span;          
            end
            
            if toggle == 0
                   mass(k) = sum(ply_mass);
            end
            if toggle == 1
                   mass(k) = sum(ply_mass);
            end
            if toggle == 2
                   mass(kk) = sum(ply_mass);
            end

        %% Build Stiffness and Compliance Matrix
        % Change diplay format
        format short e
        
        B = 1/2*B;
        D = 1/3*D;
        M_t = (1/2*M_t);

        % Assemble the stiffness Matrix
            K = [A, B;
                 B, D];
         
        % Calculate the compliance matrix
            C = K^(-1);

        %% Calculate the Deflection of the Beam
        % thickness to middle of facesheets (distance is h in Hexcel design guide)        
            if toggle == 0
                distance = (thick + cthick(k))/2;
            end
            if toggle == 1
                distance = (thick + cthick(k))/2;
            end   
            if toggle == 2
                distance = (thick + cthick)/2;
            end
            
        % wide and narrow Flexural Rigidity
            EIn = width/C(4,4); %<- panel EI for narrow
            EIw = width*K(4,4) %<- panel EI for wide

        % bend stiffness
            kbn = 48*EIn/span^3;
            kbw = 48*EIw/span^3

        % shear stiffness
            ks = width*distance*M(3,6);

        % shear deflection
            ds = 0.25*Force*span/ks;

        % bending deflection narrow and wide beam
            dbn = Force/kbn;
            dbw = Force/kbw;

        if toggle == 0
            % total deflection narrow and wide
            dtotaln(k) = ds + dbn;
            dtotalw(k) = ds + dbw;
        end        
        if toggle == 1
            % total deflection narrow and wide
            dtotaln(k) = ds + dbn;
            dtotalw(k) = ds + dbw;
        end
        if toggle == 2
            % total deflection narrow and wide
            dtotaln(kk) = ds + dbn;
            dtotalw(kk) = ds + dbw;
        end

        if toggle == 0
            % calculate specific panel stiffness
            kpaneln(k) = (kbn + ks)/mass(k);
            kpanelw(k) = (kbw + ks)/mass(k);
        end
        if toggle == 1
            % calculate specific panel stiffness
            kpaneln(k) = (kbn + ks)/mass(k);
            kpanelw(k) = (kbw + ks)/mass(k);
        end
        if toggle == 2
            % calculate specific panel stiffness
            kpaneln(kk) = (kbn + ks)/mass(kk);
            kpanelw(kk) = (kbw + ks)/mass(kk);
        end

        % calculate core shear stress
            cshear(k) = 0.5*Force/(distance*width);   % core shear using hexcel formula in psi

        if toggle == 0 
            cshear_thesis(k) = 0.5*Force/((thick + cthick(k))*width);   % core shear using simplified beam theory proposed by Allen in psi
        end
        if toggle == 1
            cshear_thesis(k) = 0.5*Force/((thick + cthick(k))*width); % core shear ising simplified beam theory proposed by Allen in psi
        end
        if toggle == 2
            cshear_thesis(kk) = 0.5*Force/((thick + cthick)*width); % core shear using simplified beam theory proposed by Allen in psi
        end

        %% Calculate elastic buckling stress
        cell_size = 1/4; % in in
        G_alu = 3.76*10^(6); % in psi
        v_alu = 0.33;

        if toggle == 0
            sigma_c_el(k) = 12.1/(1 - v_alu)*G_alu*(cthick(k)/cell_size)^3;  % equation 11
            % source: Materials Science and Engineering A380
            % "Out-of-plane compressive behaviour of metallic honeycombs"
            % F. Cote, V.S. Deshpande, N.A. Fleck, A.G. Evans
        end
        if toggle == 1
            sigma_c_el(k) = 12.1/((1 - v_alu)*G_alu*(cthick(k)/cell_size)^3);
        end
        if toggle == 2
            sigma_c_el = 12.1/(1 - v_alu)*G_alu*(cthick/cell_size)^3;      
        end

        %% Build Load Vector
        % build thermal load vector
            load_t = [N_t; M_t];
        % combine mechanical and thermal loading
            load = load_t + load_m;

        %% Calculate Strains, Curvatures, and Stresses
        % solve for strain and curvatures

        e = C*load;

        % initialize Stress/Strain Vectors
            stress12 = NaN(2*n,4);
            strain12 = NaN(2*n,4);
            stressxy = NaN(2*n,4);
            strainxy = NaN(2*n,4);
            efail = NaN(2*n,4);
        

        % loop over each ply and calculate strain
        for i = 1:n
            % loop over top and bottom of each ply starting at the top of the ply       
            for j = 1:2
                ply = i; % ply number
                loc = j; % location of the ply
                         % 1 = top, 2 = bottom

                % ply material ID
                    mi = L(i,3);
                    v21 = M(mi,3)*M(mi,5)/M(mi,2);
                    D_s = 1 - M(mi,5)*v21;

                % Q12 matrix
                    Q = [M(mi,2)/D_s,           v21*M(mi,2)/D_s,    0;
                         M(mi,5)*M(mi,3)/D_s,   M(mi,3)/D_s,        0;
                         0,                     0,                  M(mi,4)];

                % ply angle in radians
                    theta = L(i,1)*pi/180;

                % Form transformation matricies T1 for each ply
                    T1 = [(cos(theta))^2,           (sin(theta))^2,         2*sin(theta)*cos(theta);
                          (sin(theta))^2,           (cos(theta))^2,         -2*sin(theta)*cos(theta);
                          -sin(theta)*cos(theta),   sin(theta)*cos(theta),  (cos(theta))^2 - (sin(theta))^2];

                % Form T2 for thermal
                    T2 = R*T1*R^(-1);
                % Load alphas into an array
                    a = [M(mi,6); M(mi,7); 0]; % might be wrong indicies
                    e_thermal = delT*T2^(-1)*a;
                % Form Qxy
                    Qxy = T1^(-1)*Q*R*T1*R^(-1);
                % calculate distance from midsurface
                    z = h(i-1+j);
                % transform back into principle directions
                % global total strains 
                    e_tot_glo = [e(1) + z*e(4); e(2) + z*e(5); e(3) + z*e(6)];
                % local total strains
                    e_tot_loc = T2*e_tot_glo;
                % global mechanical strain
                    e_m_glo = e_tot_glo - e_thermal;
                % local mechanical strain
                    e_m_loc = e_tot_loc - e_thermal;

                % stresses only arise from mechanical strains
                % calculating stresses
                % global stresses
                    stxy = Qxy*e_m_glo;
                % local stresses
                    st12 = Q*e_m_loc;

                % uses MAX Strain criteria
                % failure index now looks at two different materials
                % check fiber direction
                
                if e_tot_loc(1) > 0
                    FI = e_tot_loc(1)/M(mi,13); %can probably clean this up a bit
                    FIF(1) = FI;
                elseif e_tot_loc(1) < 0
                    FI = abs(e_tot_loc(1))/M(mi,16);
                    FIF(1) = FI;
                end

                % check transverse direction
                if e_tot_loc(2) > 0
                    FIF(2) = e_tot_loc(2)/M(mi,14);
                elseif e_tot_loc(2) < 0
                    FIF(2) = abs(e_tot_loc(2))/M(mi,17);
                end

                % Check in-plane shear
                    FIF(3) = abs(e_tot_loc(3))/M(mi,15);

                % FIF is index on fiber failure
                % FIe is the highest failure index which could be fiber, transverse or shear

                % display stresses and stains in a matrix
                % Write local Stress for each ply on top and bottom into a matrix
                    stress12(2*i + j - 2,1) = i;                % ply number
                    stress12(2*i + j - 2,2) = L(i);             % ply angle
                    stress12(2*i + j - 2,3) = st12(1);          % stress in 1 direction
                    stress12(2*i + j - 2,4) = st12(2);          % stress in 2 direction
                    stress12(2*i + j - 2,5) = st12(3);          % Shear stress in 12
                % Write total local Strain for each ply on top and bottom into a matrix
                    strain12(2*i + j - 2,1) = i;                % ply number
                    strain12(2*i + j - 2,2) = L(i);             % ply angle
                    strain12(2*i + j - 2,3) = e_tot_loc(1);     % strain in ply 1 direction
                    strain12(2*i + j - 2,4) = e_tot_loc(2);     % strain in ply 2 direction
                    strain12(2*i + j - 2,5) = e_tot_loc(3);     % strain in ply 12 or shear strain
                % Write global Stress for each ply on top and bottom into a matrix
                    stressxy(2*i + j - 2,1) = i;                % ply number
                    stressxy(2*i + j - 2,2) = L(i);             % ply angle
                    stressxy(2*i + j - 2,3) = stxy(1);          % stress in 1 direction
                    stressxy(2*i + j - 2,4) = stxy(2);          % stress in 2 direction
                    stressxy(2*i + j - 2,5) = stxy(3);          % Shear stress in 12
                % Write total global Strain for each ply on top and bottom into a matrix
                    strainxy(2*i + j - 2,1) = i;                % ply number
                    strainxy(2*i + j - 2,2) = L(i);             % ply angle
                    strainxy(2*i + j - 2,3) = e_tot_glo(1);     % stress in 1 direction 
                    strainxy(2*i + j - 2,4) = e_tot_glo(2);     % stress in 2 direction
                    strainxy(2*i + j - 2,5) = e_tot_glo(3);     % Shear stress in 12
                % Write Failure Index into a matrix
                    efail(2*i + j - 2,1) = i;                   % ply number
                    efail(2*i + j - 2,2) = L(i);                % ply angle
                    efail(2*i + j - 2,3) = FIF(1);              % failure index longitudinal
                    efail(2*i + j - 2,4) = FIF(2);              % failure index transverse
                    efail(2*i + j - 2,5) = FIF(3);              % failure index shear

            end
        end
    end
end

%% Print out the results
disp('The deflection of the panel for a narrow beam solution: ')
fprintf('%2.3f \n', dtotaln)
disp('The deflection of the panel for a wide beam solution: ')
fprintf('%2.3f \n', dtotalw)
disp('The D16 of the Laminate: ')
fprintf('%2.3f \n', D(1,3))
disp('The D26 of the Laminate: ')
fprintf('%2.3f \n', D(2,3))

%% Plotting Core Thickness Graphs

if toggle == 1
    %% Plot Deflection and Specific Stiffness
    figure('Name', 'Stiffness Sensitivity', 'NumberTitle', 'Off')
    % yyaxis left
    plot(cthick, dtotalw, 'b')
    hold on
    xlabel('Core Thickness [in]')
    ylabel('Panel Deflection [in]')
    % yyaxis right
    % plot(cthick, kpanelw, 'r')
    % ylabel('Specific Stiffness [lb*in2/lbm]')
    % legend('Panel Deflection', 'Specific Stiffness', 'Location', 'NorthWest')
    % hold off

    %% Plot Elastic Buckling and Core Shear
    figure('Name', 'Strength Sensitivity', 'NumberTitle', 'Off')
    yyaxis left
    plot(cthick, sigma_c_el, 'b')
    hold on
    xlabel('Core Thickess [in]')
    ylabel('Core Elastic Buckling Stress [psi]')
    yyaxis right
    plot(cthick, cshear, 'r')
    ylabel('Core Shear Stress [psi]')
    legend('Elastic Buckling Stress', 'Core Shear Stress', 'Location', 'NorthWest')
    hold off
    
end

%% Plotting Ply Thickness Graph
if toggle == 2
    
    %% Plot Deflection
    figure('Name', 'Stiffness Sensitivity', 'NumberTitle', 'off')
    % yyaxis left
    plot(pthick(1,:), dtotalw, 'b')
    hold on
    xlabel('Ply Thickness [in]')
    ylabel('Panel Deflection [in]')
    % yyaxis right
    % plot(pthick, kpanelx, 'r')
    % ylabel('Specific Stiffness [lb*in2/lbm]')
    % legend('Panel Deflection', 'Specific Stiffness', 'Location', 'NorthWest')
    % hold off
    
end



    

    
    
    