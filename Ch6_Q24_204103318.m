clear; clc;

%% Given Data

rho= 0; % Density % kg/m^3
kt = 100e3; % Tortional Stiffness of Shafts % N-m/rad
Ip = [0.03, 0.02, 0.02]; % Polar Mass Moment of Inertia of the Discs
gr = [2, 2]; % Gear Ratios
%% Elements and Nodes

nele = 3; % No. of Elements

connect = [1, 1, 2
           2, 2, 3
           3, 2, 4]; % 1st column is element number and last two columns are nodes of that element

%% Boundary Conditions

bc = [1]; % for every elements
% 1: free-free ; 2:fixed-fixed ; 3:fixed-free ; 4:free-fixed

%% Mass and Stiffness Matrix Calculation

% Initializing Global Mass and Stiffness Matrices
M = zeros(nele+1, nele+1);
K = zeros(nele+1, nele+1);

% Elemental Mass and Stiffness Matrices
m_ele(:,:,1) = [Ip(1), 0; 0, 0];
m_ele(:,:,2) = [0, 0; 0, Ip(2)];
m_ele(:,:,3) = [0, 0; 0, Ip(3)];

k_ele(:,:,1) = kt * [1, -1; -1, 1];
k_ele(:,:,2) = kt * [1/gr(1)^2, 1/gr(1); 1/gr(1), 1];
k_ele(:,:,3) = kt * [1/gr(2)^2, 1/gr(2); 1/gr(2), 1];

for i = 1:nele
    nd1 = connect(i,2);
    nd2 = connect(i,3);
    % Assembly
    vec = [nd1, nd2]; % Global DOF vector for assembly
    for ii = 1:2
        for jj = 1:2
            K(vec(ii),vec(jj)) = K(vec(ii),vec(jj)) + k_ele(ii,jj,i);
            M(vec(ii),vec(jj)) = M(vec(ii),vec(jj)) + m_ele(ii,jj,i);
        end
    end
end

%% Imposing Boundary Conditions

if(bc(1) == 1) % free-free
    K_red = K;
    M_red = M;
elseif(bc(p) == 2) % fixed-fixed
    K_red = K(2:nele,2:nele);
    M_red = M(2:nele,2:nele);
elseif(bc(p) == 3) % fixed-free
    K_red = K(2:nele+1,2:nele+1);
    M_red = M(2:nele+1,2:nele+1);
elseif(bc(p) == 4) % free-fixed
    K_red = K(1:nele,1:nele);
    M_red = M(1:nele,1:nele);
end

%% Natural Frequencies and Mode Shapes

syms w

eqn = det(K_red - w^2 * M_red) == 0;

% Natiral Frequency and Mode Shape
w_nf = sort(double(unique(abs(solve(eqn, w)))))
