syms w

%% Given Data

kt = 100e3; %N-m/rad % Tortional Stiffness of Shafts

Ip_A = 0.03; % kg-m^2 % Mass Moment of Inertia of Disc-A
Ip_B = 0.02; % Mass Moment of Inertia of Disc-B
Ip_C = 0.02; % Mass Moment of Inertia of Disc-C

n_AB = 2; % Gear Ratio between A and B
n_AC = 2; % Gear Ratio between A and C

%% Overall Transfer Matrices

% For Gear System A
A = [1, 1/kt; 0, 1] * [1, 0; -Ip_A*w^2, 1]

% For Gear System B
B = [1, 1/kt; 0, 1] * [1, 0; -Ip_B*w^2, 1]

% For Gear System C
C = [1, 1/kt; 0, 1] * [1, 0; -Ip_C*w^2, 1]

%% Frequency Equation

% By writing transfer matrices and work done by trasmitted torque,
% Frequency equation for given system can be written as

Freq_Eqn = A(1,1)*B(2,2)*C(2,1)*n_AB^2 + A(1,1)*B(2,1)*C(2,2)*n_AC^2 + A(2,1)*B(2,2)*C(2,2)*(n_AB*n_AC)^2 == 0

omega = double(unique(abs(solve(Freq_Eqn, w))))