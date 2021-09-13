clear; clc;
%% Symbols

syms phi_fp Torque_fp phi_rD Torque_rD w

%% Given Data

d_d = 0.1; %m % Diameter of Disc
m = 5; %kg % Mass of Disc
d = 0.03; %m % Diameter of Shaft
d_h = 0.003; %m % Diameter of Hole
e = 0.006; %m % Eccentricity
G = 0.8e11; %N/m^2 % Modulus of Rigidity
l = 1; %m % Length

%% Calculation of Required Quantities

Ip = (1/2) * m * (d_d / 2)^2 %kg-m^2 % Mass Moment of Inertia of Disc

J_s = pi * d^4 / 32; %m^4 % Polar Moment of Inertia of Area of Shaft
J_h = pi * d_h^4 / 32; %m^4 % Polar Moment of Inertia of Area of Hole
J = J_s - (J_h + (pi*(d_h*e)^2 / 4)) %m^4 % Polar Moment of Inertia of Area of Shaft with Hole
            % Using parallel axis theorem

kt = G * J / l %N-m/rad % Tortional Stiffness of Shaft with Hole

%% Field and Point Matrices

% Field Matrix
F = [1, 1/kt; 0, 1]

% Point Matrix
P = [1, 0; -Ip*w^2, 1]

%% Overall Transfer Matrix

T = P * F

%% Equation in Matrix Form

Eqn = [phi_rD; Torque_rD] == T * [phi_fp; Torque_fp]

%% Frequency Equation And Frequencies

% Applying Boundary Conditions
phi_fp = 0 % Fixed left side at shaft
Torque_rD = 0 % Free right side of disc

% From Boundary conditions and Eqn
Freq_Eqn = T(2,2) == 0

omega = double(unique(abs(solve(Freq_Eqn, w))))