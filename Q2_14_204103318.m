%% Ch2 Q14

% Here we need to find the diameter of the shaft for particular critical...
% speed.
% As we have calculated natural frequencies in the Q2_4_204103318.m...
% here we will go in the reverse direction. We will first find I...
% and then from I, we will get diameter of the shaft, as I = (pi/64)*d^4.
% Natural frequnecies can be calculated by general formula...
% K - M*omega^2 == 0,
% where K = stiffness matrix
% and M = mass matrix

syms m I_d E d I a b a_yf a_yM a_phif a_phiM omega
%% Given Data

disp("Given data:")
a = 0.35
b = 0.7
E = 2.1e11
m = 5
omega = 5.98
%% Influence and mass matrices
% mass matrix = M
% As here told to consider mass as a point mass...
% so only one influence coefficient will be there.
a_yf = (a^2)*(a+b)/(3*E*I);
influence = a_yf

M = m
%% Stiffness Matrix
% Stiffness matrix = K = inverse of influence mstrix

K = inv(influence)
%% Second moment of area of the shaft cross-section
% from general formula

eqn = det(K - M*omega^2) == 0
I = double(solve(eqn))
%% Diameter
% from I = (pi/64)*d^4

disp("The diameter of the shaft is")
d = (64*I/pi)^(1/4)