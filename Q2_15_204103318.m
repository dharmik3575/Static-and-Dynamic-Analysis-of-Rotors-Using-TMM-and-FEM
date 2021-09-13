%% Ch2 Q15

% We are told to calculate natural frequencies of the rotor system.
% Natural frequnecies can be calculated by general formula...
% K - M*omega^2 == 0,
% where K = stiffness matrix
% and M = mass matrix

syms m I_d E I l a_yf a_yM a_phif a_phiM  omega
%% Influence and mass matrices
% mass matrix = M

a_yf = (l^3)/(48*E*I);
a_yM = 0;
a_phif = 0;
a_phiM = l/(12*E*I);
influence_matrix = [a_yf a_yM; a_phif a_phiM]

M = [m 0;0 I_d]
%% Stiffness Matrix
% Stiffness matrix = K = inverse of influence mstrix

K = inv(influence_matrix)
%% Natural Frequencies calculation
% by general formula

omega = sqrt(eig(K\M))
%% Natural Frequencies
% As negative natural frequency means nothing

omega_1 = omega(1)
omega_2 = omega(2)