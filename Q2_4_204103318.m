%% Ch2 Q4

% We are told to calculate natural frequencies of the rotor system.
% Natural frequnecies can be calculated by general formula...
% K - M*omega^2 == 0,
% where K = stiffness matrix
% and M = mass matrix

syms m I_d E I a b a_yf a_yM a_phif a_phiM omega
%% Given Data

disp("Given data:")
a = 0.3
b = 0.7
E = 2.1e11
d = .01
I = (pi/64)*d^4
I_d = 0.02
m = 5
%% Influence and mass matrices
% mass matrix = M

a_yf = (a^2)*(a+b)/(3*E*I);
a_yM = a*(3*a+2*b)/(6*E*I);
a_phif = a*(3*a+2*b)/(6*E*I);
a_phiM = (3*a+b)/(3*E*I);
influence_matrix = [a_yf a_yM; a_phif a_phiM]

M = [m 0;0 I_d]
%% Stiffness Matrix
% Stiffness matrix = K = inverse of influence mstrix

K = inv(influence_matrix)
%% Natural Frequencies calculation
% by general formula

eqn = det(K - M*omega^2) == 0
omega = double(solve(eqn))
%% Natural Frequencies
% As negative natural frequency means nothing
% Here critical speed is showed in form of natural frequencies

disp("The natural frequencies are:")
omega_1 = omega(1)
omega_2 = omega(2)