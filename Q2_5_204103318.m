%% Ch2 Q5

% We are told to calculate the bearing force and moment at different...
% speeds.
% So, first natural frequencies are calculated and then force and moment...
% are calculated.
% The natural frequencies are calculated by
% K - M*omega^2 == 0
% where K is the stiffness matrix
% and M is the mass matrix

syms m me I_d E I l a_yf a_yM a_phif a_phiM omega omg y phi
%% Given Data

disp("Given data:")
l = 0.5
E = 2.1e11
d = .01
I = (pi/64)*d^4
I_d = 0.02
m = 10
me = 25e-5
%% Influence and mass matrix(M)

a_yf = l^3/(3*E*I);
a_yM = l^2/(2*E*I);
a_phif = l^2/(2*E*I);
a_phiM = l/(E*I);
influence_matrix = [a_yf a_yM; a_phif a_phiM]

M = [m 0;0 I_d]
%% Stiffness Matrix
% Stiffness matrix = K = inverse of influence matrix

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
%% Force and Moment Calculation
% Now we will use equation...
% M*xdd + F_l = F_unb
% where, M = mass matrix
% xdd = double differentiation of x = [y; phi]
% F_l = linear forse = [f_y; M_yz]... This is what we need to fine here.
% F_unb = unbalance force = [me*omega^2; 0], as there is no moment on
% rotor.
% For each case we will first find diflection [y; phi]...
% and then from those values we will fine force and moment.

F_unb = [me*omg^2; 0]

eqn1 = (K(1,1) - (m*omg^2))*y + K(1,2)*phi == me*omg^2
eqn2 = K(2,1)*y + (K(2,2) - (I*omg^2))*phi == 0
eqn = [eqn1, eqn2]

exp = F_unb - M*[y; phi]
%% Case 1: 0.5*omega_1

new_omg = 0.5 * omega_1
eqn_1 = subs(eqn, omg, new_omg)
exp_1 = subs(exp, omg, new_omg)

sol = solve(eqn_1, [y, phi])
Force_1 = double(subs(exp_1(1), sol))
Moment_1 = double(subs(exp_1(2), sol))
%% Case 2: 0.5*(omega_1+omega_2)

new_omg = 0.5 * (omega_1 + omega_2)
eqn_2 = subs(eqn, omg, new_omg)
exp_2 = subs(exp, omg, new_omg)

sol = solve(eqn_2, [y, phi])
Force_2 = double(subs(exp_2(1), sol))
Moment_2 = double(subs(exp_2(2), sol))
%% Case 3: 1.5*omega_2

new_omg = 1.5 * omega_2
eqn_3 = subs(eqn, omg, new_omg)
exp_3 = subs(exp, omg, new_omg)

sol = solve(eqn_3, [y, phi])
Force_3 = double(subs(exp_3(1), sol))
Moment_3 = double(subs(exp_3(2), sol))