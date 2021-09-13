clear; clc;

%% Given Data

L1 = 0.0015; %mm
phiL1 = (pi/180)*(150); %rad
l1 = L1*(cos(phiL1) + 1i*sin(phiL1));
R1 = 0.0007;
phiR1 = (pi/180)*(162);
r1 = R1*(cos(phiR1) + 1i*sin(phiR1));

T_L = 29/1000; %kg
phiT_L = (pi/180)*(0);
T_l = T_L*(cos(phiT_L) + 1i*sin(phiT_L));
L2 = 0.0015;
phiL2 = (pi/180)*(110);
l2 = L2*(cos(phiL2) + 1i*sin(phiL2));
R2 = 0.0015;
phiR2 = (pi/180)*(170);
r2 = L1*(cos(phiR2) + 1i*sin(phiR2));

T_R = 29/1000;
phiT_R = (pi/180)*(0);
T_r = T_R*(cos(phiT_R) + 1i*sin(phiT_R));
L3 = 0.0025;
phiL3 = (pi/180)*(200);
l3 = L3*(cos(phiL3) + 1i*sin(phiL3));
R3 = 0.0015;
phiR3 = (pi/180)*(170);
r3 = R3*(cos(phiR3) + 1i*sin(phiR3));

%% Influence Coefficients

alpha_bR = (r3 - r1) / (T_r)
alpha_aR = (l3 - l1) / (T_r)
alpha_bL = (r2 - r1) / (T_l)
alpha_aL = (l2 - l1) / (T_l)

%% Correction Masses

C_M = inv([alpha_bR, alpha_bL; alpha_aR, alpha_aL])...
                * [r1; l1];
M_R = abs(C_M(1)) % Correction mass at right
phi_R = angle(C_M(1)) * (180/pi)
M_L = abs(C_M(2)) % Correction mass at left
phi_L = angle(C_M(2)) * (180/pi)