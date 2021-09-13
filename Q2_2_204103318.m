%% Ch2 Q2
%% Calculations for amplitude ratio
syms xi omega omega_r omega_nfd omega_nf omega_rd phi
amp_r = omega_r^2/sqrt((1- omega_r^2)^2 + (2*xi*omega_r)^2)
omega_ratio = omega/omega_nf
omega_r_d = omega/omega_nfd
omega_r_d = subs(omega_r_d, omega_nfd, omega_nf*sqrt(1 - xi^2))
omega_r_d = subs(omega_r_d, omega/omega_nf, omega_r)
amp_r = subs(amp_r, omega_r, omega_rd*sqrt(1 - xi^2))

%% Calculation of the expression for damped frequency ratio
% Differentiating the amplitude ratio with respect to the damped frequency ratio
eq_1 = diff(amp_r, omega_rd) == 0
val = solve(eq_1, omega_rd) 

% This gives us three values of the damped frequency ratio
% out of these two are the same only there is difference of sign.

%% Calculations for phase
% only the value of frequency ratio is substituted to obtained the expression
tan(phi) = (2*xi*omega_r)/(1 - omega_r^2)
tan(phi) = subs(tan(phi), omega_r, omega_rd*sqrt(1- xi^2))

%% Plotting Nondimentional Amplitude

figure(1)
hold on
fplot(@(x) NDA(x, 0.05),[0 3])
fplot(@(x) NDA(x, 0.08),[0 3])
fplot(@(x) NDA(x, 0.1),[0 3])
fplot(@(x) NDA(x, 0.5),[0 3])
fplot(@(x) NDA(x, 1),[0 3])
fplot(@(x) NDA(x, 10),[0 3])
title("Nondimensional Amplitude vs Damped Frequency Ratio")
xlabel("\omega_{rd}")
ylabel("$\frac{Y}{e}$","Interpreter","latex")
legend(["\xi = 0.05","\xi = 0.08","\xi = 0.1", "\xi = 0.5","\xi = 1","\xi = 10"])
hold off
%% Plotting Phase

figure(2)
hold on
fplot(@(x) phase(x,0.05), [0 3])
fplot(@(x) phase(x,0.08), [0 3])
fplot(@(x) phase(x,0.1), [0 3])
fplot(@(x) phase(x,0.5), [0 3])
fplot(@(x) phase(x,1), [0 3])
fplot(@(x) phase(x,0.10), [0 3])
title("Phase vs Damped Frequency Ratio")
xlabel("\omega_{rd}")
ylabel("Phase")
legend(["\xi = 0.05","\xi = 0.08","\xi = 0.1", "\xi = 0.5","\xi = 1","\xi = 10"])
hold off
%% Function used to plot data
function out = den(omega_rd, xi)
    om = omega_rd.*(sqrt(1 - xi)^2);
    part_1 = (1 - om.^2).^2;
    part_2 = (2*xi*om).^2;
    out = (part_1 + part_2).^(1/2);
end

function out = num(omega_rd, xi)
    om = omega_rd.*(sqrt(1 - xi)^2);
    out = om.^2;
end

function out = NDA(omega_rd, xi)
    out = num(omega_rd,xi)./den(omega_rd,xi);
end

function out = phase(omega_rd, xi)
    omega = omega_rd.*sqrt(1 - xi^2);
    out = atan((2 * xi * omega)./(1 - omega.^2))*(180/(pi));
end

