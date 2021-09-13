clear;
clc;

%% System Properties

m = 5; %kg
k = 10e3; %N/m
c = 5; %Ns/m

%% First Part 

k_sw = 1e3; %N/m
omega_nf = sqrt(k/m);
omega_sw = sqrt(k_sw/m);
damping_ratio = c/(2*sqrt(k/m));

% equation to evaluate whirl frequency
whirl_frequency = omega_sw^2/(2*damping_ratio*omega_nf);

stability_condition_1 = (2*damping_ratio*omega_nf > 0);
stability_condition_2 = (k_sw/k < 2*damping_ratio);

if(stability_condition_1 && stability_condition_2)
    fprintf("The system is stable for k_sw = %f N/m\n", k_sw)
else
    fprintf("The system is unstable for k_sw = %f N/m\n", k_sw)
end

fprintf("Whirl frequency = %f rad/s\n",whirl_frequency)
fprintf("Whirl frequency at boundary of stability = omega_nf = %f rad/s\n",omega_nf);

%% Second Part 

k_sw = 0.1; %N/m
omega_nf = sqrt(k/m);
omega_sw = sqrt(k_sw/m);
damping_ratio = c/(2*sqrt(k/m));

% equation to evaluate whirl frequency
whirl_frequency = omega_sw^2/(2*damping_ratio*omega_nf);

stability_condition_1 = (2*damping_ratio*omega_nf > 0);
stability_condition_2 = (k_sw/k < 2*damping_ratio);

if(stability_condition_1 && stability_condition_2)
    fprintf("\nThe system is stable for k_sw = %f N/m\n", k_sw)
else
    fprintf("The system is unstable for k_sw = %f N/m\n", k_sw)
end

fprintf("Whirl frequency = %f rad/s\n",whirl_frequency)
fprintf("Whirl frequency at boundary of stability = omega_nf = %f rad/s\n",omega_nf);
