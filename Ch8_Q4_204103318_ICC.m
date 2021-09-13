clearvars
clc

% Given data
E = 2.1e11; % N/m^2
m = 5; %kg;
Ip = 0.0; %
Id = 0.02; %m
d = 0.01; %m
L = 1; %m
a = 0.3;
b = 0.7;


% Calculation based in given data
I = pi*d^4/64; % m^4

% Cantilever beam
a_11=(a^2)*(a+b)/(3*E*I);
a_12=(3*a*a+2*a*b)/(6*E*I);
a_21=a_12;
a_22=(3*a+b)/(3*E*I);


% Simply supported
%a_11 = a^2*b^2/(3*E*I*L) % m/N
%a_12 = -(3*a^2*L-2*a^3-a*L^2)/(3*E*I*L)
%a_21 = a*b*(b-a)/(3*E*I*L)
%a_22 = -(3*a*L-3*a^2-L^2)/(3*E*I*L)

M = [m;
     Id];
A = [a_11*m a_12*Id;
     a_21*m a_22*Id];
[eig_vec,eig_val] = eig(A);
mode = [];
wnf=[];
for i=1:size(eig_val,2)
    if eig_val(i,i)~=0
        wnf = [wnf,sqrt(1/eig_val(i,i))];
        mode = [mode, eig_vec(:,i)];%/max(abs(eig_vec(:,i)))];
    end
end

for i =1:length(wnf)-1      %arranging in asceending order
    if wnf(i)>wnf(i+1)
        temp_wnf = wnf(i);
        wnf(i) = wnf(i+1);
        wnf(i+1) = temp_wnf;
        temp_mode(:,1) = mode(:,i);
        mode(:,i) = mode(:,i+1);
        mode(:,i+1) = temp_mode(:,1);
    end
end
wnf