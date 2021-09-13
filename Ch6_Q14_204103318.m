clear; clc;

%% Given Data

rho= 0; % Density % kg/m^3
G = 0.8e11; % Modulus of Rigidity % N/m^2

d = 0.03; % Diameter of Shaft % m
d_h = 0.003; % Diameter of Hole % m
e = 0.006; % Eccentricity % m
J_s = pi * d^4 / 32; %m^4 % Polar Moment of Inertia of Area of Shaft
J_h = pi * d_h^4 / 32; %m^4 % Polar Moment of Inertia of Area of Hole
J = J_s - (J_h + (pi*(d_h*e)^2 / 4)) %m^4 % Polar Moment of Inertia of Area of Shaft with Hole
            % Using parallel axis theorem

L = 1; % Length of Shaft % m

m_d = 5; % Mass of Disc % kg
d_d = 0.1; % Diameter of Disc % m
Ip = (1 / 2) * m_d * (d_d / 2)^2; % Polar Mass Moment of Inertia of the Disc

%% Elements and Nodes

nele = 1; % No. of Elements

l = L / nele; % Length of each Element
% Here, length of each element is taken same

connect = zeros(nele, 3); % % 1st column is element number and last two columns are nodes of that element
for i = 1:nele
    connect(i,:) = [i, i, i+1];
end

coord = zeros(nele+1,2); % 1st column is node number and 2nd column is coordinate of that node
for i = 1:nele+1
    coord(i,:) = [i, (i-1)*l];
end

%% Boundary Conditions

bc = [4];
% 1: free-free ; 2:fixed-fixed ; 3:fixed-free ; 4:free-fixed
% For cantilevered shaft, disc is taken at left end, so free-fixed condition.

%% Mass and Stiffness Matrix Calculation

fig = figure('Name', 'Mode Shapes');
for p = 1:length(bc)
% Initializing Global Mass and Stiffness Matrices
M = zeros(nele+1, nele+1);
K = zeros(nele+1, nele+1);

% Elemental Mass and Stiffness Matrices
for i = 1:nele
    nd1 = connect(i,2);
    nd2 = connect(i,3);
    x1 = coord(nd1,2);
    x2 = coord(nd2,2);
    k_ele(:,:,i) = (G * J / l) * [1, -1; -1, 1];
    if(bc(p) == 4 && i == 1)
        m_ele(:,:,i) = [Ip, 0; 0, 0] + (rho * J * l / 6) * [2, 1; 1, 2]; % Disc at the left free end
    else
        m_ele(:,:,i) = (rho * J * l / 6) * [2, 1; 1, 2]; % No Ip terms as there is no any disc on the element
    end
    
    % Assembly
    vec = [nd1, nd2]; % Global DOF vector for assembly
    for ii = 1:2
        for jj = 1:2
            K(vec(ii),vec(jj)) = K(vec(ii),vec(jj)) + k_ele(ii,jj,i);
            M(vec(ii),vec(jj)) = M(vec(ii),vec(jj)) + m_ele(ii,jj,i);
        end
    end
end

%% Imposing Boundary Conditions

if(bc(p) == 1) % free-free
    K_red = K;
    M_red = M;
elseif(bc(p) == 2) % fixed-fixed
    K_red = K(2:nele,2:nele);
    M_red = M(2:nele,2:nele);
elseif(bc(p) == 3) % fixed-free
    K_red = K(2:nele+1,2:nele+1);
    M_red = M(2:nele+1,2:nele+1);
elseif(bc(p) == 4) % free-fixed
    K_red = K(1:nele,1:nele);
    M_red = M(1:nele,1:nele);
end

%% Natural Frequencies and Mode Shapes

% Eigen values and vectors
D = M_red \ K_red;
[eig_vec, eig_val] = eig(D);

% Natiral Frequency and Mode Shape
w_nf = [];
mode = [];
for i = 1:size(K_red,1)
    w_nf(i) = sqrt(eig_val(i,i)); % Natural Frequency
    mode(:,i) = eig_vec(:,i) / max(abs(eig_vec(:,i))); % Normalized Mode Shape
end

[w_nf,I] = sort(w_nf);
if(bc(p) == 1) % free-free
    disp('The below natural frequencies and modes are for free-free condition of the shaft.');
elseif(bc(p) == 2) % fixed-fixed
    disp('The below natural frequencies and modes are for fixed-fixed condition of the shaft.');
elseif(bc(p) == 3) % fixed-free
    disp('The below natural frequencies and modes are for fixed-free condition of the shaft.');
elseif(bc(p) == 4) % free-fixed
    disp('The below natural frequencies and modes are for free-fixed condition of the shaft.');
end
w_nf
mode = mode(I,:)

if(bc(p) == 2)
    mode = [zeros(1,size(mode,2)); mode; zeros(1,size(mode,2))];
elseif(bc(p) == 3)
    mode = [zeros(1,size(mode,2)); mode];
elseif(bc(p) == 4)
    mode = [mode; zeros(1,size(mode,2))];
end
%% Mode Shapes Plotting

for i = 1:length(w_nf)
    phi_z=[];
    x_z = [];
    for j = 1:size(connect,1)
        nd1 = connect(j,2);  
        nd2 = connect(j,3);
        x1 = coord(nd1,2);
        x2 = coord(nd2,2);
        l = x2-x1;
        x = x1:l/10:x2;     %deviding each element into small element for smooth mode shape
        phi1 = mode(j,i);
        phi2 = mode(j+1,i);
        for k=1:length(x)
            z = x(k)-x1;
            N1 = 1-(z/l);
            N2 = z/l;
            phi(k)=[N1 N2]*[phi1; phi2];
        end
        phi_z = [phi_z,phi];
        x_z =[x_z,x];
    end
        
    mode_1(:,i)=phi_z; % /max(abs(phi_z));  % normalised
end

if(size(mode_1,2) <= 4)
    nmode=size(mode_1,2);
else
    nmode=5;
end

set(gcf, 'Position', get(0,'Screensize'));
subplot(2,1,p)
for i = 1: nmode %length(wnf)
    if(i ==1)
        plot(x_z,mode_1(:,i), '-dm', 'DisplayName',['\omega=',num2str(w_nf(i))]);
    elseif (i ==2)
         plot(x_z,mode_1(:,i), '-sc', 'DisplayName',['\omega=',num2str(w_nf(i))]);
    elseif (i ==3)
         plot(x_z,mode_1(:,i),'-*r', 'DisplayName',['\omega=',num2str(w_nf(i))]);
    elseif (i ==4)
         plot(x_z,mode_1(:,i), '-+g', 'DisplayName',['\omega=',num2str(w_nf(i))]);
    else
         plot(x_z,mode_1(:,i), '-ob', 'DisplayName',['\omega=',num2str(w_nf(i))]);
    end
    hold on;
end

grid on;
if(bc(p) == 1) % free-free
    title('Free-free','fontsize',20);
elseif(bc(p) == 2) % fixed-fixed
    title('Fixed-fixed','fontsize',20);
elseif(bc(p) == 3) % fixed-free
    title('Fixed-free','fontsize',20);
elseif(bc(p) == 4) % free-fixed
    title('Free-fixed','fontsize',20);
end
xlabel('Shaft length(m)','fontsize',16);
ylabel('Relative amplitude','fontsize',16);
ylim([-1.2 1.2]);
legend('show');
end
saveas(fig,'mode_shape','png');