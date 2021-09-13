clearvars;
clc;

%% input
E = 2.1e+11;
d = 0.015;        %shaft dia
rho = 0;        %shaft density

md1 = 5; % in kg
Id1 = 0; %md1*D1^2/16; % in kg-m^2

md2 = 8;
Id2 = 0;

mass = [1,   md1,   Id1;
        3,   md2,   Id2];  %1st column node 2nd colum mass 3rd column Id

nele = 3;
connect = [1, 1, 2;
           2, 2, 3;
           3, 3, 4];   % First Column is element number
                  % Second & Third Column are Nodes (in sequence)
coord = [1, 0;       % first column node second column x coodinate
         2, 0.15;
         3, 0.30;
         4, 0.45];
       

bc_data = [2, 1;     % First Column is Node number
           4, 1];     % Second column is prescribed local dof
                    % dof 1 for v and 2 for phi

% Mass an stifness matrix

M = zeros(2*(nele+1));
K = zeros(2*(nele+1));

for i= 1:nele
    nd1 = connect(i,2);
    nd2 = connect(i,3);
    x1  = coord(nd1,2);
    x2  = coord(nd2,2);
    l = x2-x1;
    A=pi*d^2/4;
    m(:,:,i) = rho*A*l/420*[156     22*l    54      -13*l;
                            22*l    4*l^2   13*l    -3*l^2;
                            54      13*l    156     -22*l;
                            -13*l   -3*l^2  -22*l   4*l^2];
    if size(mass)>0                    
    if any(mass(:,1)==nd1)      %checking for disc at left node of element
        ii = find(mass(:,1)==nd1);
        m(1,1,i) = m(1,1,i) + mass(ii,2);
        m(2,2,i) = m(2,2,i) + mass(ii,3);
    end
    if i==nele                  % checking for disc at right side of last element
        if any(mass(:,1)==nd2)
            ii = find(mass(:,1)==nd2);
            m(3,3,i) = m(3,3,i) + mass(ii,2);
            m(4,4,i) = m(4,4,i) + mass(ii,3);
        end
    end
    end
       
    I =pi*d^4/64;
   
    k(:,:,i) = E*I/l^3 *[12     6*l     -12     6*l;
                         6*l    4*l^2   -6*l    2*l^2;
                         -12    -6*l    12      -6*l;
                         6*l    2*l^2   -6*l    4*l^2];
                     
    %Assembly
   
    vec = [2*nd1-1, 2*nd1, 2*nd2-1, 2*nd2];
    M(vec,vec) = M(vec,vec)+ m(:,:,i);
    K(vec,vec) = K(vec,vec)+ k(:,:,i);
end

% imposition of boundary condition
supp_dof = [];      % dof to be suppressed.
for i = 1:size(bc_data,1)
    nd = bc_data(i,1);
    dof = bc_data(i,2);
    Gdof = 2*(nd-1)+dof;        %finding the global DOF of given BC
    supp_dof = [supp_dof , Gdof];
end

supp_dof = sort(supp_dof);


% Reducing Global matrix (Imposition of Boundary condition)
% ---------------------------------------------------------
K_red = K;
M_red = M;
for i = 1:size(supp_dof,2)
    dof = supp_dof(i);  %Getting the Global DOF of Given BC
    if (dof == 1)       %Loop for reducing global matrix
        K_red = K_red(dof+1:end, dof+1:end);
        M_red = M_red(dof+1:end, dof+1:end);
    elseif (dof == 2*(nele+1))
            K_red = K_red(1:dof-1, 1:dof-1);
            M_red = M_red(1:dof-1, 1:dof-1);
    else
        K_red = K_red([1:dof-1, dof+1:end],[1:dof-1, dof+1:end]);
        M_red = M_red([1:dof-1, dof+1:end],[1:dof-1, dof+1:end]);
    end
    if (i~=size(supp_dof,2))    %for considering the already reduced steps
        supp_dof(i+1:end) = supp_dof(i+1:end)-1;
    end
end

D = K_red\M_red;
[eig_vec,eig_val] = eig(D);
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
% updating boundary condition

for i = 1:size(bc_data,1)
    nd = bc_data(i,1);
    dof = bc_data(i,2);
    Gdof = 2*(nd-1)+dof;    %Finding global DOF
    if Gdof ==1
        mode = [zeros(1,size(mode,2));mode];
    else
        mode = [mode((1:Gdof-1),:);zeros(1,size(mode,2));mode((Gdof:end),:)];
    end
end

% Post processing


for i =1:length(wnf)
    x = [];
    v=[];
    phi=[];
    for ii= 1:nele
        nd1 = connect(ii,2);
        nd2 = connect(ii,3);
        x1  = coord(nd1,2);
        x2  = coord(nd2,2);
        l = x2-x1;
        x_n =x1:(x2-x1)/15:x2;
        v_n=[];
        phi_n=[];
        eta = mode(2*nd1-1:2*nd2,i);
        for j =1:length(x_n)
            z = x_n(j)-x1;
            N(1) = 1-(3*z^2/l^2)+(2*z^3/l^3);
            N(2) = z-(2*z^2/l)+(z^3/l^2);
            N(3) = (3*z^2/l^2)-(2*z^3/l^3);          
            N(4) = -(z^2/l)+(z^3/l^2);
            v_n(j) = N*eta;
            dN(1) = (-6*z/l^2)+(6*z^2/l^3);
            dN(2) = 1-(4*z/l)+(3*z^2/l^2);
            dN(3) = (6*z/l^2)-(6*z^2/l^3);
            dN(4) = (-2*z/l)+(3*z^2/l^2);
            phi_n(j) = dN*eta;
        end
        x = [x,x_n];
        v = [v,v_n];
        phi = [phi,phi_n];
    end
    V(:,i) = v'/max(abs(v));
    phi_x(:,i)  = phi'/max(abs(phi));

end

% ploting mode shape upto 4th mode

h = figure(1);
set(gcf, 'Position', get(0,'Screensize'));
subplot(1,2,1)

if(size(V,2) < 4)
    nmode = size(V,2);
else
    nmode = 4;
end


for i = 1:nmode
    if(i ==1)
        plot(x,V(:,i), '-k', 'LineWidth', 2, 'DisplayName',['\omega=',num2str(wnf(i))]);
    elseif (i ==2)
         plot(x,V(:,i), ':k', 'LineWidth', 2, 'DisplayName',['\omega=',num2str(wnf(i))]);
    elseif (i ==3)
         plot(x,V(:,i),'-.k', 'LineWidth', 2, 'DisplayName',['\omega=',num2str(wnf(i))]);
    elseif (i ==4)
         plot(x,V(:,i), '--k', 'LineWidth', 2, 'DisplayName',['\omega=',num2str(wnf(i))]);
    else
         plot(x,V(:,i), 'DisplayName',['\omega=',num2str(wnf(i))]);
    end
         
         
    hold on;
end
grid on;
%title ('Relative Transverse displacement with shaft length','fontsize',16);
ylabel('Relative transverse displacement');
xlabel('Shaft length');
legend('show');
hold on
subplot(1,2,2)
for i = 1:nmode
    if (i == 1)
        plot(x,phi_x(:,i), '-k' , 'LineWidth', 2,'DisplayName',['\omega=',num2str(wnf(i))]);
    elseif (i == 2)
        plot(x,phi_x(:,i), ':k' , 'LineWidth', 2,'DisplayName',['\omega=',num2str(wnf(i))]);
    elseif (i == 3)
        plot(x,phi_x(:,i), '-.k' , 'LineWidth', 2,'DisplayName',['\omega=',num2str(wnf(i))]);
    elseif (i == 2)
        plot(x,phi_x(:,i), '--k' , 'LineWidth', 2,'DisplayName',['\omega=',num2str(wnf(i))]);
    else
        plot(x,phi_x(:,i), ':k' , 'LineWidth', 2,'DisplayName',['\omega=',num2str(wnf(i))]);
    end
        hold on;
end
grid on;
%title ('Relative Rotational displacement with shaft legth','fontsize',16);
ylabel('Relative rotational displacement');
xlabel('Shaft length');
legend('show');
hold on
saveas(h,'Ch8_Q14_FEM','png');

% print solution
fid = fopen('Ch8_Q14_FEM.txt','w');
fprintf(fid,'Transverse Vibration - FEM\n\n');
fprintf(fid,'Number of Elements = %d\n',nele);
fprintf(fid,'Density of Shaft = %d\n',rho);

if(nele <= 6)
for i=1:nele
    fprintf(fid,'Element-%d\n',i);
    fprintf(fid,'----------\n\n');
    fprintf(fid,'Mass Matrix [M]%d\n',i);
    for ii = 1:4
        for jj=1:4
            fprintf(fid,'%.2f\t',m(ii,jj,i));
        end
        fprintf(fid,'\n');
    end
    fprintf(fid,'\n');
    fprintf(fid,'Stiffness matrix [K]%d\n',i);
    for ii = 1:4
        for jj=1:4
            fprintf(fid,'%.2f\t',k(ii,jj,i));
        end
        fprintf(fid,'\n');
    end
    fprintf(fid,'\n');
end
fprintf(fid,'Global Mass matrix [M]\n');
for ii = 1:2*(nele+1)
    for jj=1:2*(nele+1)
        fprintf(fid,'%.2e\t',M(ii,jj));
    end
    fprintf(fid,'\n');
end
fprintf(fid,'\n');
fprintf(fid,'Global Stiffness matrix [K]\n');
for ii = 1:2*(nele+1)
    for jj=1:2*(nele+1)
        fprintf(fid,'%.2e\t',K(ii,jj));
    end
    fprintf(fid,'\n');
end
fprintf(fid,'\n');
fprintf(fid,'Global Mass matrix after applying boundary condition [M]\n');
for ii = 1:size(M_red,1)
    for jj=1:size(M_red,2)
        fprintf(fid,'%.2e\t',M_red(ii,jj));
    end
    fprintf(fid,'\n');
end
fprintf(fid,'\n');
fprintf(fid,'Global Stiffness matrix after applying boundary condition [K]\n');
for ii = 1:size(K_red,1)
    for jj=1:size(K_red,2)
        fprintf(fid,'%.2e\t',K_red(ii,jj));
    end
    fprintf(fid,'\n');
end
fprintf(fid,'\n');
fprintf(fid,'A = [K]^-1 [M]\n');
for ii = 1:size(K_red,1)
    for jj=1:size(K_red,2)
        fprintf(fid,'%.2e\t',D(ii,jj));
    end
    fprintf(fid,'\n');
end
fprintf(fid,'\n');
fprintf(fid,'Natural frequencies\n');
fprintf(fid,'%.2f\n',wnf);
fprintf(fid,'Relative transverse displacement at each node\n');
for ii= 1:2:size(mode,1)
    for jj = 1:length(wnf)
        fprintf(fid,'%.4f\t',mode(ii,jj));
    end
    fprintf(fid,'\n');
end
fprintf(fid,'Relative Rotational displacement at each node\n');
for ii= 2:2:size(mode,1)
    for jj = 1:length(wnf)
        fprintf(fid,'%.4f\t',mode(ii,jj));
    end
    fprintf(fid,'\n');
end
end
fclose(fid);