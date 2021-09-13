clear;
clc;
     
% Giiven data    
E = 2.1e+11;
d =0.01;        %shaft dia in m
rho = 0; % 0 for no shaft mass        %shaft density in kg/m^3
L = 1; % shaft length in m

count =0;
nele = 2; % 4 or 6

% Connetivity matrix
connect = [1 1 2;
           2 2 3];  % Second & Third Column are Nodes (in sequence)

 % x coordinates of nodes
 coord = [1 0;       % first column node second column x coodinate
         2 0.3;
         3 1];
     
 % Boundary conditions    
 bc_data = [2, 1;     % First Column is Node number
            2, 3;
            3, 1;    % Second column is prescribed local dof
            3, 2;
            3, 3;
            3, 4;];   % dof 1 for u;2 for phiy; 3 for v; 4 for phix

 mass = [1 5 0.02 0.001 (pi/6) 0.02; % lumped masses in kg and kg-m^2
        2 0 0 0 0 0;
        3 0 0 0 0 0];  %1st column node, 2nd colum mass, 3rd column Id
                     % 4th column eccentricity, 5th column phase
                     % sixth column trial mass in kg (if any)
                   
% Mass and stifness matrices and force vector
NDOFnode=4; %size(coord, 1); DOFs per node
M = zeros(NDOFnode*(nele+1));
K = zeros(NDOFnode*(nele+1));

   
for i= 1: nele
    nd1 = connect(i,2);
    nd2 = connect(i,3);
    x1  = coord(nd1,2);
    x2  = coord(nd2,2);
    l = x2-x1;
    A = pi*d^2/4;
    rho*A*l/420;
    % Mass matrix in z-x plane
    mzx(:,:,i) = rho*A*l/420*[156     22*l    54      -13*l;
                            22*l    4*l^2   13*l    -3*l^2;
                            54      13*l    156     -22*l;
                            -13*l   -3*l^2  -22*l   4*l^2];
    fx(:,i) = [0 0 0 0];
    % Mass matrix in y-z plane
    myz(:,:,i) = rho*A*l/420*[156     -22*l    54    13*l;
                            -22*l    4*l^2   -13*l    -3*l^2;
                            54      -13*l    156     22*l;
                            13*l   -3*l^2  22*l   4*l^2];
                   
    if size(mass)>0                    
    if any(mass(:,1)==nd1)      %checking for disc at left node of element
        ii = find(mass(:,1)==nd1);
        mzx(1,1,i) = mzx(1,1,i) + mass(ii,2);
        myz(1,1,i) = myz(1,1,i) + mass(ii,2);
        mzx(2,2,i) = mzx(2,2,i) + mass(ii,3);
        % diametral mass moment of inertia
        myz(2,2,i) = myz(2,2,i) + mass(ii,3);
       
    end
   
     if i==nele                  % checking for disc at right side of last element
        if any(mass(:,1)==nd2)
            ii = find(mass(:,1)==nd2);
            mzx(3,3,i) = mzx(3,3,i) + mass(ii,2);
            myz(3,3,i) = myz(3,3,i) + mass(ii,2);
            mzx(4,4,i) = mzx(4,4,i) + mass(ii,3);
            myz(4,4,i) = myz(4,4,i) + mass(ii,3);
           
        end
    end
   
    end
       
    I =pi*d^4/64;
    E*I/l^3;
   % Stiffness matrix in z-x plane
    kzx(:,:,i) = E*I/l^3 *[12     6*l     -12     6*l;
                         6*l    4*l^2   -6*l    2*l^2;
                         -12    -6*l    12      -6*l;
                         6*l    2*l^2   -6*l    4*l^2];
 
    % Stiffness matrix in y-z plane
    kyz(:,:,i) = E*I/l^3 *[12    -6*l     -12    -6*l;
                         -6*l    4*l^2   6*l    2*l^2;
                         -12    6*l    12      6*l;
                         -6*l    2*l^2   6*l    4*l^2];
                     
                     
    %Assembly (order is u1, phiy1 v1 phix1, u2,phiy2, v2 phix2 ...
   
    vec1 = [1+(nd1-1)*4,2+(nd1-1)*4,1+(nd2-1)*4,2+(nd2-1)*4];
    M(vec1,vec1) = M(vec1,vec1)+ mzx(:,:,i);
    K(vec1,vec1) = K(vec1,vec1)+ kzx(:,:,i);
   
    vec2 = [3+(nd1-1)*4,4+(nd1-1)*4,3+(nd2-1)*4,4+(nd2-1)*4];
    M(vec2,vec2) = M(vec2,vec2)+ myz(:,:,i);
    K(vec2,vec2) = K(vec2,vec2)+ kyz(:,:,i);
 
end


% Unbalance force vector
w=500;
dw=0.01;

for jj=1:w/dw % rad/s
   
   f = zeros(size(M(:,1))); % unbalance force vector

    if (jj ==1)
         omega(jj) = dw;
    else
        omega(jj) = omega(jj-1) + dw;
    end
    for i= 1: nele
    nd1 = connect(i,2);
    nd2 = connect(i,3);
   
    % force vector in z-x plane
    fx(:,i) = [0 0 0 0];
   
    % force vector in y-z plane
    fy(:,i) = [0 0 0 0];
                   
    if size(mass)>0                    
    if any(mass(:,1)==nd1)      %checking for disc at left node of element
        ii = find(mass(:,1)==nd1);
        fx(1,i) = fx(1,i) + mass(ii,6)*mass(ii, 4)*omega(jj)^2*exp(-j*mass(ii,5)); % unbalance force mew^2exp(-jtheta)
        fy(1,i) = fy(1,i) -j*fx(1,i);
        %moment
        %fx(2,i) = fx(2,i) + mass(ii,2)*mass(ii, 4)*omega(jj)^2*exp(-j*mass(ii,5));
        % diametral mass moment of inertia
        %fy(2,i) = fy(2,i) -j*fx(2,i);
    end
   
     if i==nele                  % checking for disc at right side of last element
        if any(mass(:,1)==nd2)
            ii = find(mass(:,1)==nd2);
            fx(3,i) = fx(3,i) + mass(ii,6)*mass(ii, 4)*omega(jj)^2*exp(-j*mass(ii,5));
            fy(3,i) = fy(3,i) -j*fx(3,i);
            %fx(4,i) = fx(4,i) + mass(ii,2)*mass(ii, 4)*omega(jj)^2*exp(-j*mass(ii,5));
            %fy(4,i) = fy(4,i) -j*fx(4,i);
        end
    end
   
    end
   
    % imposition of boundary condition
    supp_dof = [];      % dof to be suppressed.
for iii = 1:size(bc_data,1)
    nd = bc_data(iii,1);
    dof = bc_data(iii,2);
    Gdof = 4*(nd-1)+dof;        %finding the global DOF of given BC
    supp_dof = [supp_dof , Gdof];
   
end

supp_dof = sort(supp_dof);
                 
                     
    %Assembly (order is u1, phiy1 v1 phix1, u2,phiy2, v2 phix2 ...
   
    vec1 = [1+(nd1-1)*4,2+(nd1-1)*4,1+(nd2-1)*4,2+(nd2-1)*4];
    f(vec1) = f(vec1) + fx(:,i);
   
    vec2 = [3+(nd1-1)*4,4+(nd1-1)*4,3+(nd2-1)*4,4+(nd2-1)*4];
    f(vec2) = f(vec2) + fy(:,i);
   
end

% Reducing Global matrix (Imposition of Boundary condition)
% ---------------------------------------------------------
K_red = K;
M_red = M;
f_red = f; % force vector
for i = 1:size(supp_dof,2)
    dof = supp_dof(i);  %Getting the Global DOF of Given BC
    if (dof == 1)       %Loop for reducing global matrix
        K_red = K_red(dof+1:end, dof+1:end);
        M_red = M_red(dof+1:end, dof+1:end);
        f_red = f_red(dof+1:end);
    elseif (dof == 4*(nele+1))
            K_red = K_red(1:dof-1, 1:dof-1);
            M_red = M_red(1:dof-1, 1:dof-1);
            f_red = f_red(1:dof-1);
    else
        K_red = K_red([1:dof-1, dof+1:end],[1:dof-1, dof+1:end]);
        M_red = M_red([1:dof-1, dof+1:end],[1:dof-1, dof+1:end]);
        f_red = f_red([1:dof-1, dof+1:end]);
    end
    if (i~=size(supp_dof,2))    %for considring the already reduced steps
        supp_dof(i+1:end) = supp_dof(i+1:end)-1;
    end
end

% Forced vibration analysis (without damping)
X(:,jj)= (K_red - omega(jj)^2 * M_red)\f_red;

% With Rayleigh damping from Example 9.9
%X(:,jj)= (K_red - omega(jj)^2 * M_red + j*omega(jj)*(0.1897*M_red+0.0005*K_red))\f_red;

end

figure (2)
subplot(2,1,1)
semilogy(omega, abs(X(3,:)), 'k', 'LineWidth', 2)
xlabel('Frequency (rad/s)')
ylabel('Linear displ amplitude x1 (m)')
subplot(2,1,2)
plot(omega, angle(X(3,:)), 'k', 'LineWidth', 2)
xlabel('Frequency (rad/s)')
ylabel('Linear displ phase x1 (rad)')
%hold

figure (3)
subplot(2,1,1)
semilogy(omega, abs(X(6,:)), 'k', 'LineWidth', 2)
xlabel('Frequency (rad/s)')
ylabel('Amplitude (m)')
xlabel('Angular displ amplitude (rad/s)')
ylabel('Angular displ amplitude  phi1 (rad)')
subplot(2,1,2)
plot(omega, angle(X(3,:)), 'k', 'LineWidth', 2)
xlabel('Frequency (rad/s)')
ylabel('Angular displ phase phi1 (rad)')



% Free vibration analysis

D=M_red/K_red;
[eig_vec,eig_val] = eig(D);
mode = [];
wnf=[];
for i=1:1:size(eig_val,2)
    if eig_val(i,i)~=0
        wnf = [wnf,sqrt(1/eig_val(i,i))];
        mode = [mode, eig_vec(:,i)];
    end
end

for i =1:2      %arranging in asceending order
    ii = [1:length(wnf)/2]+(i-1)*length(wnf)/2;
    for j =ii
        for jj = j:ii(end)
            if wnf(j)>wnf(jj)
                temp_wnf = wnf(j);
                wnf(j) = wnf(jj);
                wnf(jj) = temp_wnf;
                temp_mode(:,1) = mode(:,j);
                mode(:,j) = mode(:,jj);
                mode(:,jj) = temp_mode(:,1);
            end
        end
    end
    wnf_sort(i:2:length(wnf)) = wnf(ii);
    mode_sort(:,i:2:length(wnf)) = mode(:,ii);
end
wnf = wnf_sort;
mode = mode_sort;
unique(wnf)

% updating boundary condition

for i = 1:size(bc_data,1)
    nd = bc_data(i,1);
    dof = bc_data(i,2);
    Gdof = 4*(nd-1)+dof;    %Finding global DOF
    if Gdof ==1
        mode = [zeros(1,size(mode,2));mode];
    else
        mode = [mode((1:Gdof-1),:);zeros(1,size(mode,2));mode((Gdof:end),:)];
    end
end

% Post processing


for i =1:2:length(wnf)
    x = [];
    v=[];
    phi=[];
    for ii= 1:nele
        nd1 = connect(ii,2);
        nd2 = connect(ii,3);
        x1  = coord(nd1,2);
        x2  = coord(nd2,2);
        l = x2-x1;
        x_n =x1:(x2-x1)/50:x2; % to get smoother mode shape curve
        v_n=[];
        phi_n=[];
        eta = mode([1+(nd1-1)*4,2+(nd1-1)*4,1+(nd2-1)*4,2+(nd2-1)*4],i);
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
    [~,jj]=max(abs(v));
    if(v(jj)~=0)
    v  =normalize(v,'scale',v(jj));
    end
    [~,jj]=max(abs(phi));
    if phi(jj~=0)
    phi  =normalize(phi,'scale',phi(jj));
    end
    V(:,i) = v';
    phi_x(:,i)  = phi';
end

% ploting mode shape

h = figure(1);
set(gcf, 'Position', get(0,'Screensize'));
subplot(1,2,1)

if(size(wnf,2) < 8)
    nmode=size(wnf,2);
else
    nmode=8;
end
 


for i =1:2:nmode %length(wnf)
   % plot(x,V(:,i),'DisplayName',['\omega=',num2str(wnf(i))]);
       if(i ==1)
        plot(x,V(:,i), '-k', 'LineWidth', 2, 'DisplayName',['\omega=',num2str(wnf(i))]);
    elseif (i ==3)
         plot(x,V(:,i), ':k', 'LineWidth', 2, 'DisplayName',['\omega=',num2str(wnf(i))]);
    elseif (i ==5)
         plot(x,V(:,i),'-.k', 'LineWidth', 2, 'DisplayName',['\omega=',num2str(wnf(i))]);
    elseif (i ==7)
         plot(x,V(:,i), '--k', 'LineWidth', 2, 'DisplayName',['\omega=',num2str(wnf(i))]);
    else
         %plot(x,V(:,i), 'DisplayName',['\omega=',num2str(wnf(i))]);
    end

    hold on;
end
grid on;
title ('Relative Transverse displacement with shaft legth','fontsize',16);
ylabel('Relative Transverse displacement');
xlabel('Shaft length');
legend('show');
hold on
subplot(1,2,2)

for i = 1:2:nmode %length(wnf)
    %plot(x,phi_x(:,i),'DisplayName',['\omega=',num2str(wnf(i))]);
    if (i == 1)
        plot(x,phi_x(:,i), '-k' , 'LineWidth', 2,'DisplayName',['\omega=',num2str(wnf(i))]);
    elseif (i == 3)
        plot(x,phi_x(:,i), ':k' , 'LineWidth', 2,'DisplayName',['\omega=',num2str(wnf(i))]);
    elseif (i == 5)
        plot(x,phi_x(:,i), '-.k' , 'LineWidth', 2,'DisplayName',['\omega=',num2str(wnf(i))]);
    elseif (i == 7)
        plot(x,phi_x(:,i), '--k' , 'LineWidth', 2,'DisplayName',['\omega=',num2str(wnf(i))]);
    else
        plot(x,phi_x(:,i), ':k' , 'LineWidth', 2,'DisplayName',['\omega=',num2str(wnf(i))]);
    end

    hold on;
end
grid on;
title ('Relative Rotational displacement with shaft legth','fontsize',16);
ylabel('Relative Rotational displacement');
xlabel('Shaft length');
legend('show');
hold on
saveas(h,'Ch9_Q1_FEM','png');

% print solution
fid = fopen('Ch9_Q1_FEM','w');
fprintf(fid,'Transverse Vibration Forced Vibration- FEM\n\n');
fprintf(fid,'Number of Elements = %d\n',nele);
fprintf(fid,'Density of shaft = %d\n',rho);
for i=1:nele
    fprintf(fid,'Element-%d\n',i);
    fprintf(fid,'----------\n\n');
   
    fprintf(fid,'Mass Matrix [Mzx]%d\n',i);
    for ii = 1:4
        for jj=1:4
            fprintf(fid,'%.2f\t',mzx(ii,jj,i));
        end
        fprintf(fid,'\n');
    end

    fprintf(fid,'Mass Matrix [Myz]%d\n',i);
    for ii = 1:4
        for jj=1:4
            fprintf(fid,'%.2f\t',myz(ii,jj,i));
        end
        fprintf(fid,'\n');
    end

    fprintf(fid,'\n');
   
    fprintf(fid,'Stiffness matrix [Kzx]%d\n',i);
    for ii = 1:4
        for jj=1:4
            fprintf(fid,'%.2f\t',kzx(ii,jj,i));
        end
        fprintf(fid,'\n');
    end
   
    fprintf(fid,'Stiffness matrix [Kyz]%d\n',i);
    for ii = 1:4
        for jj=1:4
            fprintf(fid,'%.2f\t',kyz(ii,jj,i));
        end
        fprintf(fid,'\n');
    end
   
    fprintf(fid,'\n');
   
   
    fprintf(fid,'\n');
   
    fprintf(fid,'Unbalance force vector [fx]%d\n',i);
    %for ii = 1:4
        for jj=1:4
            fprintf(fid,'%.2f\t',fx(jj));
        end
        fprintf(fid,'\n');
    %end
   
    fprintf(fid,'Unbalance force vector [fy]%d\n',i);
    %for ii = 1:4
        for jj=1:4
            fprintf(fid,'%.2f\t',fy(jj));
        end
        fprintf(fid,'\n');
    %end
   
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


fprintf(fid,'Global force vector{f}\n');
%for ii = 1:2*(nele+1)
    for jj=1:2*(nele+1)
        fprintf(fid,'%.2e\t',f(jj));
    end
    fprintf(fid,'\n');
%end
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

fprintf(fid,'Global force vector after applying boundary condition {f_red}\n');
%for ii = 1:size(K_red,1)
    for jj=1:size(f_red,1)
        fprintf(fid,'%.2e\t',f_red(jj));
    end
    fprintf(fid,'\n');
%end
fprintf(fid,'\n');

fprintf(fid,'For eidenvaue solution A = [K](^-1) [M]\n');
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

fclose(fid);