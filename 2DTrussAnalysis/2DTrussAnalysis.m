%This program is developed for the analysis of 2D truss structures given
%the number of nodes, their locations in the x-y plane, their corresponding
%element number, cross-sectional area, Young's Modulus, applied forces, and
%known displacements
%COE 321K Spring 2018
%Joseph Tate
%**************************************************************************
clc;  clear all; hold off;

%Read data from input file
finput = fopen('project_input.txt');

nnode = fscanf(finput, ['nodes: %d\n %*s %*s\n']);        % # of nodes
size_nc = [2 nnode];
node_coor = fscanf(finput,'%f',size_nc);          % coordinates of each node
node_coor = transpose(node_coor);
node_coor = 12.*node_coor;                         % convert to inches

nelem = fscanf(finput, ['\n elements: %d\n %*s %*s %*s %*s\n']);% # of elements
size_ne = [4 nelem];
elemdata = fscanf(finput,'%f',size_ne); % node1, node2, Area, Young's Modulus
elemdata = transpose(elemdata);

nforce = fscanf(finput, ['\n force_BCs: %d\n %*s %*s %*s\n']);% # of force BCs
size_nf = [3 nforce];
forcedata = fscanf(finput,'%f',size_nf);    % node, degree of freeedom, value
forcedata = transpose(forcedata);

ndisp = fscanf(finput, ['\n displacement_BCs: %d\n %*s %*s %*s\n']);% number of displacement BCs
size_nf = [3 ndisp];
dispdata = fscanf(finput,'%f',size_nf); % node, degree of freedom, value
dispdata = transpose(dispdata);

fclose(finput);                                      % close the input file

% %**************************************************************************

% Create global stiffness matrix K and force vector P
k = [];
connec = [];
for n = 1:nelem
    %build k matrix
    node1 = elemdata(n,1);
    node2 = elemdata(n,2);
    node1x = node_coor(node1, 1);
    node1y = node_coor(node1, 2);
    node2x = node_coor(node2, 1);
    node2y = node_coor(node2, 2);
    L = sqrt((node2x - node1x)^2 + (node2y - node1y)^2);
    c = (node2x - node1x)/L;
    s = (node2y - node1y)/L;
    k(:,:,n) = elemdata(n,4).*elemdata(n,3)./L.*[c^2 c*s -c^2 -c*s; 
                                                 c*s s^2 -c*s -s^2;
                                                -c^2 -c*s c^2 c*s; 
                                                -c*s -s^2 c*s s^2];
    % build connectivity vector
    connec(:,:,n) = [2*(node1 - 1) + 1
                     2*(node1 - 1) + 2
                     2*(node2 - 1) + 1
                     2*(node2 - 1) + 2];
end

% Build K matrix
K = zeros(nnode*2);
for n = 1:nelem
    for i = 1:4
        iglobal = connec(i,1,n);
        for j = 1:4
            jglobal = connec(j,1,n);
            K(iglobal,jglobal) = K(iglobal,jglobal) + k(i,j,n);
        end
    end
end

% build P vector
P = zeros(nnode*2,1);             %zeros accounts for Penalty Method BCs
for n = 1:nforce
    node = forcedata(n,1);
    dof = forcedata(n,2);
    nglobal = 2*(node - 1) + dof;
    P(nglobal) = forcedata(n,3);
end

%**************************************************************************

%Modify K to account for boundary conditions
M = 100000000*elemdata(1,4);               % set number of penalty method
for n = 1:ndisp
    node = dispdata(n,1);
    dof = dispdata(n,2);
    nglobal = 2*(node - 1) + dof;
    K(nglobal,nglobal) = K(nglobal,nglobal) + M;
end

%**************************************************************************

%Solve for displacements   
U=K\P;

%**************************************************************************
%Solve for reactions Ps and internal forces F
Ps = zeros(17,2);
for n = 1:ndisp
    node = dispdata(n,1);
    dof = dispdata(n,2);
    nglobal = 2*(node - 1) + dof;
    if dof == 2
        Ps(nglobal/2,dof) = -M*U(nglobal);
    else
        Ps(nglobal/2 + 0.5,dof) = -M*U(nglobal);
    end
end

F = zeros(nelem,1);
Sigma = zeros(nelem,1);
Sigmax = 0;
Sigmin = 0;
Sigmaxelem = 1;
Sigminelem = 1;
for n = 1:nelem
    E = elemdata(n,4);
    A = elemdata(n,3);
    node1 = elemdata(n,1);
    node2 = elemdata(n,2);
    nglobal1x = node1*2 - 1;
    nglobal2x = node2*2 - 1;
    nglobal1y = node1*2;
    nglobal2y = node2*2;
    node1x = node_coor(node1, 1);
    node1y = node_coor(node1, 2);
    node2x = node_coor(node2, 1);
    node2y = node_coor(node2, 2);
    L = sqrt((node2x - node1x)^2 + (node2y - node1y)^2);
    c = (node2x - node1x)/L;
    s = (node2y - node1y)/L;
    e = [c s];
    ua = [U(nglobal1x) U(nglobal1y)];
    ub = [U(nglobal2x) U(nglobal2y)];
    u = ub - ua;
    F(n) = E*A/L*dot(u,e);
    Sigma(n) = F(n)/A;
    if Sigma(n) > Sigmax
        Sigmax = Sigma(n);
        Sigmaxelem = n;
    elseif Sigma(n) < Sigmin
            Sigmin = Sigma(n);
            Sigminelem = n;
    end
end



%**************************************************************************

%Give output 
disp('     DISPLACEMENT RESULTS (inches)       ')
disp('     Node        x-dir(u)        y-dir(v)')
disp('                                         ')
Up=zeros(nnode,3);
for i=1:nnode
   Up(i,1)=i;
   Up(i,2)=U(2*i-1);
   Up(i,3)=U(2*i);
end
fprintf('     %i %18.3E %15.3E\n',transpose(Up))
%***********************
disp('                                         ')
disp('                                         ')
disp('     REACTION RESULTS (lbs)              ')
disp('     Node        x-dir(u)        y-dir(v)')
disp('                                         ')
Pp=zeros(nnode,3);
for i=1:nnode
   Pp(i,1)=i;
   Pp(i,2)=Ps(i,1);
   Pp(i,3)=Ps(i,2);
end
fprintf('     %i %18.3f %15.3f\n',transpose(Pp))
%***********************
disp('                                         ')
disp('                                         ')
disp('     MEMBER FORCES AND STRESSES          ')
disp('     Elem.       Force(lbs)   Stress(psi)')
disp('                                         ')
Fp=zeros(nelem,3);
for i=1:nelem
   Fp(i,1)=i;
   Fp(i,2)=F(i);
   Fp(i,3)=Sigma(i);
end
fprintf('     %i %18.3f %15.3f\n',transpose(Fp))
disp('                                         ')
disp('                                         ')

%*************************************************************************

%Show deformed structure (magnified k times)
k = 10000;
limits=[min(node_coor(:,1))-1, max(node_coor(:,1))+1, min(node_coor(:,2))-1, max(node_coor(:,2))+1];
Coord1 = zeros(length(elemdata),4);
for i=1:length(elemdata)
   Coord1(i,1:2) = node_coor(elemdata(i,1),1:2);
   Coord1(i,3:4) = node_coor(elemdata(i,2),1:2);
end
X1 = [Coord1(:,1) Coord1(:,3)];
Y1 = [Coord1(:,2) Coord1(:,4)];

node_coor = node_coor+k*Up(:,2:3)/12;
for i=1:length(elemdata)
   Coord2(i,1:2) = node_coor(elemdata(i,1),1:2);
   Coord2(i,3:4) = node_coor(elemdata(i,2),1:2);
end
X2 = [Coord2(:,1) Coord2(:,3)];
Y2 = [Coord2(:,2) Coord2(:,4)];

plot(X1(1,:),Y1(1,:),'b','LineWidth',5),hold on
axis(limits)
plot(X2(1,:),Y2(1,:),'r--','LineWidth',2)
for i=2:length(Coord1)
    plot(X1(i,:),Y1(i,:),'b','LineWidth',5)
    plot(X2(i,:),Y2(i,:),'r--','LineWidth',2)
end
stitle = sprintf('Undeformed & Deformed Structure (magnification factor: %d)',k);
title(stitle)



