% Joseph Tate
% jnt664
% This code plots the location of n-bodies in an orbit around each other

clc; clear;

% input n bodies and organize data into array of cells
n = input('number of bodies: ');
traj = cell(n, 1);

for i = 1:n
    if i > 9
        traj{i} = load(sprintf('fort.10%d',i));
    else
        traj{i} = load(sprintf('fort.100%d',i));
    end
end

size = size(load('fort.1001'));
rows = size(1);

mass = [];
% read initial conditions file for mass data
startline1 = 9 - 1;
endline1 = 9 + n*6 - 2;
startline2 = endline1 + 2;
mass = dlmread('initialConditions.out', '\t',startline2, 0);

% compute center of mass
r = [];
M = sum(mass);
mx = 0;
my = 0;
mz = 0;

% x center of mass
for i = 1:n
    mx = mx + mass(i)* traj{i}(:,1);
end
r(:,1) = mx./M;

% y center of mass
for i = 1:n
    my = my + mass(i)* traj{i}(:,2);
end
r(:,2) = my./M;

% z center of mass
for i = 1:n
    mz = mz+ mass(i)* traj{i}(:,3);
end
r(:,3) = mz./M;

% calculate angular momentum sum
h = zeros(rows,3);
for i = 1:n
    for j = 1:rows
        h(j,:) = h(j,:) + mass(i)*cross(r(j,:),traj{i}(j,4:6));
    end
end

% calculate kinetic energy
T = zeros(rows,1);
for i = 1:n
    for j = 1:rows
        T(j) = T(j) + 1/2*mass(i)*(norm(traj{i}(j,4:6)))^(2);
    end
end

% calculate potential energy
U = zeros(rows,1);
for i = 1:n
    for j = 1:n
        if i == j
            continue
        else
            for k = 1:rows
                U(k) = U(k) + 1/2*mass(i)*mass(j)/(norm(traj{i}(k,1:3) - traj{j}(k,1:3)));
            end
        end
    end
end

% calculate energy total
c = T - U;

% plot center of mass components
time = traj{1}(:,7);

figure
hold on
plot(time,r(:,1), '-')
plot(time,r(:,2), '--')
plot(time,r(:,3), '-.')
title('Center of mass coordinates')
xlabel('time (TU)')
ylabel('distance (LU)')
legend('x component','y component','z component')
hold off

% plot angular momentum components
figure
hold on
plot(time,h(:,1), '-')
plot(time,h(:,2), '--')
plot(time,h(:,3), '-.')
title('Total angular momentum')
xlabel('time (TU)')
ylabel('angular momentum (MU*LU^2/TU)')

% plot and calculate angular momentum mean
hmean = zeros(rows, 1);
for i = 1:rows
    hmean(i) = norm(h(i,:));
end
plot(time,hmean(:,1), ':')
legend('x component','y component','z component','average')
hold off

% plot energy values
figure
hold on
plot(time,T(:,1), '-')
plot(time,U(:,1), '--')
plot(time,c(:,1), '-.')
title('Total energy')
xlabel('time (TU)')
ylabel('energy (MU*LU^2/TU^2)')
legend('Kinetic','Potential','Total (T - U)')
hold off

% plot trajectories and center of mass
figure
hold on
plot3(r(:,1),r(:,2),r(:,3), '-')
for i = 1:n
    plot3(traj{i}(:,1),traj{i}(:,2),traj{i}(:,3))
end
view(0,90)
axis equal
title('Body trajectories')
xlabel('x (LU)')
ylabel('y (LU)')
zlabel('z (LU)')
legendinfo = [];
for i = 0:n
    if i ==0
        legendinfo{i+1} = ['center of mass'];
    else
        a = ['body ' num2str(i)];
        legendinfo{i+1} = a;
    end
end
legend(legendinfo)
hold off

% % calculate and plot body 3 relative to 5
% r35 = traj{3}(:,1:3) - traj{5}(:,1:3);
% 
% figure
% hold on
% plot3(r35(:,1),r35(:,2),r35(:,3))
% axis equal
% view(0,90)
% title('Body 3 relative to body 5')
% xlabel('x (LU)')
% ylabel('y (LU)')
% zlabel('z (LU)')