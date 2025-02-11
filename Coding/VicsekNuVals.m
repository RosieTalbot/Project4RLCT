%% Vicsek Nu Values

%% Loop

NuVals = -5:0.05:5;
OrdParam = zeros(1,length(NuVals));
for i = 1:length(NuVals)
    nu = NuVals(i);
    OrdParam(i) = DoVicsek(nu);
end

figure(5)
plot(NuVals,OrdParam,'LineWidth',4)
xlabel('Nu')
ylabel('Average W')
xticks(NuVals)
title('Plot of W whilst varying Nu')


function averageW = DoVicsek(nu)
%% Variables
N = 100; % number of particles
L = 20; % size of domain
v0 = 0.1; % speed of all particles
U = 0.5; % interacting radius
tMax = 1000;
tVals = 1:1:tMax;
dt = 1;
cat = zeros(2*N,length(tVals)); % to add future time
inta = -pi;
intb = pi;
psi = nu*(inta + (intb-inta).*rand(N,tMax)); % Noise
Z = zeros(N,length(tVals));
M = N*pi*U*U/(L*L)

% position of all particles at time 1 = r_t1 = x1,...,xN,y1,...,yN as a 2Nby1 matrix
% position of ith particle at time 1 = r_i_t1 = r([i,N+i],;)
% dirn/angle of velocities for particles = theta_i
% Circular Nbhds = U
% each centered at a x

%% Initial Conditions
% position of all particles 
r = L*rand(2*N,1);
r = horzcat(r,cat);
% angle of velocity for all particles
theta = 2*pi*rand(N,1);
% velocity of all particles
v = v0*[cos(theta); sin(theta)];
v = horzcat(v,cat);
%v = v0*exp(1i*theta); - which is best to use - matrix or scalar

%% Figure to test initial conditions
% figure(1)
% r = mod(r,L);
% % plot(r(1:N,i),r(N+1:2*N,i),'.b','MarkerSize',15)
% xlim([0 L]);
% ylim([0 L]);
% quiver(r(1:N,1),r(N+1:2*N,1),v(1:N,1),v(N+1:2*N,1),'MarkerSize',15)
% axis square
    

%% Time loop

for i = 2 : tMax
    
%% Calculating new thetas

    a = [r(1:N,i) r(N+1:2*N,i)];
    dist = pdist(a);
    D = squareform(dist);

    D1 = D<=U; % & D~=0;

    distsin = D1*sin(theta);
    distcos = D1*cos(theta);

    newtheta = atan2(distsin,distcos);

%% Updating values    
    % Angle and noise
    for j = 1:N
        if newtheta(j) ~= 0
            theta(j) = newtheta(j) + psi(j,i);
        else
            theta(j) = theta(j) + psi(j,i);
        end
        % Order parameter
        Z(j,i) = (1/N)*exp(1i*theta(j));
        % Velocity
        v(1:N,i) = v0*cos(theta);
        v(N+1:2*N,i) = v0*sin(theta);
        % Position
        r(1:N,i) = r(1:N,i-1) + dt*v(1:N,i);
        r(N+1:2*N,i) = r(N+1:2*N,i-1) + dt*v(N+1:2*N,i);
    end
    
    %% Order parameter
    Ztot = sum(Z,1);
    W = abs(Ztot);
    
    %% Plotting

    % figure(2)
    % r = mod(r,L);
    % plot(r(1:N,i),r(N+1:2*N,i),'.k','MarkerSize',12)
    % % colororder('glow12')
    % xlim([0 L]);
    % ylim([0 L]);
    % % quiver(r(1:N,i),r(N+1:2*N,i),v(1:N,i),v(N+1:2*N,i),'MarkerSize',15)
    % axis square
    % % hold on
    
end
averageW = sum(W)/(length(tVals));
figure(3)
plot(tVals,W,'.','MarkerSize',10)
ylim([0 1]);
hold on
end