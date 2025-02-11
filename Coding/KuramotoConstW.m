% Kuramoto Model

%% Variables

KVals = 0:0.1:10;

N = 100; % Number of particles
K = 5;
dt = 0.01;
times1 = 0:0.01:100;
times2 =  0:dt:1000;

% omega1 = zeros(1,N/2)+pi;% Natural frequency of particles
% omega2 = zeros(1,N/2)+1;
% omega = [omega1,omega2];

omega = rand(1,N);

% figure(8)
% histogram(omega,50)

% A is the matrix of the coupling values
A = rand(N)>0.1; A = A+A'; A(A ~=0) = 1; % A = A/N;
A = A*K;
%A = 20*ones(N,N);

%% ODE

[~,Us] = ode45(@(t,U) RHS(U,A,omega,N), times1, mod(2*pi*rand(1,N),2*pi));
[T,U] = ode45(@(t,U) RHS(U,A,omega,N), times2, Us(end,:));
M = mod(U,2*pi);

LastU = M(length(times2),:);
LastT = T(length(times2),:);

% Ushift = U;
% Ushift(1,:) = [];
% newcolm = zeros(1,1000);
% Ushifted = [Ushift;newcolm];
% diff = Ushifted - U;
% diff(100001,:) = [];

SMany = sin(U);
CMany = cos(U);

ordparam = 1/N*sum(exp(1i*M'));
firstordparam = 1/N*sum(exp(1i*M(1,:)));
lastordparam = 1/N*sum(exp(1i*M(length(times2),:)));

r = abs(ordparam);
psi = atan2(imag(ordparam),real(ordparam));

%% plotting
% figure(1);
% hold on;
% plot(T,M,'.','MarkerSize',5)
% xlabel('Time')
% ylabel('Omega')
% title('Phase against Time')
% fontsize(24,"points")


%% Moving circle plot
close all
for i = 1:length(times2)
    PlotCircle(ordparam(i),M(i,:))
    pause(0.01)
    hold off
end


% PlotR(times2,r,K)

% figure(2);
% hold on;
% plot(C,S,'*b','MarkerSize',10)
% hold on
% plot(cos(th),sin(th),'-r','LineWidth', 3)
% hold on
% quiver(0,0,x,y,0,'.m')
% hold on
% plot(0,0,'.k','MarkerSize',30)
% hold on
% plot(x,y,'.m','MarkerSize',40)
% axis square
% xlim([-1,1])
% ylim([-1,1])
% title('Angles of Particles at tMax')
% legend('Particles','Unit Circle','','Origin','OrderParameter')
% fontsize(24,"points")
% 
% figure(3);
% hold on
% plot(times2, r, '-', 'LineWidth', 4,'DisplayName',strcat('K = ',num2str(K)))
% ylim([0,1])
% xlabel('Time')
% ylabel('Order Parameter Size')
% title('Order Parameter Size over Time')
% legend('show')
% fontsize(24,"points")
% 
% figure(4);
% histogram(LastU,100)
% xlim([0,2*pi])
% hold on

%% displaying stuff
disp(['The initial r is: ',num2str(sqrt((real(firstordparam)).^2+(imag(firstordparam)).^2))])
disp(['The final r is: ', num2str(sqrt((real(lastordparam)).^2+(imag(lastordparam)).^2))])

%% du = change in theta
function du = RHS(theta,A,omega,N)
du = zeros(N,1);
for i = 1:N
    %keyboard;
    du(i) = omega(i) + A(i,:)*sin(theta-theta(i));
end

end