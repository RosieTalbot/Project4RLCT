% Kuramoto Model
rng(1)

%% Variables

WVals = 1:0.1:3.2;
r = zeros(1,length(WVals));
for i = 1:length(WVals)
    W = WVals(i);
    r(i) = DoKuramoto(W);
end
% A = 3
% a = DoKuramoto(A)
figure(5)
plot(WVals,r,'LineWidth',4)
xlabel('W')
ylabel('Average r')
xticks(WVals)
title('Plot of r whilst varying K')
ylim([0,1])

function averager = DoKuramoto(W)
K = 1;
N = 500; % Number of particles
times1 = 0:0.01:100;
times2 =  0:0.01:1000;
omega = W+zeros(1,N);
%omega = ones(1,N)+pi; % Natural frequency of particles
th = 0:pi/50:2*pi;

% A is the matrix of the coupling values
A = rand(N)>0.1; A = A+A'; A(A ~=0) = 1; A = A/N;
A = A*K;
%A = 20*ones(N,N);

%% ODE

[~,Us] = ode45(@(t,U) RHS(U,A,omega,N), times1, mod(rand(1,N)*2*pi,2*pi));
[T,U] = ode45(@(t,U) RHS(U,A,omega,N), times2, Us(end,:));
M = mod(U,2*pi);

LastU = M(length(times2),:);
LastT = T(length(times2),:);

S = sin(LastU);
C = cos(LastU);

ordparam = 1/N*sum(exp(1i*M'));
r = abs(ordparam);
psi = atan2(imag(ordparam),real(ordparam));
averager = sum(r)/(length(times2));

% Mean = mean(LastU);  
% Stdev = std(LastU);
% Scale = (Stdev - 1.7919)/(-1.7919);

%% plotting
% figure(1);
% hold on;
% plot(T,M)

% figure(2);
% hold on;
% plot(C,S,'*')
% hold on
% plot(cos(th),sin(th))
% hold on
% plot(0,0,'.','MarkerSize',20)
% hold on
% quiver(0,0,cos(Mean),sin(Mean),Scale)
% axis square
% xlim([-1,1])
% ylim([-1,1])

% figure(3);
% hold on;
% plot(times,std(M'),'-*','LineWidth',3);
% axis tight;

figure(3);
hold on
plot(times2, r, '-', 'LineWidth', 4,'DisplayName',strcat('W = ',num2str(W)))
ylim([0,1])
xlabel('Time')
ylabel('Order Parameter Size')
title('Order Parameter Size over Time')
legend('show')

%% displaying stuff
disp(['W is: ',num2str(W)])
disp(['The initial r is: ',num2str(abs(1/N*sum(exp(i*M(1,:)))))])
disp(['The final r is: ',num2str(abs(1/N*sum(exp(i*M(length(times2),:)))))])

%% du = change in theta
function du = RHS(theta,A,omega,N)
du = zeros(N,1);
for i = 1:N
    %keyboard;
    du(i) = omega(i) + A(i,:)*sin(theta-theta(i));
end
end
end