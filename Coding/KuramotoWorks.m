% Kuramoto Model

%% Variables

KVals = 0:0.1:10;

N = 500; % Number of particles
K = 1.7;
times1 = 0:0.01:100;
times2 =  0:0.01:1000;

omega = randn(1,N); % Natural frequency of particles
th = 0:pi/50:2*pi;

% A is the matrix of the coupling values
A = rand(N)>0.1; A = A+A'; A(A ~=0) = 1; A = A/N;
A = A*K;
%A = 20*ones(N,N);

%% ODE

[T,Us] = ode45(@(t,U) RHS(U,A,omega,N), times1, mod(pi+randn(1,N),2*pi));
[T,U] = ode45(@(t,U) RHS(U,A,omega,N), times2, Us(end,:));
M = mod(U,2*pi);

LastU = M(length(times2),:);
LastT = T(length(times2),:);

S = sin(LastU);
C = cos(LastU);

ordparam = 1/N*sum(exp(1i*M'));
firstordparam = 1/N*sum(exp(1i*M(1,:)));
lastordparam = 1/N*sum(exp(1i*M(length(times2),:)));

% r = sqrt((real(ordparam)).^2+(imag(ordparam)).^2);

r = abs(ordparam);
psi = atan2(imag(ordparam),real(ordparam));

x = real(lastordparam);
y = imag(lastordparam);

%% plotting
figure(1);
hold on;
plot(T,M)
xlabel('Time')
ylabel('Omega')
title('Phase against Time')
fontsize(24,"points")

figure(3);
hold on
plot(times2, r, '-', 'LineWidth', 4)
ylim([0,1])
xlabel('Time')
ylabel('Order Parameter Size')
title('Order Parameter Size over Time')
legend(num2str(K))
fontsize(24,"points")

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













