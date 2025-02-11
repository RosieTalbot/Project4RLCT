
rng(1);

N = 250;

omega = randn(1,N);

A = rand(N)>0.1; A = A+A'; A(A ~=0) = 1;

A = A*0.02;

[T,U] = ode15s(@(t,U)RHS(U,A,omega,N),0:0.01:500,randn(1,N));
plot(T,mod(U,2*pi))

function du = RHS(theta,A,omega,N)
du = zeros(N,1);
for i = 1:N
    du(i) = omega(i) + A(i,:)*sin(theta-theta(i));
end
end