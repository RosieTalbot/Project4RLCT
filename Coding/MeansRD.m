
N = 50000;
omega = randn(1,N);
M = mod(omega,2*pi);
root = abs(omega);
mean(root)
std(root)

n = 1/N*exp(1i*M);
r = abs(n);

n = 1/sqrt(N);
n^2;
histogram(root,100)

