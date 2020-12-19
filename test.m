m=1024;
n=1024;
A=randn(m,n);
b=randn(m,1);
mu=0.001;
y=ones(1024,1);
0.5*norm(A*x-b)^2+mu*norm(x,1)