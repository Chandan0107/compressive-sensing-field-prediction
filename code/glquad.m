function [ y ] = glquad( fun, n)
%n-point Gauss Legendre quadrature rule for integration from 0 to 1
%std wgts,pts are for int from -1 to 1
%we change to int from 0 to 1, giving w' = w/2, x' = (x+1)/2
%In general it is better to use a rule like Clenshaw Curtis or Gauss
%Kronrad (matlab uses this)to estimate the accuracy of the quadrature
w = zeros(n,1); x = w;
if n == 1
    w(1) = 1; x(1) = 1;
elseif n == 2
    x(1) = (1-1.0/sqrt(3))/2; x(2) = (1+1.0/sqrt(3))/2;
    w(1) = 0.5; w(2) = 0.5;
end
y = 0;
for i=1:n
    y = y + w(i)*fun(x(i));
end
end
