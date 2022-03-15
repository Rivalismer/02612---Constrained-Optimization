function [x,lambda] = LUSolver(KKT,g,b,n)

z = zeros(2*n + 1, 1);
x = zeros(n+1,1);
lambda = zeros(n,1);
rhs = [g;b];

[L,U,p] = lu(KKT,'vector');

z(p) = U\(L\rhs(p));
x = -z(1:n+1);
lambda = z(n+2:size(z));