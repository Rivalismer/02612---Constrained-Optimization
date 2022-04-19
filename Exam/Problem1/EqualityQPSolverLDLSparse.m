function [x, lambda] = EqualityQPSolverLDLSparse(H,g,A,b)

[n,m] = size(A);
KKT = [H, -A; A', zeros(m)];

z = zeros(2*n + 1, 1);
x = zeros(n+1,1);
lambda = zeros(n,1);
rhs = [-g;b];
KKT = sparse(KKT);

[L,D,p] = ldl(KKT,'lower','vector');
z(p) = L'\(D\(L\rhs(p)));

x = z(1:n+1);
lambda = z(n+2:size(z));
