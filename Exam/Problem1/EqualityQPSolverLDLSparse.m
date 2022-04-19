function [x, lambda] = EqualityQPSolverLDLSparse(H,g,A,b)

[n,m] = size(A);
KKT = [H, -A; -A', zeros(m)];

z = zeros(m + n, 1);
rhs = -[g;b];
KKT = sparse(KKT);

[L,D,p] = ldl(KKT,'lower','vector');
z(p) = L'\(D\(L\rhs(p)));

x = z(1:n);
lambda = z(n+1:size(z));
