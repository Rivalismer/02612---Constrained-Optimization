function [x,lambda] = EqualityQPSolverLUSparse(H,g,A,b)

% Insert some information about the QP solver
[n,m] = size(A);
KKT = [H, -A; -A', zeros(m)];

z = zeros(n + m, 1);
rhs = -[g;b];
KKT = sparse(KKT);

[L,U,p] = lu(KKT,'vector');

z = U\(L\rhs(p));
x = z(1:n);
lambda = z(n+1:size(z));