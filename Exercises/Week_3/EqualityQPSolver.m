function [x,lambda] = EqualityQPSolver(H,g,A,b)

[m, n] = size(A);

x = zeros(n,1);
lambda = zeros(m,1);
temp = [H -A; -A' zeros(n)];
g_b_vec = -[g; b];

sol_vec = temp\g_b_vec;

x = sol_vec(1:n);
lambda = sol_vec(n:(n+m));