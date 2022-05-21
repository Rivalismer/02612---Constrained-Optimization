function [x,lambda] = EqualityQPSolverLUSparse(H,g,A,b)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This functions solves equality constrained problems on the form 
%   min     f(x) = 1/2*x'Hx + g'x
%    x
%   s.t.    A'x = b
%  
%   It solves the system quite directly through the use of the
%   lower-upper factorization (LU). This function assumes
%   the KKT matrix can be represented sparsely.
%   With the output being a set of optimal values, x, and the
%   corresponding langrange multipliers, lambda.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up the sparse KKT system
[n,m] = size(A);
KKT = [H, -A; -A', zeros(m)];
rhs = -[g;b];
KKT = sparse(KKT);

% Compute solution
[L,U,p] = lu(KKT,'vector');
z = U\(L\rhs(p)); % Back substitution matching pivot on rhs

% Extract solution
x = z(1:n);
lambda = z(n+1:size(z));