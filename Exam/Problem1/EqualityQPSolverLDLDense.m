function [x, lambda] = EqualityQPSolverLDLDense(H,g,A,b)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This functions solves equality constrained problems on the form 
%   min     f(x) = 1/2*x'Hx + g'x
%    x
%   s.t.    A'x = b
%  
%   It solves the system quite directly through the use of the
%   square-root-free cholesky factorization (LDL). This function assumes
%   the KKT matrix is represented densely.
%   With the output being a set of optimal values, x, and the
%   corresponding langrange multipliers, lambda.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Setup KKT
[n,m] = size(A);
KKT = [H, -A; -A', zeros(m)];
z = zeros(n+m, 1);
rhs = -[g;b];

% Conpute solution
[L,D,p] = ldl(KKT,'lower','vector');
z(p) = L'\(D\(L\rhs(p))); % Back-substitution matching the pivots

% Extract solution
x = z(1:n);
lambda = z(n+1:size(z));
