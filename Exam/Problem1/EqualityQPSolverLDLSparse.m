function [x, lambda] = EqualityQPSolverLDLSparse(H,g,A,b)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This functions solves equality constrained problems on the form 
%   min     f(x) = 1/2*x'Hx + g'x
%    x
%   s.t.    A'x = b
%  
%   It solves the system quite directly through the use of the
%   square-root-free cholesky factorization (LDL). This function assumes
%   the KKT matrix can be represented sparsely.
%   With the output being a set of optimal values, x, and the
%   corresponding langrange multipliers, lambda.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Setup the sparse KKT matrix
[n,m] = size(A);
KKT = [H, -A; -A', zeros(m)];
z = zeros(m + n, 1);
rhs = -[g;b];
KKT = sparse(KKT);

% Compute solution
[L,D,p] = ldl(KKT,'lower','vector');
z(p) = L'\(D\(L\rhs(p))); % Match the pivot elements

% Extract solution
x = z(1:n);
lambda = z(n+1:size(z));
