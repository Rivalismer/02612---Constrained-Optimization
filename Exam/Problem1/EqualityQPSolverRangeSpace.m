function [x, lambda] = EqualityQPSolverRangeSpace(H,g,A,b)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This functions solves equality constrained problems on the form 
%   min     f(x) = 1/2*x'Hx + g'x
%    x
%   s.t.    A'x = b
%  
%   Specifically the RangeSpace solver should be employed when:
%   1. H is conditioned and easy to invert
%   2. H^-1 is known explicitly
%   3. The number of constraints is small
%   With the output being a set of optimal values, x, and the
%   corresponding langrange multipliers, lambda.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute the lower triangular matrix from cholesky factorization
[L, info] = chol(H,'lower');

% Check if matrix is positive definite
if info == 0
    v = L'\(L\g); % Use H = LL'
    Ha = A'*inv((L*L'))*A; % Use inverse without inverting H directly

    % Compute solutions
    lambda = Ha\(b + A'*v);
    x = (L*L')\(A*lambda - g);
    
    lambda = -lambda;
% If not then there is not solution with this solver
else
    x = [];
    lambda = [];
end