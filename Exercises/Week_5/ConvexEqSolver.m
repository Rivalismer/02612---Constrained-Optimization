function [x, lambda, f, c] = ConvexEqSolver(H,g,A,b)

% General solver for equality constrained quadratic programs
% Solves anything with the general form:
%       min     f(x) = 1/2 * x'Hx + g.*x
%       x
%       s.t.    c(x) = A'x + b = 0

%Get dimensions
[n,m] = size(A);

% Setup KKT
K = [H A; A' zeros(m,m)];
rhs = [g ; b];

% Factorize and solve with LDL
z = zeros(n+m,1);
[L,D,p] = ldl(K,'vector');
z(p) = L'\(D\(L\rhs(p)));

% Get solution
x = -z(1:n,1);
lambda = z((n+1):(n+m),1);

% Check for output arguments
if nargout > 2
    f = 0.5*x'*H*x + g'.*x;
end

if nargout > 3
    c = A'*x + b;
end
   