function [x,lambda] = EqualityQPSolverNullSpace(H,g,A,b)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This functions solves equality constrained problems on the form 
%   min     f(x) = 1/2*x'Hx + g'x
%    x
%   s.t.    A'x = b
%  
%   The null-space solver does not require non-singularity in H, only that
%   the constraint matrix A has full row rank. However, the price for this
%   is that it can be more costly than other solvers when H is
%   non-singular.
%   The null-space solver is useful when the degrees of freedom is a small
%   number
%   The output is a set of optimal values, x, and the
%   corresponding langrange multipliers, lambda.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[n,~] = size(A);

% Set up the system with a QR factorization
[Q,Rbar] = qr(A);
m1 = size(Rbar,2);
Q1 = Q(:,1:m1);
Q2 = Q(:,m1+1:n);
R = Rbar(1:m1,1:m1);

% Compute x (indirectly with x_y and x_z) and lambda
x_y = R'\b;
x_z = (Q2'*H*Q2)\(-Q2'*(H*Q1*x_y + g));
x = (Q1*x_y + Q2*x_z);
lambda = (R\Q1'*(H*x+g));

% Invert solution
x = -x;
lambda = -lambda;