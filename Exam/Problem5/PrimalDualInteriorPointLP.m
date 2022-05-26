function [x_opt, iter, converged] = PrimalDualInteriorPointLP(x0, mu0, lambda0, g, A, b, l, u)
%{
    This function is designed to solve linear programming minimization 
    problems wrt. x vector in form of:

    min f(x) = g'x

    s.t      = A'x = b
               l <= x <= u

    Input notes:
            -- A matrix in input is in fact A' matrix, what is transformed
               at very beginning of function

    Output:
            -- x         - reached optimal value of x
            -- iter      - number of iterations computed
            -- converged - 1 if algorithm did converge, 0 otherwise
%}

% Solver settings
max_iter = 1000;         % Maximal iteration
tol = 1e-5;              % Allowed tolerance
eta = 0.995;             % Correction parameter for iteration updates
% Not display linprog solver options
options = optimoptions('linprog','Display','off');   

% Setting initial number of iterations
iter = 0;

%% Initialization:
    
% Converting input to standard form of Linear Problem
A = A';
[m, n] = size(A);
g = [g; zeros(2*length(x0), 1)];
b = [b; l; u];
x = [x0; x0-l; u-x0];
A = [A, zeros(m, n), zeros(m, n)
    eye(n), -eye(n), zeros(n, n) 
    eye(n), zeros(n, n), eye(n)];

mu = [mu0; ones(2*length(x0), 1)];
lambda = [lambda0; lambda0; lambda0];

% Preparing X, huge LAMBDA and e matrices and vectors
X = diag(x);
LAMBDA = diag(lambda);
e = ones(length(x), 1);
n = length(x);
lm = length(mu);

% Computing initial residuals
rL = g - A'*mu - lambda;
rA = A*x - b;
rXL = X*LAMBDA*e;
s = x'*lambda/n;

%% Stopping criteria 

% Checking if algorithm has already converged, without need of entering
% main loop, by comparing residuals to set tolerance
converged = (norm(rL, inf) <= tol) & ...
            (norm(rA, inf) <= tol) & ...
            (abs(s) <= tol);
         
%% Main loop

while ~converged && (iter < max_iter)
    
    % Iteration number
    iter = iter + 1;
    
    [m, n] = size(X);
    
    % Defining correctors step Newton's matrix
    Ls = [zeros(m, n), -A', -eye(m, n);
         A, zeros(lm, lm), zeros(lm, n);
         LAMBDA, zeros(n, lm), X];
    Rs = -[rL; rA; rXL];
    %Ls = sparse(Ls);
    
    % LU 
    [L,U,p] = lu(Ls,'vector');
    sols = U\(L\Rs(p));
  
    % Solving system of equations
    dxaff = (sols(1:n));
    dmuaff = (sols(n+1:n+lm));
    dlambdaaff = (sols(n+lm+1:end));
    
    % Calculating largest alpha and beta for duality gap and centering parameter 
    alpha = linprog(-1, -dxaff, x, [], [], [], [], options);
    beta = linprog(-1, -dlambdaaff, lambda, [], [], [], [], options);
    % Check whether alpha or beta isn't equal to infinity and prevent it
    if isempty(alpha) 
        alpha = 1;
    end
    if isempty(beta) 
        beta = 1;
    end
    
    % Updating affine iterations
    xaff = x + alpha*dxaff;
    lambdaaff = lambda + beta*dlambdaaff;
    saff = xaff'*lambdaaff/n;
    sigma = (saff/s)^3;
    s = x'*lambda/n;
    
    % Based on affine directions creating new diagonal matrices
    dXaff = diag(dxaff);
    dLAMBDAaff = diag(dlambdaaff);
    
    % Recalculating complementarity
    rXL = rXL + dXaff*dLAMBDAaff*e - sigma*s*e;
    Rs = -[rL; rA; rXL];
    
    % Solving system of equations
    sols = U\(L\Rs(p));
    dx = (sols(1:n));
    dmu = (sols(n+1:n+lm));
    dlambda = (sols(n+lm+1:end));
    
    % Calculating largest alpha and beta for duality gap and centering parameter 
    alpha = linprog(-1, -dx, x, [], [], [], [], options);
    beta = linprog(-1, -dlambda, lambda, [], [], [], [], options);
    % Check whether alpha or beta isn't equal to infinity and prevent it
    if isempty(alpha) 
        alpha = 1;
    end
    if isempty(beta) 
        beta = 1;
    end
    
    % Updating parameters
    x = x + eta*alpha*dx;
    mu = mu + eta*beta*dmu;
    lambda = lambda + eta*beta*dlambda;
    X = diag(x);
    LAMBDA = diag(lambda);
    
    % Residuals calculation for convergence check
    rL = g - A'*mu - lambda;
    rA = A*x - b;
    rXL = X*LAMBDA*e;
    s = x'*lambda/n;
    
    % Checking for convergence
    converged = (norm(rL, inf) <= tol) & ...
                (norm(rA, inf) <= tol) & ...
                (abs(s) <= tol);
            
end

% Extracting optimizer x from transformed problem
x_opt = x(1:length(x0));