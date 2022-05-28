function [x, stats] = SQP_TrustRegion(objective, const, x0, lambda0)
%{
    This function is designed to solve nonlinear programming minimization 
    problems wrt. x vector in form of:

    min f(x)

    s.t      = gl <= g(x) <= gu
               xl <=  x   <= xu
    
    Output:
            -- x         - reached optimal value of x
            -- stats     - struct with iterations number, number of
                           function calls and iteration sequence
%}

% Solver settings and info
maxit = 100*length(x0);
tol   = 1.0e-5;
options = optimoptions('quadprog', 'Display', 'none');

% Initializing statistics         
stats.x(:, 1) = x0;

% Initializing variables
x = x0;
n = length(x);
lambda = lambda0;
B = eye(n);
it = 0;
mu = 20;     % penalty value
tr = 1;    % initial trust region 

% Evaluating function and constraints
[~, df] = objective(x);
[c, dc] = const(x); 
nc = length(c);
fcall = 1;

% Optimality conditions matrix
dL = df - dc*lambda;
F = [dL; c];

% Creating trust region B ang df matrices and constraints
Bhat = [B, zeros(n, nc); zeros(nc, n), zeros(nc)];
g = [df; mu*ones(nc, 1)];
dc = [dc', eye(nc)];

l = [-tr*ones(n, 1); ones(nc, 1)];
u = [tr*ones(n, 1); inf*ones(nc, 1)];
while ((it < maxit) && (norm(F(1:n),'inf') > tol))
    
    it = it + 1;
    
    % Solving quadratic sub-problem
    [sols, ~, ~, ~, lambda_qp] = quadprog(Bhat, g, -dc, c, [], [], l, u, [], options);

    % Iteration step 
    x = x + sols(1:n);
    lambda = lambda_qp.ineqlin;
    
    % Saving x
    stats.x(:, it+1) = x;
    
    % Evaluating function on updated parameters
    [~, df] = objective(x);
    [c, dc] = const(x);
    fcall = fcall + 1;
    
    % new Lagrangian gradient for BFGS check
    dLnew = df - dc*lambda;
    
    % NEW MATRICES FOR TRUST REGION
    g = [df; mu*ones(nc, 1)];
    dc = [dc', eye(nc)];
    
    % BFGS update
    p = x - stats.x(:, it);
    q = dLnew - dL;
    if p'*q >= 0.2*p'*B*p 
        theta = 1;
    else
        theta = (0.8*p'*B*p)/(p'*B*p-p'*q);
    end
    r = theta*q+(1-theta)*B*p;
    B = B + (r*r')/(p'*r) - (B*p)*(B*p)'/(p'*B*p);
    
    % Updating iterate
    Bhat = [B, zeros(n, nc); zeros(nc, n), zeros(nc)];
    dL = dLnew;
    F = [dL; c];
end
stats.iter = it;
stats.fcalls = fcall;
stats.lambda(:, 1) = lambda;