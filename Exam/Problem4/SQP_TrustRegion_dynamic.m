function [x, stats] = SQP_TrustRegion_dynamic(objective, const, x0, lambda0)
% This function is not used in assignment sqp trust region algorithm with
% adaptive change of trust region in each iteration

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
tr = 10;    % initial trust region 
stats.tr(1) = tr;

% Evaluating function and constraints
[f, df] = objective(x);
[c, dc] = const(x); 
nc = length(c);
fcall = 1;

% Optimality conditions matrix
dL = df - dc*lambda;
F = [dL; c];

% Creating trust region B ang df matrices and constraints
Bhat = [B, zeros(n, nc); zeros(nc, n), zeros(nc)];
g = [df; mu*ones(nc, 1)];
dc = [dc', zeros(nc)];


while ((it < maxit) && (norm(F(1:n),'inf') > tol))
    
    it = it + 1;
    stats.tr(it+1) = tr;
    
    % Trust region constraints
    l = [-tr*ones(n, 1); ones(nc, 1)];
    u = [tr*ones(n, 1); inf*ones(nc, 1)];
    
    % Solving quadratic sub-problem
    [sols, ~, ~, ~, lambda_qp] = quadprog(Bhat, g, -dc, c, [], [], l, u, [], options);

    % Iteration step 
    x = x + sols(1:n);
    lambda = lambda_qp.ineqlin;
    
    % Saving x
    stats.x(:, it+1) = x;
    
    fold = f;
    dfold = df;
    % Evaluating function on updated parameters
    [f, df] = objective(x);
    [c, dc] = const(x);
    fcall = fcall + 1;
    
    % Update trust region
    p = sols(1:n);
    if p==zeros(n, 1)
        rho = 1
    else
        rho = (f-fold)/(0.5*p'*B*p+dfold'*p);
    end
    if rho <0.25
        corr = 0.25;
    else
        if rho <=0.75
            corr = 1;
        else
            corr = 2;
        end
    end
    if rho < 0 
        tr = corr*norm(p, inf);
    else
        tr = corr*tr;
    end
    
    
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