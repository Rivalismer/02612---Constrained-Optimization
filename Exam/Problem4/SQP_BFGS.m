function [x, stats] = SQP_BFGS(objective, const, x0, lambda0)

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

% Evaluating function and constraints
[~, df] = objective(x);
[c, dc] = const(x); 
fcall = 1;
ncons = length(c);

% Optimality conditions matrix
dL = df - dc*lambda;
F = [dL; c];

while ((it < maxit) && (norm(F(1:length(x)),'inf') > tol))
    
    it = it + 1;
    
    % Solving quadratic sub-problem
    %[sols, ~, ~,  = [B -dc; dc' zeros(ncons)]\[-dL; -c];
    [sols, ~, ~, ~, lambda_qp] = quadprog(B, df, -dc', c, [], [], [], [], [], options);
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
    
    % BFGS update
    p = x - stats.x(:, it);
    q = dLnew - dL;
    if p'*q >= 0.2*p'*B*p 
        theta = 1;
    else
        theta = (0.8*p'*B*p)/(p'*B*p-p'*q);
    end
    r = theta*q+(1-theta)*B*p;
    B = B + (r*r')/(p'*r) - (B*p)*(B*p)'/(p'*B*p)
    
    dL = dLnew;
    F = [dL; c];
end
stats.iter = it;
stats.fcalls = fcall;
stats.lambda(:, 1) = lambda;