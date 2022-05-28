function [x, stats] = SQP_BFGS_LineSearch(objective, const, x0, lambda0)
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
options = optimoptions('quadprog','Algorithm','active-set','Display', 'none');
%options = optimoptions('quadprog', 'Display', 'none');
% For backtracking line search
rho = 0.25;              % has to be in (0,1)
tau = 0.8;               % has to be in (0,1)
eta = 0.1;               % has to be in (0,1), set to 0.1 as in lecture slides

% Initializing statistics         
stats.x(:, 1) = x0;

% Initializing variables
x = x0;
n = length(x);
lambda = lambda0;
B = eye(n);
it = 0;
p = ones(n,1);

% Evaluating function and constraints
[f, df] = objective(x);
[c, dc] = const(x); 
fcall = 1;

% Optimality conditions matrix
dL = df - dc*lambda;
F = [dL; c];

while ((it < maxit) && (norm(F(1:length(x)),'inf') > tol))
    
    it = it + 1;
    
    % Solving quadratic sub-problem
    [sols, ~, ~, ~, lambda_qp] = quadprog(B, df, -dc', c, [], [], [], [],p, options);
    
    % Backtracking line search
    plambda = lambda_qp.ineqlin - lambda;
    alpha = 1;
    pk = sols(1:n);
    mu = (df'*pk+0.5*pk'*B*pk)/((1-rho)*norm(c, 1));
    % Function initialization for conditions check
    f_ls = objective(x+alpha*pk); 
    fcall = fcall + 1;
    c_ls = const(x+alpha*pk);
    % Checking 
    while f_ls + mu*norm(c_ls, 1) > f + mu*norm(c,1) + eta*alpha*(df'*pk-mu*norm(c,1))
       % Updating alpha parameter
       alpha = tau*alpha; 
       f_ls = objective(x+alpha*pk);
       c_ls = const(x+alpha*pk);
    end
    
    % Updating iterates
    x = x + alpha*pk;
    lambda = lambda + alpha*plambda;
    
    % Saving x
    stats.x(:, it+1) = x;
    
    % Evaluating function on updated parameters
    [f, df] = objective(x);
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
    B = B + (r*r')/(p'*r) - (B*p)*(B*p)'/(p'*B*p);
    
    dL = dLnew;
    F = [dL; c];
end
stats.iter = it;
stats.fcalls = fcall;
stats.lambda(:, 1) = lambda;