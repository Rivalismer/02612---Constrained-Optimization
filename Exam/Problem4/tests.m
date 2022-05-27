clear all, clc, close all
objective = @himmelblau;
const = @himmel_const;
x0 = [1; 3];
lambda0 = [0.2; 0.2];

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

while ((it < maxit) && (norm(F,'inf') > tol))
    
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
end
stats.iter = it;
stats.fcalls = fcall;

%% see
objfun = @himmelblau;
const = @himmel_const;
x0 = [1; 3];
lambda0 = [0.2; 0.2];

[~,g] = objfun(x);
[c,A] = const(x);

dxL = g-A'*lambda; % gradient of the Lagrangian  
F = [dxL; c];

k = 1;
n = length(x);
x(:,k) = x;

B = eye(n); % Hessian approximation initial value

while ((k < maxit) && (norm(F,'inf') > tol))
    
    % solve the system to get the search direction
    z = [B -A'; A zeros(2)]\[-dxL; -c];
    % obtain next point in the iteration sequence
    p = z(1:n);
    x(:,k+1) = x(:,k) + p;
    % update lambda value
    lambda = lambda + z(n+1:end);
    k = k+1;
    % evaluate function and constraints on the new point
    [~,g] = objfun(x(:,k));
    [c,A] = const(x(:,k));
    dxLn = g-A'*lambda; % new gradient of the Lagrangian
    % check KKT conditions
    F = [dxLn; c];
    % damped BFGS updating
    s = x(:,k) - x(:,k-1);
    y = dxLn - dxL;
    dxL = dxLn;
     % ensure possitive definiteness
    if s'*y >= 0.2*s'*B*s 
        theta = 1;
    else
        theta = (0.8*s'*B*s)/(s'*B*s-s'*y);
    end
    r = theta*y+(1-theta)*B*s;
    % update Bs
    B = B - (B*(s*s')*B)/(s'*B*s) + r*r'/(s'*r); 
end
% store number of iterations and iteration sequence
info.iter = k; 
info.seq = x;
x = x(:,k); % return last point as solution