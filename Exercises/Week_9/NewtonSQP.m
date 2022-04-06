function [x,lambda] = NewtonSQP(fun, gfun, hfun, cfun, dcfun, hcfun, x0,...
    lambda0)

% Need to implement Newton's algorithm

g = feval(gfun,x0);
f = feval(fun,x0);
h = feval(hfun,x0);
c = feval(cfun,x0);
dc = feval(dcfun,x0);
hc = feval(hcfun,x0);
maxit = 100;
iter = 0;

dxk = -h\g;
xk = x0 + dxk;
lambda = lambda0;

options = optimoptions('quadprog', )
Converged = (norm(g - dc*lambda,'inf') > 1.e-8) && (norm(c,'inf') > 1.e-8);

while ~Converged && (iter < maxit)
    iter = iter + 1;
    H = h - sum(lambda*hc);
    [p, lambda] = quadprog(H,f,dc,-c);
    
    % Take step
    xk = xk + p;
    
    % Function Evaluation
    g = feval(gfun,x0);
    f = feval(fun,x0);
    h = feval(hfun,x0);
    c = feval(cfun,x0);
    dc = feval(dcfun,x0);
    hc = feval(hcfun,x0);
    
    % Convergence
    Converged = (norm(g - dc*lambda,'inf') > 1.e-8) && (norm(c,'inf') > 1.e-8);
end

x = xk;