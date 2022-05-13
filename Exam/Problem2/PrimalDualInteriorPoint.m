function [x, iter, converged] = PrimalDualInteriorPoint(x0, y0, s0, z0, H, g, A, b, l, u)
%{
    This function is designed to solve quadratic programming minimization 
    problems wrt. x vector in form of:

    min f(x) = 1/2 x'Hx + g'x

    s.t      = A'x = b
               l <= x <= u
    
    Output:
            -- x         - reached optimal value of x
            -- iter      - number of iterations computed
            -- converged - 1 if algorithm did converge, 0 otherwise
%}

% Solver settings
max_iter = 100;         % Maximal iteration
tol = 1e-5;             % Allowed tolerance
eta = 0.995;             % Alpha correction parameter for iteration updates
% Not display linprog solver options
options = optimoptions('linprog','Display','off');   

% Setting initial number of iterations
iter = 0;

%% Initialization:

% Converting lower and upper boundaries to non-equality constraint in form Cx >= d

C = [diag(ones(1, length(x0))), diag(-ones(1, length(x0)))];
d = [l; -u];

% Preparing mc, S, Z, and e matrices and vectors
mc = length(s0);
S = diag(s0);
Z = diag(z0);
e = ones(mc, 1);

% Assigning variables from input
x = x0;
y = y0;
s = s0;
z = z0;

% Computing initial residuals
rL = H*x + g - A*y - C*z;
rA = b - A'*x;
rC = s + d - C'*x;
rSZ = S*Z*e;
mu = z'*s/mc;

%% Initial point

% LDL factorization of transformed KKT conditions and its preparation
[n, m] = size(A);
Hhat = H + C*(inv(S)*Z)*C';
KKT = [Hhat, -A; -A', zeros(m, m)];
[L, D, p] = ldl(KKT, 'lower', 'vector');

% Computation of initial affine direction
rLhat = rL - C*(inv(S)*Z)*(rC-inv(Z)*rSZ);
rs = -[rLhat; rA];
aff(p) = L'\(D\(L\rs(p)));
dxaff = (aff(1:length(x0)))';
dyaff = (aff(length(x0)+1:end))';
dzaff = -(inv(S)*Z)*C'*dxaff + (inv(S)*Z)*(rC - inv(Z)*rSZ);
dsaff = -inv(Z)*rSZ - inv(Z)*S*dzaff;

% Initial z and s vectors
z = max(1, abs(z + dzaff));
s = max(1, abs(s + dsaff));

%% Stopping criteria 

% Calculation of new S and Z matrices
S = diag(s);
Z = diag(z);

% Residuals estimation
rL = H*x + g - A*y - C*z;
rA = -A'*x + b;
rC = -C'*x + s + d;
rSZ = S*Z*e;
mu0 = mu;
mu = z'*s/mc;

% Checking if algorithm has already converged, without need of entering
% main loop, by comparing residuals to set tolerance
converged = (norm(rL, 'inf') <= tol * max(1, norm([H, g, A, C], 'inf'))) & ...
            (norm(rA, 'inf') <= tol * max(1, norm([A', b], 'inf'))) & ...
            (norm(rC, 'inf') <= tol * max(1, norm([eye(length(d)), d, C'], 'inf'))) & ...
            (mu <= tol * 10^(-2) * mu0);
         
%% Main loop

while ~converged && (iter < max_iter)
    
    % Iteration number
    iter = iter + 1;
    
    % Main loop LDL factorization of modified KKT conditions
    Hhat = H + C*(inv(S)*Z)*C';
    KKT = [Hhat, -A; -A', zeros(m, m)];
    [L, D, p] = ldl(KKT, 'lower', 'vector');
    
    % Main loop affine direction 
    rLhat = rL - C*(inv(S)*Z)*(rC-inv(Z)*rSZ);
    rs = -[rLhat; rA];
    aff(p) = L'\(D\(L\rs(p)));
    dxaff = (aff(1:length(x0)))';
    dyaff = (aff(length(x0)+1:end))';
    dzaff = -(inv(S)*Z)*C'*dxaff + (inv(S)*Z)*(rC - inv(Z)*rSZ);
    dsaff = -inv(Z)*rSZ - inv(Z)*S*dzaff;
    
    % Calculating largest alpha for duality gap and centering parameter 
    Aalph = -[dzaff; dsaff];
    balph = [z; s];
    alpha = linprog(-1, Aalph, balph, [], [], [], [], options);
    
    % Duality gap and centering parameter updates
    % Computation of the affine duality gap
    muaff = (z+alpha*dzaff)'*(s+alpha*dsaff)/mc;
    % Centering parameter
    sigma = (muaff/mu)^3;
    
    % Calculation of affine centering correction direction
    dSaff = diag(dsaff);
    dZaff = diag(dzaff);
    rSZhat = rSZ + dSaff*dZaff*e - sigma*mu*e;
    rLhat = rL - C*(inv(S)*Z)*(rC - inv(Z)*rSZhat);
    
    % Solving modified KKT conditions with respect to newly calcualated
    % rLhat value from above
    rs = -[rLhat; rA];
    sols(p) = L'\(D\(L\rs(p)));
    dx = (sols(1:length(x0)))';
    dy = (sols(length(x0)+1:end))';
    dz = -(inv(S)*Z)*C'*dx+(inv(S)*Z)*(rC-inv(Z)*rSZhat);
    ds = -inv(Z)*rSZhat - inv(Z)*S*dz;
    
    % Computing largest alpha parameter for duality gap calculation
    Aalph = -[dz; ds];
    balph = [z; s];
    alpha = linprog(-1, Aalph, balph, [], [], [], [], options);
    
    % Updating iteration parameters
    alphahat = eta*alpha;
    x = x + alphahat*dx;
    y = y + alphahat*dy;
    z = z + alphahat*dz;
    s = s + alphahat*ds;
    S = diag(s);
    Z = diag(z);
    
    % Residuals calculation for convergence check
    rL = H*x + g - A*y - C*z;
    rA = b - A'*x;
    rC = s + d - C'*x;
    rSZ = S*Z*e;
    mu = (z'*s)/mc;
    
    % Checking for convergence
    converged = (norm(rL, 'inf') <= tol * max(1, norm([H, g, A, C], 'inf'))) & ...
             (norm(rA, 'inf') <= tol * max(1, norm([A', b], 'inf'))) & ...
             (norm(rC, 'inf') <= tol * max(1, norm([eye(length(d)), d, C'], 'inf'))) & ...
             (mu <= tol * 10^(-2) * mu0);
end