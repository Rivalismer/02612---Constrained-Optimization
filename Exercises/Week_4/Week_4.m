A = [4 -10];
b = 0;
x0 = [0; 0];
Aeq = [];
beq = [];
ub = [];
lb = [];
nonlcon = @circleon;

%% Standard

x = fmincon(@himmelblau, x0, A, b, Aeq, beq, ub, lb, nonlcon);

%% With Gradient

options = optimoptions('fmincon', 'SpecifyObjectiveGradient', true);

tic
[xg, fval, exitflag, output_g] = fmincon(@HimmelWithGrad, x0, A, b, Aeq, beq, ub, lb, nonlcon, options); 
toc

%% With Gradient and Hessian
%tolerance = 1.0e-6;
options = optimoptions('fmincon','Algorithm','sqp','SpecifyObjectiveGradient',true,'HessianFcn',@HimmelHessian);
tic
[xh, fval, exitflag, output_h] = fmincon(@HimmelWithGrad, x0, A, b, Aeq, beq, ub, lb, nonlcon, options); 
toc

%% Quadprog tests

H = [1 -1; -1 2];
f = [-2; -6];
A = [1 1; -1 2; 2 1];
b = [2; 2; 3];

[x,fval,exitflag,output,lambda] = quadprog(H,f,A,b);

%% Quadprog with options and printing

q_options = optimoptions('quadprog','Algorithm', 'interior-point-convex', 'Display','iter');
x0 = [];
Aeq = [];
beq = [];
ub = [];
lb = [];

[x,fval,exitflag,output,lambda] = quadprog(H,f,A,b, Aeq, beq, lb, ub, x0,q_options);

disp(lambda.ineqlin);
disp(lambda.eqlin);
disp(lambda.lower);
disp(lambda.upper);

%% Linprog test

% These are the constraints
A = [1 1
    1 1/4
    1 -1
    -1/4 -1
    -1 -1
    -1 1];

b = [2 1 2 1 -1 2];

% Objective function
f = [-1 -1/3];

x = linprog(f,A,b);

%% Linprog with options and printing

l_options = optimoptions('linprog', 'Algorithm', 'interior-point', 'Display', 'iter');

[x,fval,exitflag,output,lambda] = linprog(f, A, b, Aeq, beq, lb, ub, l_options);

disp(lambda.ineqlin);
disp(lambda.eqlin);
disp(lambda.lower);
disp(lambda.upper);