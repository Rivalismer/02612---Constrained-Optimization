close all, clear all, clc
% This program is investigaring all implemented methods along with
% commercial solvers to solve 3 nonlinear optimization problems, each for 3
% different starting points

import casadi.*

%% HIMMELBLAU PROBLEM

% Initial points
x_init(1).x0 = [1; 2];
x_init(2).x0 = [-3; 0.75];
x_init(3).x0 = [5; 5];
lambda0 = [0.2; 0.2; 0.2; 0.2; 0.2; 0.2];

%################################
%           Damped BFGS
%################################
for i = [1, 2, 3]
    tic
    [x, stats] = SQP_BFGS(@himmelblau, @himmel_const, x_init(i).x0, lambda0); 
    cputime = toc;
    prob1.cpu(i, 1) = cputime;
    prob1.iter(i, 1) = stats.iter;
    prob1.fcalls(i, 1) = 3*stats.fcalls;
    prob1.xopt(i, 1:2) = x;
end

%################################
%           BFGS - line search
%################################
for i = [1, 2, 3]
    tic
    [x, stats] = SQP_BFGS_LineSearch(@himmelblau, @himmel_const, x_init(i).x0, lambda0); 
    cputime = toc;
    prob1.cpu(i, 2) = cputime;
    prob1.iter(i, 2) = stats.iter;
    prob1.fcalls(i, 2) = 4*stats.fcalls;
    prob1.xopt(i, 3:4) = x;
end

%################################
%           Trust Region
%################################
for i = [1, 2, 3]
    tic
    [x, stats] = SQP_TrustRegion(@himmelblau, @himmel_const, x_init(i).x0, lambda0); 
    cputime = toc;
    prob1.cpu(i, 3) = cputime;
    prob1.iter(i, 3) = stats.iter;
    prob1.fcalls(i, 3) = 3*stats.fcalls;
    prob1.xopt(i, 5:6) = x;
end

%################################
%           FMINCON
%################################
% Defining function for fmincon solver
himmelblau_fm = @(x)(x(1)^2 + x(2) - 11)^2 + (x(1) + x(2)^2 - 7)^2;
% Solving
for i=[1, 2, 3]
    tic
    [x, ~, ~, output] = fmincon(himmelblau_fm, x_init(i).x0, [4,-10], [0], [], [], [-5;-5], [5;5], @HimmCons);
    cputime = toc;
    prob1.cpu(i, 5) = cputime;
    prob1.iter(i, 5) = output.iterations;
    prob1.fcalls(i, 5) = output.funcCount;
    prob1.xopt(i, 9:10) = x;
end

%################################
%           CasADI
%################################
% Defining problem for CasADI
x1 = MX.sym('x1');
x2 = MX.sym('x2');
f = (x1^2 + x2 - 11)^2 + (x1 + x2^2 - 7)^2;
c1 = (x1+2)^2 - x2;
c2 = -4*x1 + 10*x2;
cons = [c1; c2];
nlp = struct;       % NLP declaration
nlp.x = [x1;x2];    % decision variables
nlp.f = f;          % objective
nlp.g = cons;       % constraints
F = nlpsol('F','ipopt',nlp);
for i = [1, 2, 3]
    cas_himmelblau = F('x0',x_init(i).x0,'ubg',1e8,'lbg',0,'lbx',[-5;-5],'ubx',[5;5]);
end
    
%% CUSTOM PROBLEM

% Initial points
x_init(1).x0 = [0.2; 0.2; 0.2];
x_init(2).x0 = [-0.7; -0.7; -0.7];
x_init(3).x0 = [1; -0.2; -0.2];
lambda0 = [0.2; 0.2; 0.2; 0.2; 0.2; 0.2; 0.2; 0.2];

%################################
%           Damped BFGS
%################################
for i = [1, 2, 3]
    tic
    [x, stats] = SQP_BFGS(@custom, @custom_const, x_init(i).x0, lambda0); 
    cputime = toc;
    prob2.cpu(i, 1) = cputime;
    prob2.iter(i, 1) = stats.iter;
    prob2.fcalls(i, 1) = 3*stats.fcalls;
    prob2.xopt(i, 1:3) = x;
end

%################################
%           BFGS - line search
%################################
for i = [1, 2, 3]
    try
        tic
        [x, stats] = SQP_BFGS_LineSearch(@custom, @custom_const, x_init(i).x0, lambda0); 
        cputime = toc;
        prob2.cpu(i, 2) = cputime;
        prob2.iter(i, 2) = stats.iter;
        prob2.fcalls(i, 2) = 4*stats.fcalls;
        prob2.xopt(i, 4:6) = x;
    catch
        prob2.cpu(i, 2) = NaN;
        prob2.iter(i, 2) = NaN;
        prob2.fcalls(i, 2) = NaN;
        prob2.xopt(i, 4:6) = [NaN, NaN, NaN];
    end
end

%################################
%           Trust region
%################################
for i = [1, 2, 3]
    tic
    [x, stats] = SQP_TrustRegion(@custom, @custom_const, x_init(i).x0, lambda0); 
    cputime = toc;
    prob2.cpu(i, 3) = cputime;
    prob2.iter(i, 3) = stats.iter;
    prob2.fcalls(i, 3) = 3*stats.fcalls;
    prob2.xopt(i, 7:9) = x;
end

%################################
%           FMINCON
%################################
% Defining function for fmincon solver
custom_fm = @(x)x(1)*x(2)+x(2)*x(3);
% Solving
for i=[1, 2, 3]
    tic
    [x, ~, ~, output] = fmincon(custom_fm, x_init(i).x0, [], [], [], [], [0;-1;0], [3;3;3], @CustomCons);
    cputime = toc;
    prob2.cpu(i, 5) = cputime;
    prob2.iter(i, 5) = output.iterations;
    prob2.fcalls(i, 5) = output.funcCount;
    prob2.xopt(i, 13:15) = x;
end

%################################
%           CasADI
%################################
% Defining problem for CasADI
x1 = MX.sym('x1');
x2 = MX.sym('x2');
x3 = MX.sym('x3');
f = x1*x2+x2*x3;
c1 = x1^2-x2^2+x3^2;
c2 = x1^2+x2^2+x3^2;
cons = [c1; c2];
nlp = struct;       % NLP declaration
nlp.x = [x1;x2;x3]; % decision variables
nlp.f = f;          % objective
nlp.g = cons;       % constraints
Cust = nlpsol('Cust','ipopt',nlp);
for i = [1, 2, 3]
    cas_custom = Cust('x0',x_init(i).x0,'ubg',[2;10],'lbg',[-10;-10],'lbx',[0;-1;0],'ubx',[3;3;3]);
end

%% CONSTRAINED ROSENBROCK

% Initial points
x_init(1).x0 = [0; 0];
x_init(2).x0 = [1; -0.8];
x_init(3).x0 = [-1; 2];
lambda0 = [0.2; 0.2; 0.2; 0.2; 0.2; 0.2; 0.2; 0.2];

%################################
%           Damped BFGS
%################################
for i = [1, 2, 3]
    tic
    [x, stats] = SQP_BFGS(@rosenbrock, @rosen_const, x_init(i).x0, lambda0); 
    cputime = toc;
    prob3.cpu(i, 1) = cputime;
    prob3.iter(i, 1) = stats.iter;
    prob3.fcalls(i, 1) = 3*stats.fcalls;
    prob3.xopt(i, 1:2) = x;
end

%################################
%           BFGS - line search
%################################
for i = [1, 2, 3]
    try
        tic
        [x, stats] = SQP_BFGS_LineSearch(@rosenbrock, @rosen_const, x_init(i).x0, lambda0); 
        cputime = toc;
        prob3.cpu(i, 2) = cputime;
        prob3.iter(i, 2) = stats.iter;
        prob3.fcalls(i, 2) = 4*stats.fcalls;
        prob3.xopt(i, 3:4) = x;
    catch
        prob3.cpu(i, 2) = NaN;
        prob3.iter(i, 2) = NaN;
        prob3.fcalls(i, 2) = NaN;
        prob3.xopt(i, 3:4) = [NaN, NaN];
    end
end

%################################
%           Trust region
%################################
for i = [1, 2, 3]
    try
        tic
        [x, stats] = SQP_TrustRegion(@rosenbrock, @rosen_const, x_init(i).x0, lambda0); 
        cputime = toc;
        prob3.cpu(i, 3) = cputime;
        prob3.iter(i, 3) = stats.iter;
        prob3.fcalls(i, 3) = 3*stats.fcalls;
        prob3.xopt(i, 5:6) = x;
    catch
        prob3.cpu(i, 3) = NaN;
        prob3.iter(i, 3) = NaN;
        prob3.fcalls(i, 3) = NaN;
        prob3.xopt(i, 5:6) = [NaN, NaN];
    end
end

%################################
%           FMINCON
%################################
% Defining function for fmincon solver
rosenbrock_fm = @(x)100*(x(2)-x(1)^2)^2 + (1-x(1))^2;
% Solving
for i=[1, 2, 3]
    tic
    [x, ~, ~, output] = fmincon(rosenbrock_fm, x_init(i).x0, [-4 10; 4 -10], [2; 1], [], [], [-1;-1], [2;2], @RosenCons); 
    cputime = toc;
    prob3.cpu(i, 5) = cputime;
    prob3.iter(i, 5) = output.iterations;
    prob3.fcalls(i, 5) = output.funcCount;
    prob3.xopt(i, 9:10) = x;
end

%################################
%           CasADI
%################################
% Defining problem for CasADI
x1 = MX.sym('x1');
x2 = MX.sym('x2');
f = (1-x1)^2+100*(x2-x1^2)^2;
c1 = (x1+2)^2 - x2;
c2 = -4*x1 + 10*x2;
cons = [c1; c2];
nlp = struct;       % NLP declaration
nlp.x = [x1;x2];    % decision variables
nlp.f = f;          % objective
nlp.g = cons;       % constraints
Rosen = nlpsol('Rosen','ipopt',nlp);
for i = [1, 2, 3]
    cas_rosenbrock = Rosen('x0',x_init(i).x0,'ubg',2,'lbg',-1,'lbx',[-1;-1],'ubx',[2;2]);
end

% Writing results for CasADI
prob1.cpu(:, 4) = [0.029; 0.017; 0.006];
prob1.iter(:, 4) = [10; 9; 10];
prob1.fcalls(:, 4) = [54; 56; 53];
prob1.xopt(:, 7:8) = prob1.xopt(:, 9:10);

prob2.cpu(:, 4) = [0.012; 0.004; 0.001];
prob2.iter(:, 4) = [8; 8; 8];
prob2.fcalls(:, 4) = [44; 44; 44];
prob2.xopt(:, 10:12) = prob2.xopt(:, 13:15);

prob3.cpu(:, 4) = [0.012; 0.019; 0.018];
prob3.iter(:, 4) = [9; 10; 11];
prob3.fcalls(:, 4) = [49; 44; 49];
prob3.xopt(:, 7:8) = prob3.xopt(:, 9:10);