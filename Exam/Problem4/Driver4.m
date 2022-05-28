%% Question 4
clear all; close all; clc;
%{
    This program is creating contour plot of Himmelblau problem, with 
    saddle points and minimia location
%}
% Creating grid
x=-5:0.05:5;
[X,Y]=meshgrid(x);
Z=(X.^2+Y-11).^2+(X+Y.^2-7).^2;
v=[0:2:10 10:10:100 100:20:200];

% Creating contour plot
[c, h]=contour(X, Y, Z, v, 'linewidth', 2);
colorbar, hold on;
yc1 = (x+2).^2;
yc2 = (4*x)/10;
hold on
    title('Himmelblau problem')
    fill(x,yc1,[0.7 0.7 0.7],'facealpha',0.4);
    fill([x x(end) x(1)],[yc2 -5 -5],[0.7 0.7 0.7],'facealpha',0.4);
    xlim([-5,5]);
    ylim([-5,5]);
    scatter([-3.65, -3.55, -0.3, 3], [2.72, -1.42, 2.89, 2], 50, 'filled', 'red');
    scatter([-3.1, -1.4, -0.55, 0.0867, 3.25], [-0.1, 0.36, -0.22, 2.8843, 1.3], 50, 'filled', 'black');

%% Question 5
clear all; close all; clc;
%{
    This program is solving three defined nonlinear programming problems
    with use of commercial symbolic solvers - CadADI and fmincon
%}
% Using CasADI solver


%% Himmelblau problem - CasADI
import casadi.*
% Symbols/expressions
x1 = MX.sym('x1');
x2 = MX.sym('x2');

% Objective function
f = (x1^2 + x2 - 11)^2 + (x1 + x2^2 - 7)^2;

% Constraints
c1 = (x1+2)^2 - x2;
c2 = -4*x1 + 10*x2;
cons = [c1; c2];

% Model declaration
nlp = struct;       % NLP declaration
nlp.x = [x1;x2];    % decision variables
nlp.f = f;          % objective
nlp.g = cons;       % constraints

% Create solver instance
F = nlpsol('F','ipopt',nlp);

% Solve the problem using a guess
cas_himmelblau = F('x0',[0.0 0.0],'ubg',1e8,'lbg',0,'lbx',[-5;-5],'ubx',[5;5]);

cas_himmelblau.x

%% Custom problem - CasADI
import casadi.*
% Symbols/expressions
x1 = MX.sym('x1');
x2 = MX.sym('x2');
x3 = MX.sym('x3');

% Objective function
f = x1*x2+x2*x3;

% Constraints
c1 = x1^2-x2^2+x3^2;
c2 = x1^2+x2^2+x3^2;
cons = [c1; c2];

% Model declaration
nlp = struct;       % NLP declaration
nlp.x = [x1;x2;x3]; % decision variables
nlp.f = f;          % objective
nlp.g = cons;       % constraints

% Create solver instance
Cust = nlpsol('Cust','ipopt',nlp);

% Solve the problem using a guess
cas_custom = Cust('x0',[0.0 0.0 0.0],'ubg',[2;10],'lbg',[-10;-10],'lbx',[0;-1;0],'ubx',[3;3;3]);

cas_custom.x

%% Rosenbrock problem - CasADI
import casadi.*
% Symbols/expressions
x1 = MX.sym('x1');
x2 = MX.sym('x2');

% Objective function
f = (1-x1)^2+100*(x2-x1^2)^2;

% Constraints
c1 = (x1+2)^2 - x2;
c2 = -4*x1 + 10*x2;
cons = [c1; c2];

% Model declaration
nlp = struct;       % NLP declaration
nlp.x = [x1;x2];    % decision variables
nlp.f = f;          % objective
nlp.g = cons;       % constraints

% Create solver instance
Rosen = nlpsol('Rosen','ipopt',nlp);

% Solve the problem using a guess
cas_rosenbrock = Rosen('x0',[0.0 0.0],'ubg',2,'lbg',-1,'lbx',[-1;-1],'ubx',[2;2]);

cas_rosenbrock.x

%% Himmelblau problem - fmincon

% Defining function for fmincon solver
himmelblau = @(x)(x(1)^2 + x(2) - 11)^2 + (x(1) + x(2)^2 - 7)^2;

% Solving 
x_himm = fmincon(himmelblau, [0; 0], [4,-10], [0], [], [], [-5;-5], [5;5], @HimmCons);   

%% Custom problem - fmincon

% Defining function for fmincon solver
custom = @(x)x(1)*x(2)+x(2)*x(3);

% Solving 
x_cust = fmincon(custom, [0;0;0], [], [], [], [], [0;-1;0], [3;3;3], @CustomCons);


%% Rosenbrock problem - fmincon

rosenbrock = @(x)100*(x(2)-x(1)^2)^2 + (1-x(1))^2;

x_rose = fmincon(rosenbrock, [0;0], [-4 10; 4 -10], [2; 1], [], [], [-1;-1], [2;2], @RosenCons);

%% Printing results
disp('Himmelblau problem:')
disp('x_opt')
disp('   CasADI')
disp(cas_himmelblau.x)
disp('   fmincon')
disp(x_himm')

disp('Custom problem:')
disp('x_opt')
disp('   CasADI')
disp(cas_custom.x)
disp('   fmincon')
disp(x_cust')

disp('Rosenbrock problem:')
disp('x_opt')
disp('   CasADI')
disp(cas_rosenbrock.x)
disp('   fmincon')
disp(x_rose')

%% Question 6
clear all; close all; clc;
% This program is investigating optimization of Himmelblau problem with use 
% of damped BFGS Lagrangian approximation method.

% initial guesses
x01 = [1;2];
x02 = [-3;0.75];
x03 = [5;5];
lambda0 = 0.2*ones(6, 1);

% Calculating minimizers
[x_s1, stats_s1] = SQP_BFGS(@himmelblau, @himmel_const, x01, lambda0);
[x_s2, stats_s2] = SQP_BFGS(@himmelblau, @himmel_const, x02, lambda0);
[x_s3, stats_s3] = SQP_BFGS(@himmelblau, @himmel_const, x03, lambda0);

x=-5:0.05:5;
[X,Y]=meshgrid(x);
Z=(X.^2+Y-11).^2+(X+Y.^2-7).^2;
v=[0:2:10 10:10:100 100:20:200];

% Plotting iteration sequence
[c, h]=contour(X, Y, Z, v, 'linewidth', 2);
colorbar, hold on;
yc1 = (x+2).^2;
yc2 = (4*x)/10;
    title('Himmelblau problem')
    fill(x,yc1,[0.7 0.7 0.7],'facealpha',0.4);
    fill([x x(end) x(1)],[yc2 -5 -5],[0.7 0.7 0.7],'facealpha',0.4);
    xlim([-5,5]);
    ylim([-5,5]);
    pl1 = scatter(stats_s1.x(1, :), stats_s1.x(2, :), 'r', 'filled');
    pl2 = scatter(stats_s2.x(1, :), stats_s2.x(2, :), 'm', 'filled');
    pl3 = scatter(stats_s3.x(1, :), stats_s3.x(2, :), 'b', 'filled');
    legend({'', '', '', 'x_0 = [1  2]''', 'x_0 = [-3  0.75]''', 'x_0 = [5  5]'''}, 'location', 'south');
hold off
 
%% Question 7
clear all; close all; clc;
% This program is investigating optimization of Himmelblau problem with use 
% of damped BFGS Lagrangian with backtracking line search approximation method.

% initial guesses
x01 = [1;2];
x02 = [-3;0.75];
x03 = [5;5];
lambda0 = 0.2*ones(6, 1);

% Calculating minimizers
[x_s1, stats_s1] = SQP_BFGS_LineSearch(@himmelblau, @himmel_const, x01, lambda0);
[x_s2, stats_s2] = SQP_BFGS_LineSearch(@himmelblau, @himmel_const, x02, lambda0);
[x_s3, stats_s3] = SQP_BFGS_LineSearch(@himmelblau, @himmel_const, x03, lambda0);

x=-5:0.05:5;
[X,Y]=meshgrid(x);
Z=(X.^2+Y-11).^2+(X+Y.^2-7).^2;
v=[0:2:10 10:10:100 100:20:200];

% Plotting iteration sequence
[c, h]=contour(X, Y, Z, v, 'linewidth', 2);
colorbar, hold on;
yc1 = (x+2).^2;
yc2 = (4*x)/10;
hold on
    title('Himmelblau problem')
    fill(x,yc1,[0.7 0.7 0.7],'facealpha',0.4);
    fill([x x(end) x(1)],[yc2 -5 -5],[0.7 0.7 0.7],'facealpha',0.4);
    xlim([-5,5]);
    ylim([-5,5]);
    pl1 = plot(stats_s1.x(1, :), stats_s1.x(2, :), 'r', 'Linewidth', 2);
    pl2 = plot(stats_s2.x(1, :), stats_s2.x(2, :), 'm', 'Linewidth', 2);
    pl3 = plot(stats_s3.x(1, :), stats_s3.x(2, :), 'b', 'Linewidth', 2);
    legend({'', '', '', 'x_0 = [1  2]''', 'x_0 = [-3  0.75]''', 'x_0 = [5  5]'''}, 'location', 'south');
hold off

%% Question 8
clear all, close all, clc
% This program is investigating optimization of Himmelblau problem with use 
% of trust region with BFGS hessian aproximation method.

% initial guesses
x01 = [1;2];
x02 = [-3;0.75];
x03 = [5;5];
lambda0 = 0.2*ones(6, 1);

% Calculating minimizers
[x_s1, stats_s1] = SQP_TrustRegion(@himmelblau, @himmel_const, x01, lambda0);
[x_s2, stats_s2] = SQP_TrustRegion(@himmelblau, @himmel_const, x02, lambda0);
[x_s3, stats_s3] = SQP_TrustRegion(@himmelblau, @himmel_const, x03, lambda0);

x=-5:0.05:5;
[X,Y]=meshgrid(x);
Z=(X.^2+Y-11).^2+(X+Y.^2-7).^2;
v=[0:2:10 10:10:100 100:20:200];

% Plotting iteration sequence
[c, h]=contour(X, Y, Z, v, 'linewidth', 2);
colorbar, hold on;
yc1 = (x+2).^2;
yc2 = (4*x)/10;
hold on
title('Himmelblau problem')
fill(x,yc1,[0.7 0.7 0.7],'facealpha',0.4);
fill([x x(end) x(1)],[yc2 -5 -5],[0.7 0.7 0.7],'facealpha',0.4);
xlim([-5,5]);
ylim([-5,5]);
pl1 = plot(stats_s1.x(1, :), stats_s1.x(2, :), 'r', 'Linewidth', 2);
pl2 = plot(stats_s2.x(1, :), stats_s2.x(2, :), 'm', 'Linewidth', 2);
pl3 = plot(stats_s3.x(1, :), stats_s3.x(2, :), 'b', 'Linewidth', 2);
legend({'', '', '', 'x_0 = [1  2]''', 'x_0 = [-3  0.75]''', 'x_0 = [5  5]'''}, 'location', 'south');

%% Question 9
close all, clear all, clc
% This program is investigaring all implemented methods along with
% commercial solvers to solve 3 nonlinear optimization problems, each for 3
% different starting points


%% HIMMELBLAU PROBLEM
import casadi.*

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
import casadi.*

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
import casadi.*

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

%% Plotting - Iterations
figure
subplot(1,3,1)
hold on;
scatter(1, prob1.iter(1, 1), 30,  'red', 'filled');
scatter(1, prob1.iter(1, 2), 30, 'blue', 'filled');
scatter(1, prob1.iter(1, 3), 30, 'cyan', 'filled');
scatter(1, prob1.iter(1, 4), 30, 'yellow', 'filled');
scatter(1, prob1.iter(1, 5), 30, 'black', 'filled');
xlabel('Set of initial values')
ylabel('Iterations')
set(gca, 'FontSize', 12)
hold off
legend({'BFGS', 'Line Search', 'Trust Region', 'CasADI', 'fmincon'}, 'location', 'northeast')

subplot(1,3,2)
hold on
scatter(2, prob1.iter(2, 1), 30,  'red', 'filled');
scatter(2, prob1.iter(2, 2), 30, 'blue', 'filled');
scatter(2, prob1.iter(2, 3), 30, 'cyan', 'filled');
scatter(2, prob1.iter(2, 4), 30, 'yellow', 'filled');
scatter(2, prob1.iter(2, 5), 30, 'black', 'filled');
title('Number of iterations for the methods - Problem 1')
xlabel('Set of initial values')
ylabel('Iterations')
set(gca, 'FontSize', 12)
hold off

subplot(1,3,3)
hold on
scatter(3, prob1.iter(3, 1), 30,  'red', 'filled');
scatter(3, prob1.iter(3, 2), 30, 'blue', 'filled');
scatter(3, prob1.iter(3, 3), 30, 'cyan', 'filled');
scatter(3, prob1.iter(3, 4), 30, 'yellow', 'filled');
scatter(3, prob1.iter(3, 5), 30, 'black', 'filled');
xlabel('Set of initial values')
ylabel('Iterations')
set(gca, 'FontSize', 12)
hold off

figure
subplot(1,3,1)
hold on;
scatter(1, prob2.iter(1, 1), 30,  'red', 'filled');
scatter(1, prob2.iter(1, 2), 30, 'blue', 'filled');
scatter(1, prob2.iter(1, 3), 30, 'cyan', 'filled');
scatter(1, prob2.iter(1, 4), 30, 'yellow', 'filled');
scatter(1, prob2.iter(1, 5), 30, 'black', 'filled');
xlabel('Set of initial values')
ylabel('Iterations')
set(gca, 'FontSize', 12)
hold off
legend({'BFGS', 'Line Search', 'Trust Region', 'CasADI', 'fmincon'}, 'location', 'northeast')

subplot(1,3,2)
hold on
scatter(2, prob2.iter(2, 1), 30,  'red', 'filled');
scatter(2, prob2.iter(2, 2), 30, 'blue', 'filled');
scatter(2, prob2.iter(2, 3), 30, 'cyan', 'filled');
scatter(2, prob2.iter(2, 4), 30, 'yellow', 'filled');
scatter(2, prob2.iter(2, 5), 30, 'black', 'filled');
title('Number of iterations for the methods - Problem 2')
xlabel('Set of initial values')
ylabel('Iterations')
set(gca, 'FontSize', 12)
hold off

subplot(1,3,3)
hold on
scatter(3, prob2.iter(3, 1), 30,  'red', 'filled');
scatter(3, prob2.iter(3, 2), 30, 'blue', 'filled');
scatter(3, prob2.iter(3, 3), 30, 'cyan', 'filled');
scatter(3, prob2.iter(3, 4), 30, 'yellow', 'filled');
scatter(3, prob2.iter(3, 5), 30, 'black', 'filled');
xlabel('Set of initial values')
ylabel('Iterations')
set(gca, 'FontSize', 12)
hold off

figure
subplot(1,3,1)
hold on;
scatter(1, prob3.iter(1, 1), 30,  'red', 'filled');
scatter(1, prob3.iter(1, 2), 30, 'blue', 'filled');
scatter(1, prob3.iter(1, 3), 30, 'cyan', 'filled');
scatter(1, prob3.iter(1, 4), 30, 'yellow', 'filled');
scatter(1, prob3.iter(1, 5), 30, 'black', 'filled');
xlabel('Set of initial values')
ylabel('Iterations')
set(gca, 'FontSize', 12)
hold off
legend({'BFGS', 'Line Search', 'Trust Region', 'CasADI', 'fmincon'}, 'location', 'northeast')

subplot(1,3,2)
hold on
scatter(2, prob3.iter(2, 1), 30,  'red', 'filled');
scatter(2, prob3.iter(2, 2), 30, 'blue', 'filled');
scatter(2, prob3.iter(2, 3), 30, 'cyan', 'filled');
scatter(2, prob3.iter(2, 4), 30, 'yellow', 'filled');
scatter(2, prob3.iter(2, 5), 30, 'black', 'filled');
title('Number of iterations for the methods - Problem 3')
xlabel('Set of initial values')
ylabel('Iterations')
set(gca, 'FontSize', 12)
hold off

subplot(1,3,3)
hold on
scatter(3, prob3.iter(3, 1), 30,  'red', 'filled');
scatter(3, prob3.iter(3, 2), 30, 'blue', 'filled');
scatter(3, prob3.iter(3, 3), 30, 'cyan', 'filled');
scatter(3, prob3.iter(3, 4), 30, 'yellow', 'filled');
scatter(3, prob3.iter(3, 5), 30, 'black', 'filled');
xlabel('Set of initial values')
ylabel('Iterations')
set(gca, 'FontSize', 12)
hold off
%% Plotting - fcalls
figure
subplot(1,3,1)
hold on;
scatter(1, prob1.fcalls(1, 1), 30,  'red', 'filled');
scatter(1, prob1.fcalls(1, 2), 30, 'blue', 'filled');
scatter(1, prob1.fcalls(1, 3), 30, 'cyan', 'filled');
scatter(1, prob1.fcalls(1, 4), 30, 'yellow', 'filled');
scatter(1, prob1.fcalls(1, 5), 30, 'black', 'filled');
xlabel('Set of initial values')
ylabel('Calls')
set(gca, 'FontSize', 12)
hold off
legend({'BFGS', 'Line Search', 'Trust Region', 'CasADI', 'fmincon'}, 'location', 'northeast')

subplot(1,3,2)
hold on
scatter(2, prob1.fcalls(2, 1), 30,  'red', 'filled');
scatter(2, prob1.fcalls(2, 2), 30, 'blue', 'filled');
scatter(2, prob1.fcalls(2, 3), 30, 'cyan', 'filled');
scatter(2, prob1.fcalls(2, 4), 30, 'yellow', 'filled');
scatter(2, prob1.fcalls(2, 5), 30, 'black', 'filled');
title('Number of function calls for the methods - Problem 1')
xlabel('Set of initial values')
ylabel('Calls')
set(gca, 'FontSize', 12)
hold off

subplot(1,3,3)
hold on
scatter(3, prob1.fcalls(3, 1), 30,  'red', 'filled');
scatter(3, prob1.fcalls(3, 2), 30, 'blue', 'filled');
scatter(3, prob1.fcalls(3, 3), 30, 'cyan', 'filled');
scatter(3, prob1.fcalls(3, 4), 30, 'yellow', 'filled');
scatter(3, prob1.fcalls(3, 5), 30, 'black', 'filled');
xlabel('Set of initial values')
ylabel('Calls')
set(gca, 'FontSize', 12)
hold off

figure
subplot(1,3,1)
hold on;
scatter(1, prob2.fcalls(1, 1), 30,  'red', 'filled');
scatter(1, prob2.fcalls(1, 2), 30, 'blue', 'filled');
scatter(1, prob2.fcalls(1, 3), 30, 'cyan', 'filled');
scatter(1, prob2.fcalls(1, 4), 30, 'yellow', 'filled');
scatter(1, prob2.fcalls(1, 5), 30, 'black', 'filled');
xlabel('Set of initial values')
ylabel('Calls')
set(gca, 'FontSize', 12)
hold off
legend({'BFGS', 'Line Search', 'Trust Region', 'CasADI', 'fmincon'}, 'location', 'northeast')

subplot(1,3,2)
hold on
scatter(2, prob2.fcalls(2, 1), 30,  'red', 'filled');
scatter(2, prob2.fcalls(2, 2), 30, 'blue', 'filled');
scatter(2, prob2.fcalls(2, 3), 30, 'cyan', 'filled');
scatter(2, prob2.fcalls(2, 4), 30, 'yellow', 'filled');
scatter(2, prob2.fcalls(2, 5), 30, 'black', 'filled');
title('Number of function calls for the methods - Problem 2')
xlabel('Set of initial values')
ylabel('Calls')
set(gca, 'FontSize', 12)
hold off

subplot(1,3,3)
hold on
scatter(3, prob2.fcalls(3, 1), 30,  'red', 'filled');
scatter(3, prob2.fcalls(3, 2), 30, 'blue', 'filled');
scatter(3, prob2.fcalls(3, 3), 30, 'cyan', 'filled');
scatter(3, prob2.fcalls(3, 4), 30, 'yellow', 'filled');
scatter(3, prob2.fcalls(3, 5), 30, 'black', 'filled');
xlabel('Set of initial values')
ylabel('Calls')
set(gca, 'FontSize', 12)
hold off

figure
subplot(1,3,1)
hold on;
scatter(1, prob3.fcalls(1, 1), 30,  'red', 'filled');
scatter(1, prob3.fcalls(1, 2), 30, 'blue', 'filled');
scatter(1, prob3.fcalls(1, 3), 30, 'cyan', 'filled');
scatter(1, prob3.fcalls(1, 4), 30, 'yellow', 'filled');
scatter(1, prob3.fcalls(1, 5), 30, 'black', 'filled');
xlabel('Set of initial values')
ylabel('Calls')
set(gca, 'FontSize', 12)
hold off
legend({'BFGS', 'Line Search', 'Trust Region', 'CasADI', 'fmincon'}, 'location', 'northeast')

subplot(1,3,2)
hold on
scatter(2, prob3.fcalls(2, 1), 30,  'red', 'filled');
scatter(2, prob3.fcalls(2, 2), 30, 'blue', 'filled');
scatter(2, prob3.fcalls(2, 3), 30, 'cyan', 'filled');
scatter(2, prob3.fcalls(2, 4), 30, 'yellow', 'filled');
scatter(2, prob3.fcalls(2, 5), 30, 'black', 'filled');
title('Number of function calls for the methods - Problem 3')
xlabel('Set of initial values')
ylabel('Calls')
set(gca, 'FontSize', 12)
hold off

subplot(1,3,3)
hold on
scatter(3, prob3.fcalls(3, 1), 30,  'red', 'filled');
scatter(3, prob3.fcalls(3, 2), 30, 'blue', 'filled');
scatter(3, prob3.fcalls(3, 3), 30, 'cyan', 'filled');
scatter(3, prob3.fcalls(3, 4), 30, 'yellow', 'filled');
scatter(3, prob3.fcalls(3, 5), 30, 'black', 'filled');
xlabel('Set of initial values')
ylabel('Calls')
set(gca, 'FontSize', 12)
hold off

%% Plotting - CPU time

figure
subplot(1,3,1)
hold on;
scatter(1, prob1.cpu(1, 1), 30,  'red', 'filled');
scatter(1, prob1.cpu(1, 2), 30, 'blue', 'filled');
scatter(1, prob1.cpu(1, 3), 30, 'cyan', 'filled');
scatter(1, prob1.cpu(1, 4), 30, 'yellow', 'filled');
scatter(1, prob1.cpu(1, 5), 30, 'black', 'filled');
xlabel('Set of initial values')
ylabel('Time (s)')
set(gca, 'FontSize', 12)
hold off
legend({'BFGS', 'Line Search', 'Trust Region', 'CasADI', 'fmincon'}, 'location', 'northeast')

subplot(1,3,2)
hold on
scatter(2, prob1.cpu(2, 1), 30,  'red', 'filled');
scatter(2, prob1.cpu(2, 2), 30, 'blue', 'filled');
scatter(2, prob1.cpu(2, 3), 30, 'cyan', 'filled');
scatter(2, prob1.cpu(2, 4), 30, 'yellow', 'filled');
scatter(2, prob1.cpu(2, 5), 30, 'black', 'filled');
title('Computation time for the methods - Problem 1')
xlabel('Set of initial values')
ylabel('Time (s)')
set(gca, 'FontSize', 12)
hold off

subplot(1,3,3)
hold on
scatter(3, prob1.cpu(3, 1), 30,  'red', 'filled');
scatter(3, prob1.cpu(3, 2), 30, 'blue', 'filled');
scatter(3, prob1.cpu(3, 3), 30, 'cyan', 'filled');
scatter(3, prob1.cpu(3, 4), 30, 'yellow', 'filled');
scatter(3, prob1.cpu(3, 5), 30, 'black', 'filled');
xlabel('Set of initial values')
ylabel('Time (s)')
set(gca, 'FontSize', 12)
hold off

figure
subplot(1,3,1)
hold on;
scatter(1, prob2.cpu(1, 1), 30,  'red', 'filled');
scatter(1, prob2.cpu(1, 2), 30, 'blue', 'filled');
scatter(1, prob2.cpu(1, 3), 30, 'cyan', 'filled');
scatter(1, prob2.cpu(1, 4), 30, 'yellow', 'filled');
scatter(1, prob2.cpu(1, 5), 30, 'black', 'filled');
xlabel('Set of initial values')
ylabel('Time (s)')
set(gca, 'FontSize', 12)
hold off
legend({'BFGS', 'Line Search', 'Trust Region', 'CasADI', 'fmincon'}, 'location', 'northeast')

subplot(1,3,2)
hold on
scatter(2, prob2.cpu(2, 1), 30,  'red', 'filled');
scatter(2, prob2.cpu(2, 2), 30, 'blue', 'filled');
scatter(2, prob2.cpu(2, 3), 30, 'cyan', 'filled');
scatter(2, prob2.cpu(2, 4), 30, 'yellow', 'filled');
scatter(2, prob2.cpu(2, 5), 30, 'black', 'filled');
title('Computation time for the methods - Problem 2')
xlabel('Set of initial values')
ylabel('Time (s)')
set(gca, 'FontSize', 12)
hold off

subplot(1,3,3)
hold on
scatter(3, prob2.cpu(3, 1), 30,  'red', 'filled');
scatter(3, prob2.cpu(3, 2), 30, 'blue', 'filled');
scatter(3, prob2.cpu(3, 3), 30, 'cyan', 'filled');
scatter(3, prob2.cpu(3, 4), 30, 'yellow', 'filled');
scatter(3, prob2.cpu(3, 5), 30, 'black', 'filled');
xlabel('Set of initial values')
ylabel('Time (s)')
set(gca, 'FontSize', 12)
hold off

figure
subplot(1,3,1)
hold on;
scatter(1, prob3.cpu(1, 1), 30,  'red', 'filled');
scatter(1, prob3.cpu(1, 2), 30, 'blue', 'filled');
scatter(1, prob3.cpu(1, 3), 30, 'cyan', 'filled');
scatter(1, prob3.cpu(1, 4), 30, 'yellow', 'filled');
scatter(1, prob3.cpu(1, 5), 30, 'black', 'filled');
xlabel('Set of initial values')
ylabel('Time (s)')
set(gca, 'FontSize', 12)
hold off
legend({'BFGS', 'Line Search', 'Trust Region', 'CasADI', 'fmincon'}, 'location', 'northeast')

subplot(1,3,2)
hold on
scatter(2, prob3.cpu(2, 1), 30,  'red', 'filled');
scatter(2, prob3.cpu(2, 2), 30, 'blue', 'filled');
scatter(2, prob3.cpu(2, 3), 30, 'cyan', 'filled');
scatter(2, prob3.cpu(2, 4), 30, 'yellow', 'filled');
scatter(2, prob3.cpu(2, 5), 30, 'black', 'filled');
title('Computation time for the methods - Problem 3')
xlabel('Set of initial values')
ylabel('Time (s)')
set(gca, 'FontSize', 12)
hold off

subplot(1,3,3)
hold on
scatter(3, prob3.cpu(3, 1), 30,  'red', 'filled');
scatter(3, prob3.cpu(3, 2), 30, 'blue', 'filled');
scatter(3, prob3.cpu(3, 3), 30, 'cyan', 'filled');
scatter(3, prob3.cpu(3, 4), 30, 'yellow', 'filled');
scatter(3, prob3.cpu(3, 5), 30, 'black', 'filled');
xlabel('Set of initial values')
ylabel('Time (s)')
set(gca, 'FontSize', 12)
hold off

%% Solutions Problem 1
figure
subplot(1,3,1)
hold on;
scatter(prob1.xopt(1,1), prob1.xopt(1,2), 30,  'red', 'filled');
scatter(prob1.xopt(1,3), prob1.xopt(1,4), 30, 'blue', 'filled');
scatter(prob1.xopt(1,5), prob1.xopt(1,6), 30, 'cyan', 'filled');
scatter(prob1.xopt(1,7), prob1.xopt(1,8), 30, 'yellow', 'filled');
scatter(prob1.xopt(1,9), prob1.xopt(1,10), 30, 'black', 'filled');
xlabel('x_1')
ylabel('x_2')
xlim([2.99 3.01])
ylim([1.99 2.01])
set(gca, 'FontSize', 12)
hold off
legend({'BFGS', 'Line Search', 'Trust Region', 'CasADI', 'fmincon'}, 'location', 'northeast')

subplot(1,3,2)
hold on;
scatter(prob1.xopt(2,1), prob1.xopt(2,2), 30,  'red', 'filled');
scatter(prob1.xopt(2,3), prob1.xopt(2,4), 30, 'blue', 'filled');
scatter(prob1.xopt(2,5), prob1.xopt(2,6), 30, 'cyan', 'filled');
scatter(prob1.xopt(2,7), prob1.xopt(2,8), 30, 'yellow', 'filled');
scatter(prob1.xopt(2,9), prob1.xopt(2,10), 30, 'black', 'filled');
xlabel('x_1')
ylabel('x_2')
xlim([-3.7 -3.6])
ylim([2.7 2.8])
title('Solutions for problem 1')
set(gca, 'FontSize', 12)
hold off

subplot(1,3,3)
hold on;
scatter(prob1.xopt(3,1), prob1.xopt(3,2), 30,  'red', 'filled');
scatter(prob1.xopt(3,3), prob1.xopt(3,4), 30, 'blue', 'filled');
scatter(prob1.xopt(3,5), prob1.xopt(3,6), 30, 'cyan', 'filled');
scatter(prob1.xopt(3,7), prob1.xopt(3,8), 30, 'yellow', 'filled');
scatter(prob1.xopt(3,9), prob1.xopt(3,10), 30, 'black', 'filled');
xlabel('x_1')
ylabel('x_2')
xlim([2.99 3.01])
ylim([1.99 2.01])
set(gca, 'FontSize', 12)
hold off
