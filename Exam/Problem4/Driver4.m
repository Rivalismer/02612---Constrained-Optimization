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