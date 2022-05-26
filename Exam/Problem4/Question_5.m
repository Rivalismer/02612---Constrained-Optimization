clear all, close all, clc
%{
    This program is solving three defined nonlinear programming problems
    with use of commercial symbolic solvers - CadADI and fmincon
%}
% Using CasADI solver
import casadi.*

%% Himmelblau problem - CasADI

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

