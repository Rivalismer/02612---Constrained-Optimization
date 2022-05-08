% clc, clear, close all

% Setting quadprog not to show notifications
options_qp = optimoptions('quadprog', 'Display', 'off');

%% Input data - H, g, A, b, l, u matrices

H = [5.0000 1.8600 1.2400 1.4800 -0.4600
    1.8600 3.0000 0.4400 1.1200 0.5200
    1.2400 0.4400 3.8000 1.5600 -0.5400
    1.4800 1.1200 1.5600 7.2000 -1.1200
    -0.4600 0.5200 -0.5400 -1.1200 7.8000];

g = [-16.1000 
    -8.5000
    -15.7000
    -10.0200
    -18.6800];

A = [16.1000 1.0000 
    8.5000 1.0000
    15.7000 1.0000
    10.0200 1.0000
    18.6800 1.0000];

b = [15
    1];

l = zeros(5, 1);
u = ones(5, 1);

%% Initial guesses for estimated variables
x0 = [0.5; 0.5; 0.5; 0.5; 0.5];
y0 = [1; 1];
z0 = [0.5; 0.5; 0.5; 0.5; 0.5; 0.5; 0.5; 0.5; 0.5; 0.5];
s0 = [0.5; 0.5; 0.5; 0.5; 0.5; 0.5; 0.5; 0.5; 0.5; 0.5];

%% Calculating optimal value with use of different methods

% Primal dual interior point implemented by us
disp('Primal dual interior point')
tic
[x_primaldual, iter, converged] = PrimalDualInteriorPoint(x0, y0, s0, z0, H, g, A, b, l, u);
toc

% MATLAB quadoprog function without initial guess
disp('quadprog')
tic
[x_quadprog, ~, ~, stats] = quadprog(H, g, [], [], A', b, l, u, [], options_qp);
toc

% MATLAB quadoprog function with same initial guess as interior point
disp('quadprog')
tic
[x_quadprog, ~, ~, stats] = quadprog(H, g, [], [], A', b, l, u, x0, options_qp);
toc