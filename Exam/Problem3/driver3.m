clc, clear, close all
%{
    This program is performint test for implemented primal-dual interior
    point algorithm for linear programming. Test is performed on described 
    linear programming problem.
%}
% Turn off linprog notifications
options = optimoptions('linprog', 'Display', 'off');

%% Input data: g, A, b, l, u matrices

g = [-16.1000
     -8.5000
     -15.7000
     -10.0200
     -18.6800];
 
A = ones(5, 1);

b = 1;

l = zeros(5, 1);

u = ones(5, 1);

%% Initial guesses for estimated variables
x0 = [0.5; 0.5; 0.5; 0.5; 0.5];
mu0 = 0.5;
lambda0 = [0.5; 0.5; 0.5; 0.5; 0.5];

%% Calculating actual minimizer with use of linprog function
x_linprog = linprog(g, [], [], A', b, l, u, options);

[x, iter, converged] = PrimalDualInteriorPoint(x0, mu0, lambda0, g, A, b, l, u);