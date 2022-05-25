clc, clear, close all
%{
    This program is performing test for implemented primal-dual interior
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
%x0 = [0.5; 0.5; 0.5; 0.5; 0.5];
%x0 = [1.5; 1.5; 1.5; 1.5; 1.5];
%x0 = [0.5; 0.5; 0.5; 0.5; 1.1];
x0 = [0.2; 0.2; 0.2; 0.2; 1.0];
mu0 = 0.5;
lambda0 = [0.5; 0.5; 0.5; 0.5; 0.5];

%% Testing implemented algorithm

% Calculating actual minimizer with use of linprog function
x_linprog = linprog(g, [], [], A', b, l, u, options);

% Calculiting minimizer with use of implemented fuction
tic
[x_opt, iter, converged] = PrimalDualInteriorPoint(x0, mu0, lambda0, g, A, b, l, u);
cputime = toc;

if converged
    disp('Algorithm converged to minimizer:')
else
    disp('No convergence was obtained.')
end
disp(x_opt)
disp('Performance statistics:')
fprintf('Number of iterations: %i \n', iter)
fprintf('Computation time: %f seconds', cputime)