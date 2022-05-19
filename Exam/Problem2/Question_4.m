clc, clear, close all
%{
    This program is testing implemented primal-dual interior point
    algorithm for qudratic problems. Example was taken from matlab quadprog
    function documentation. Url:
 
    https://www.mathworks.com/help/optim/ug/quadprog.html
    Section: Quadratic Minimization with Linear Constraints and Bounds

    Minimizer that should be obtained is:
    
    x = [0, 0.5, 0]'
%}
%% Defining quadratic problem: H, g, A, b, l, u

H = [1 -1 1
    -1 2 -2
    1 -2 4];

g = [2; -3; 1];

A = ones(3,1); 
b = 1/2;

l = zeros(3, 1);
u = ones(3, 1);

%% Initial guesses for estimated variables
x0 = [0.2; 0.2; 0.2];
y0 = [0.5];
s0 = [0.5; 0.5; 0.5; 0.5; 0.5; 0.5];
z0 = [0.5; 0.5; 0.5; 0.5; 0.5; 0.5];

%% Running primal-dual interior point algorithm
tic
[x_opt, iter, converged] = PrimalDualInteriorPoint(x0, y0, s0, z0, H, g, A, b, l, u);
cputime = toc;

%% Printing out the results
if converged
   disp('Algorithm converged to minimizer:')
else
    disp('No convergence was obtained.')
end
disp(x_opt)
disp('Performance statistics:')
fprintf('Number of iterations: %i \n', iter)
fprintf('Computation time: %f seconds', cputime)
