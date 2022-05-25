%clc, clear, close all
%{
    This program is testing implemented primal-dual interior point algorithm
    for equaly constrained linear program with lower and upper bounds. Linear
    problem is defined, solution is calculated with linprog and than
    with pdip algorithm.
%}

%% Defining testing problem

% Problem 1
prob(1).g = [-1
     -1/3];
 
prob(1).A = [1
     0.25];

prob(1).b = 0.5;

prob(1).l = -ones(2, 1);

prob(1).u = 2*ones(2, 1);

% Problem 2 
% random problem with 10 variables
prob(2).g = randi([-10, 15], [10, 1]);

prob(2).A = round(randn(10, 5), 2);

prob(2).b = round(randn(5, 1), 1);

prob(2).l = zeros(10, 1);
prob(2).u = ones(10, 1);

%% Initial guesses for PDIP algorithm
prob(1).x0 = [1.2; 1.2];
prob(1).mu0 = 0.5;
prob(1).lambda0 = [0.5; 0.5];

prob(2).x0 = [0.5; 0.5; 0.5; 0.5; 0.5; 0.5; 0.5; 0.5; 0.5; 0.5];
prob(2).mu0 = [0.5; 0.5; 0.5; 0.5; 0.5];
prob(2).lambda0 = [0.5; 0.5; 0.5; 0.5; 0.5; 0.5; 0.5; 0.5; 0.5; 0.5];

%% Calculating solutions
for i = [1, 2]
    % With linprog
    x_linprog = linprog(prob(i).g, [], [], prob(i).A', prob(i).b, prob(i).l, prob(i).u);

    % Calculiting minimizer with use of implemented fuction
    tic
    [x_opt, iter, converged] = PrimalDualInteriorPoint(prob(i).x0, prob(i).mu0, prob(i).lambda0, prob(i).g, prob(i).A, prob(i).b, prob(i).l, prob(i).u);
    cputime = toc;
    if converged
        disp('Algorithm converged to minimizer:')
    else
        disp('No convergence was obtained.')
    end
fprintf('Problem: %i \n', i);
disp(x_opt')
disp('Performance statistics:')
fprintf('Number of iterations: %i \n', iter)
fprintf('Computation time: %f seconds', cputime)
end
