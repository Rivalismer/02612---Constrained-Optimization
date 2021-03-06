%% Problem 3 - Question 4
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
    fprintf('Linprog solution:\n');
    disp(x_linprog')
    % Calculating minimizer with use of implemented fuction
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
    fprintf('Computation time: %f seconds\n', cputime)
end


%% Question 5
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
test.x0(1,:) = [0.5; 0.5; 0.5; 0.5; 0.5];
test.x0(2,:) = [1.5; 1.5; 1.5; 1.5; 1.5];
test.x0(3,:) = [0.5; 0.5; 0.5; 0.5; 1.1];
test.x0(4,:) = [0.2; 0.2; 0.2; 0.2; 1.0];
mu0 = 0.5;
lambda0 = [0.5; 0.5; 0.5; 0.5; 0.5];

%% Testing implemented algorithm

% Calculating actual minimizer with use of linprog function
x_linprog = linprog(g, [], [], A', b, l, u, options);

% Calculiting minimizer with use of implemented fuction
for i=1:4
    try
        tic
        [x_opt, iter, converged] = PrimalDualInteriorPoint(test.x0(i,:)', mu0, lambda0, g, A, b, l, u);
        test.cputime(i) = toc;
        test.iter(i) = iter;
        test.xopt(i,:) = x_opt';
        if converged
            disp('Algorithm converged to minimizer:')
        else
            disp('No convergence was obtained.')
        end
    catch
        fprintf('Could not compute for i = %d\n', i);
    end
end


%% Question 6
clc, clear all, close all
%{
    This driver is comparing optimization methods: linprog, gurobi and
    primal dual interior point for linear programming. Two problems are
    created and solvers are going to be used to solve them.

    For this program to work, cvx libary needs to be configured for matlab
%}
% Setting linprog not to show notifications
options_lp = optimoptions('linprog', 'Display', 'off');
% Setting cvx solver for Gurobi
methods = {'PDIP', 'linprog', 'CVX'};

%% Defining problems

%Problem 1
prob(1).g = [-16.1000
     -8.5000
     -15.7000
     -10.0200
     -18.6800];
 
prob(1).A = ones(5, 1);

prob(1).b = 1;

prob(1).l = zeros(5, 1);
prob(1).u = ones(5, 1);

% Problem 2
prob(2).g = [-5; 13; 7; 2; 13; -8; 9; 9; 4; -6];

prob(2).A = [0.3300   -1.0900   -1.3400    1.2000    1.0800
   -0.5100   -1.4300    0.0300   -0.5900    0.9700
   -0.9000   -1.0100    0.8500   -0.4700   -0.5700
   -1.2000   -0.2100    0.4000    0.8900    0.8100
    1.0400   -0.3300   -0.7000   -1.3900    0.1700
   -0.8500    1.9400   -1.6300   -1.9600   -0.5100
   -0.1700   -0.5700    1.4600    0.4200   -1.1900
   -1.2100   -0.2500    2.0500    0.4000    0.6500
   -0.3000   -1.5700    0.1200    0.1000   -0.3500
   -3.2300   -0.4800   -0.9900    0.5000    0.0500];

prob(2).b = [-0.8; -1.6; 0.2; -0.1;  1.2];

prob(2).l = zeros(10, 1);
prob(2).u = ones(10, 1);

%% Initial guesses for estimated variables
prob(1).x0 = [0.5; 0.5; 0.5; 0.5; 0.5];
prob(2).x0 = [0.5; 0.5; 0.5; 0.5; 0.5; 0.5; 0.5; 0.5; 0.5; 0.5];

prob(1).mu0 = 0.5;
prob(2).mu0 = [0.5; 0.5; 0.5; 0.5; 0.5];

prob(1).lambda0 = [0.5; 0.5; 0.5; 0.5; 0.5];
prob(2).lambda0 = [0.5; 0.5; 0.5; 0.5; 0.5; 0.5; 0.5; 0.5; 0.5; 0.5];

% Loop for algorithms performance investigation on considered problem, the
% loop will continue to run even though algorithm would show error
for i = [1,2]
    % Primal dual interior point algorithm
    try    
        tic
        [x_pdip, iter_pdip] = PrimalDualInteriorPoint(prob(i).x0, prob(i).mu0, prob(i).lambda0, prob(i).g, prob(i).A, prob(i).b, prob(i).l, prob(i).u);
        cputime_pdip = toc;
    catch
        fprintf('Primal-dual interior point\nError occured for problem %i', i)
        x_pdip = NaN;
        iter_pdip = NaN; 
        cputime_pdip = NaN;
    end
    % Linprog
    try
        tic
        [x_lp, ~, ~, stats_lp] = linprog(prob(i).g, [], [], prob(i).A', prob(i).b, prob(i).l, prob(i).u, options_lp);
        cputime_lp = toc;
    catch 
        fprintf('Linprogprog\nError occured for problem %i', i)
        x_lp = NaN;
        stats_lp.iterations = NaN; 
        cputime_lp = NaN;
    end
    % CVX
    ni = length(prob(i).x0);
    try
        tic
        cvx_begin
            variable x_cvx(ni)
            minimize(prob(i).g'*x_cvx)
            subject to
                prob(i).A' * x_cvx == prob(i).b
                x_cvx <= prob(i).u
                x_cvx >= prob(i).l
        cvx_end 
        cputime_cvx = toc;
    catch
        printf('CVX\nError occured for problem %i', i)
        x_cvx = NaN;
        iter_cvx = NaN; 
        cputime_cvx = NaN;
    end
    % Saving results for all algorithms
    opts(i).x_pdip = x_pdip;    % primal-dual minimizer
    opts(i).x_lp = x_lp;        % quadprog minimizer
    opts(i).x_cvx = x_cvx;        % cvx minimizer
    stats.iter(i, :) = [iter_pdip, stats_lp.iterations];  % iterations number 
    stats.cpu(i, :) = [cputime_pdip, cputime_lp]; % cpu time - for cvx read from console
end
