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
y0 = 0.5;
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

%% Question 5
clc, clear, close all
%{
    This program is comparing implemented primal-dual interior point
    algorithm for QP with Matlab function quadprog. 
%}

% Setting quadprog not to show notifications
options_qp = optimoptions('quadprog', 'Display', 'off');
methods = {'PDIP', 'quadprog', 'quadprog_x0'};

%% Defining quadratic problems: H, g, A, b, l, u

% Problem 1
qp(1).H = [1 -1 1
           -1 2 -2
           1 -2 4];
qp(1).g = [2; -3; 1];

qp(1).A = ones(3,1); 
qp(1).b = 1/2;

qp(1).l = zeros(3, 1);
qp(1).u = ones(3, 1);

% Problem 2
qp(2).H = [1 -1
           -1 2];
qp(2).g = [-2; -6];

qp(2).A = [1; 1];
qp(2).b = 0;

qp(2).l = [-1; -1];
qp(2).u = [1; 1];

% Problem 3
H_temp = rand(10, 10);
qp(3).H = tril(H_temp)+tril(H_temp)'; % Symmetric
qp(3).g = rand(10, 1);

qp(3).A = rand(10, 4);
qp(3).b = rand(4, 1);

qp(3).l = -1*ones(10, 1);
qp(3).u = ones(10, 1);

%% Initial guesses for estimated variables

% For problem 1:
qp(1).x0 = [0.2 0.1 0.6
            0.2 0.2 0.4
            0.2 0.2 -0.1];
qp(1).y0 = 0.5;
qp(1).s0 = [0.5; 0.5; 0.5; 0.5; 0.5; 0.5];
qp(1).z0 = [0.5; 0.5; 0.5; 0.5; 0.5; 0.5];

% For problem 2:
qp(2).x0 = [0.2 0.5 0.1
            0.2 0.8 -0.4];
qp(2).y0 = 0.5;
qp(2).s0 = [0.5; 0.5; 0.5; 0.5];
qp(2).z0 = [0.5; 0.5; 0.5; 0.5];

% For problem 3;
qp(3).x0 = [0.5*ones(10, 1), 0.2*ones(10, 1), randi([-1, 1], 10)*ones(10, 1)/10];
qp(3).y0 = [1; 1; 1; 1];
qp(3).z0 = 0.5*ones(20, 1);
qp(3).s0 = qp(3).z0;

%% Calling algorithms and saving performance for optimization problems

% Loop for testing algorithms and saving results on both problems
for i = [1, 2, 3]
   for j = [1, 2, 3]
      % continue the loop even if algorithm brakes
      % Primal dual interior point 
      try 
          tic
          [x_pdip, iter_pdip, converged] = PrimalDualInteriorPoint(qp(i).x0(:, j), qp(i).y0, ... 
                                                                    qp(i).s0, qp(i).z0, qp(i).H, qp(i).g, ...
                                                                    qp(i).A, qp(i).b, qp(i).l, qp(i).u);
          cputime_pdip = toc;
      catch
          fprintf('Primal-dual interior point\nError occured in problem %i, for starting point %i \n', i, j)
          x_pdip = zeros(length(qp(i).x0(:, j)),1);
          iter_pdip = 0; 
          cputime_pdip = 0;
      end
      % quadprog with initial guess the same as for primal-dual interior point
      try 
          tic
          [x_qp, ~, ~, stats_qp] = quadprog(qp(i).H, qp(i).g, [], [], qp(i).A', qp(i).b, ...
                                            qp(i).l, qp(i).u, qp(i).x0(:, j), options_qp);
          cputime_qp = toc;
      catch 
          fprintf('Quadprog with initial point\nError occured in problem %i, for starting point %i \n', i, j)
          x_qp = zeros(length(qp(i).x0(:, j)),1);
          stats_qp.iterations = 0; 
          cputime_qp = 0;
      end
      % No initial quess quadprog
      tic 
        [x_qpn, ~, ~, stats_qpn] = quadprog(qp(i).H, qp(i).g, [], [], qp(i).A', qp(i).b, ...
                                            qp(i).l, qp(i).u, [], options_qp);
      cputime_qpn = toc;
      % Save results in struct
      stats(i, j).iterations = [methods; num2cell([iter_pdip, stats_qpn.iterations, stats_qp.iterations])];
      stats(i, j).cputime = [methods; num2cell([cputime_pdip, cputime_qpn, cputime_qp])];
      stats(i, j).xopt = [methods; num2cell([x_pdip, x_qpn, x_qp])];
   end
end

%% Question 6
clc, clear, close all
%{
    This program is investigating set quadratic problem, while changing one
    value in  b  matrix of equality constraints. It is using primal-dual
    interiot point, quadprog and cvx gurobi optimization algorithms.

    For this program to work, cvx libary needs to be configured for matlab
%}
% Setting quadprog not to show notifications
options_qp = optimoptions('quadprog', 'Display', 'off');
% Setting cvx solver for Gurobi
cvx_solver Gurobi
methods = {'PDIP', 'quadprog', 'Gurobi'};

%% Input data: H, g, A, b, l, u matrices

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

% Vector of investigated b values
bs = [8.5 : 0.2 : 18.68, 18.68];

%% Initial guesses for estimated variables
x0 = [0.5; 0.5; 0.5; 0.5; 0.5];
y0 = [1; 1];
z0 = [0.5; 0.5; 0.5; 0.5; 0.5; 0.5; 0.5; 0.5; 0.5; 0.5];
s0 = [0.5; 0.5; 0.5; 0.5; 0.5; 0.5; 0.5; 0.5; 0.5; 0.5];

%% Calculating optimal value with use of different methods

% Loop for algorithms performance investigation on considered problem, the
% loop will continue to run even though algorithm would show error
for i = 1:1:length(bs)
    b(1) = bs(i);
    % Primal dual interior point algorithm
    try    
        tic
        [x_pdip, iter_pdip] = PrimalDualInteriorPoint(x0, y0, s0, z0, H, g, A, b, l, u);
        cputime_pdip = toc;
    catch
        fprintf('Primal-dual interior point\nError occured for b = %f', bs(i))
        x_pdip = zeros(5,1);
        iter_pdip = 0; 
        cputime_pdip = 0;
    end
    % Quadprog
    try
        tic
        [x_qp, ~, ~, stats_qp] = quadprog(H, g, [], [], A', b, l, u, [], options_qp);
        cputime_qp = toc;
    catch 
        fprintf('Quadprog\nError occured for b = %f', bs(i))
        x_qp = zeros(5,1);
        stats_qp.iterations = 0; 
        cputime_qp = 0;
    end
    % CVX Gurobi
    try
        tic
        cvx_begin
            variable x_gb(5)
            minimize(1/2*x_gb'*H*x_gb + g'*x_gb)
            subject to
                A' * x_gb == b
                x_gb <= u
                x_gb >= l
        cvx_end 
        cputime_gb = toc;
    catch
        printf('Gurobi\nError occured for b = %f', bs(i))
        x_gb = zeros(5,1);
        iter_gb = 0; 
        cputime_gb = 0;
    end
    % Saving results for all algorithms
    stats.x_pdip(:, i) = x_pdip;    % primal-dual minimizer
    stats.x_qp(:, i) = x_qp;        % quadprog minimizer
    stats.x_gb(:, i) = x_gb;        % gurobi minimizer
    stats.iter(i, :) = [iter_pdip, stats_qp.iterations];  % iterations number 
    stats.cpu(i, :) = [cputime_pdip, cputime_qp, cputime_gb]; % cpu time
    stats.opt_val(i) = cvx_optval;  % optimum value
end

%% Plotting results

figure(1);
scatter(bs, stats.opt_val, 30, 'filled');
title('Optimal objective function value');
xlabel('b(1) value');
ylabel('Optimal value');
set(gca, "FontSize", 12);

figure(2);
scatter(bs, vecnorm(stats.x_gb), 30, 'filled');
title('L-2 norm of reached minimizer');
xlabel('b(1) value');
ylabel('x_o_p_t norm');
set(gca, "FontSize", 12);

figure(3);
hold on;
scatter(bs, stats.cpu(:, 1), 30, 'r', 'filled');
scatter(bs, stats.cpu(:, 2), 30, 'm', 'filled');
scatter(bs, stats.cpu(:, 3), 30, 'g', 'filled');
title('Computation time for all methods');
xlabel('b(1) value');
ylabel('Time (s)');
set(gca, "FontSize", 12);
legend({'primal-dual','quadprog', 'gurobi'}, 'Location', 'best');
hold off;

figure(4);
hold on;
scatter(bs, stats.iter(:, 1), 30, 'filled');
scatter(bs, stats.iter(:, 2), 30, 'filled');
title('Number of iterations for selected methods');
xlabel('b(1) value');
ylabel('Iterations number');
set(gca, "FontSize", 12);
legend({'primal-dual','quadprog'}, 'Location', 'best')
hold off;

% Plotting separately all x dimensions
figure(5);
hold on;
% x1
subplot(2, 6, [1,2]);
plot(bs, stats.x_qp(1, :));
title('x_1');
xlabel('b(1)');
ylabel('x_1');
% x2
subplot(2, 6, [3,4]);
plot(bs, stats.x_qp(2, :));
title('x_2');
xlabel('b(1)');
ylabel('x_2');
% x3
subplot(2, 6, [5,6]);
plot(bs, stats.x_qp(3, :));
title('x_3');
xlabel('b(1)');
ylabel('x_3');
% x4
subplot(2, 6, [8,9]);
plot(bs, stats.x_qp(4, :));
title('x_4');
xlabel('b(1)');
ylabel('x_4');
% x5
subplot(2, 6, [10,11]);
plot(bs, stats.x_qp(5, :));
title('x_5');
xlabel('b(1)');
ylabel('x_5');
