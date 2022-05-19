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
methods = {'PDIP', 'quadprog', 'Gurobi'}

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