clc, clear, close all
%{
    Purpose of this program is to test all implemented equality constrained
    qp solvers on a defined problem. It is also calculating sensitivities
    for considered problem and using them to approximate values of
    optimizer in specified range.
%}

% Not showing quadprog displays
options_qp = optimoptions('quadprog', 'Display', 'off');
% List of solvers
solvers = {'LUDense', 'LUSparse', 'LDLDense', 'LDLSparse', 'RangeSpace', 'NullSpace'};

%% Input data - H, g, A, b matrices

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

%% Finding minimizer with quadprog

% Using quadprog to find actual minimizer
L = [H -A; -A' zeros(2,2)];
R = - [g ; b];
x_backslash = L\R;
[x_quadprog, ~, ~, ~, lambda_qp] = quadprog(H, g, [], [], A', b, [], [], [], options_qp);

%% Implemented solvers tests for single value

% Loop through each solver and retrieve solution
for i = 1:1:length(solvers)
    [x_opt, lambda] = EqualityQPSolver(H,g,A,b,solvers{i});
    x(:, i) = x_opt;
    lambda(:, i) = lambda;
end

% Solving obtained minimizers 
test_solvers = [solvers; num2cell(x)];

%% Calculating solutions for varying b(1) value (without sensitivities)

% Defining values
bs = [8.5:0.3:18.68, 18.68];

% Looping over all solvers and investigated values, even if some solvers
% won't reach mimimum, the loop will continue
for i = 1:1:length(solvers)
    for j = 1:1:length(bs)
        b(1) = bs(j);
        try
            [x_opt, lambda] = EqualityQPSolver(H,g,A,b,solvers{i});
        catch
            x_opt = NaN(5, 1);
            lambda = NaN(2, 1);
        end
        optim_solvers(:, j, i) = x_opt;
    end
end

% quadprog confimation of results
for i = 1:1:length(bs)
    b(1) = bs(i);
    true_minimizers(:, i) = quadprog(H, g, [], [], A', b, [], [], [], options_qp);
end

%% Calculating minimizer values with use of sensitivities

% Sizes variables
[m, ~] = size(H);
n = length(b);

% Sensitivities matrix including both x and lambda sensitivities
sensitivities = [zeros(n, m), -eye(n)]* ...
                inv([H, -A; -A', zeros(n)]);

% Extracting x and lambda sensitivities
x_sens = sensitivities(:, 1:m);
lambda_sens = sensitivities(:, m+1:end);

% Initial opzimizer for b(1) = 8.5
xp0 = true_minimizers(:, 1);

% Calculating optimal x values
for i = 1:1:length(bs)
    b(1) = bs(i);
    sens_optim(:, i) = xp0 + x_sens'*(b-[8.5; 1]);
end

