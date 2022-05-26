%% Markowitz portfolio as Quadratic Program
R = 12; % Our desired return

H = [2.50 0.93 0.62 0.74 -0.23;
    0.93 1.50 0.22 0.56 0.26;
    0.62 0.22 1.90 0.78 -0.27;
    0.74 0.56 0.78 3.60 -0.56;
    -0.23 0.26 -0.27 -0.56 3.90];


A = [16.10; 8.50; 15.70; 10.02; 18.68]';
beq = [1; R];
Aeq = [ones(5,1)'; A];
lb = zeros(5,1);
ub = ones(5,1);
f = zeros(5,1);

options = optimoptions('quadprog','Algorithm','interior-point-convex',...
    'Display', 'off', 'TolFun', 1e-10);

[x,fval,exitflag,output,lambda] = quadprog(H,f,[],[],Aeq,beq,lb,ub,[],options);
disp(x)
disp(2*fval) % Risk

%% Efficient Frontier
% Min return is 0, but it can't compute that as it makes the problem
% infeasible
% Max return is max(A)
Rmin = min(A);
is = zeros((max(A) - Rmin)/0.001, 1);
fs = zeros((max(A) - Rmin)/0.001, 1);
xs = zeros((max(A) - Rmin)/0.001, 5);
cnt = 1;
for i = Rmin:0.001:max(A)
    beq = [1; i];
    [x,fval,exitflag,output,lambda] = quadprog(H,f,[],[],Aeq,beq,lb,ub,[],options);
    is(cnt) = i;
    % quadprog considers problem in form 1/2x'Hx, thus optimal value
    % extracted is 2 times smaller than actual value of risk, which is x'Hx
    fs(cnt) = 2*fval;
    xs(cnt, :) = x;
    cnt = cnt + 1;
end
% Locating optimal portfolio
opt_return = is(find(fs==min(fs)));

% Preparing plots
figure(1)
hold on;
plot(is, fs, "LineWidth", 2)
scatter(opt_return, min(fs), 30, 'filled')
xlim([Rmin max(A)])
title('Efficient Frontier')
xlabel('Return')
ylabel('Risk')
set(gca, "FontSize", 12);
hold off;
legend(["Efficient frontier", "Optimal portfolio"], 'location', 'best')

figure
subplot(2,3,1)
plot(is,xs(:,1));
xlim([Rmin max(A)])
ylim([0 1])
title('Asset 1')
xlabel('Return')
ylabel('Portfolio % (in decimals)')

subplot(2,3,2)
plot(is,xs(:,2));
xlim([Rmin max(A)])
ylim([0 1])
title('Asset 2')
xlabel('Return')
ylabel('Portfolio % (in decimals)')

subplot(2,3,3)
plot(is,xs(:,3));
xlim([Rmin max(A)])
ylim([0 1])
title('Asset 3')
xlabel('Return')
ylabel('Portfolio % (in decimals)')

subplot(2,3,4)
plot(is,xs(:,4));
xlim([Rmin max(A)])
ylim([0 1])
title('Asset 4')
xlabel('Return')
ylabel('Portfolio % (in decimals)')

subplot(2,3,5)
plot(is,xs(:,5));
xlim([Rmin max(A)])
ylim([0 1])
title('Asset 5')
xlabel('Return')
ylabel('Portfolio % (in decimals)')

%% Bi-criterion problem
clear all, close all, clc
%{
    This driver is used for bi-criterion optimization of portfolio. Whole
    variety of solvers will be investigated to solve problem, with alpha
    parameter change
%}
% Options for solvers
options = optimoptions('quadprog','Algorithm','interior-point-convex',...
    'Display', 'off', 'TolFun', 1e-10);

%% Define problem

% Covariance matrix
H = [2.50 0.93 0.62 0.74 -0.23;
    0.93 1.50 0.22 0.56 0.26;
    0.62 0.22 1.90 0.78 -0.27;
    0.74 0.56 0.78 3.60 -0.56;
    -0.23 0.26 -0.27 -0.56 3.90];

% Returns vector
g = [16.10; 8.50; 15.70; 10.02; 18.68]';

% Constraints
lb = zeros(5,1);
ub = ones(5,1);
alpha = linspace(0,1,101);
beq = 1;
Aeq = ones(5,1)';

%################################################
%       DIFFERENT SOLVERS INVESTIGATION
%################################################

%% Interior-Point Linear
% Calculation only for alpha value equal to 1 since problem becomes linear
x0 = 0.2*ones(5,1);
lambda0 = 0.5*ones(5,1);
mu0 = 0.5;
tic
g = -g';
[x,iter,converged] = PrimalDualInteriorPointLP(x0,mu0,lambda0,g,Aeq',beq,lb,ub);
cputime = toc;
% Printing results
disp('Primal-dual interior point for linear programming:')
disp(x);
fprintf('Number of iterations: %i \n', iter)
fprintf('Computation time: %f s \n', cputime)

%% LU
for i=1:1:length(alpha)-1
    H_bi = (1-alpha(i))*H;
    g_bi = -alpha(i)*g;
    tic
    [x,lambda] = EqualityQPSolver(H_bi,g_bi, Aeq',beq,'LUDense');
    cputime = toc;
    stats.x_lu(i, :) = x;
    stats.return(i, 1) = -g'*x;
    stats.risk(i, 1) = x'*H*x;
    stats.cputime(i, 1) = cputime;
end

%% LDL
for i=1:1:length(alpha)-1
    H_bi = (1-alpha(i))*H;
    g_bi = -alpha(i)*g;
    tic
    [x,lambda] = EqualityQPSolver(H_bi,g_bi, Aeq',beq,'LDLDense');
    cputime = toc;
    stats.x_ldl(i, :) = x;
    stats.return(i, 2) = -g'*x;
    stats.risk(i, 2) = x'*H*x;
    stats.cputime(i, 2) = cputime;
end

%% Null-Space
for i=1:1:length(alpha)-1
    H_bi = (1-alpha(i))*H;
    g_bi = -alpha(i)*g;
    tic
    [x,lambda] = EqualityQPSolver(H_bi,g_bi, Aeq',beq,'NullSpace');
    cputime = toc;
    stats.x_rs(i, :) = x;
    stats.return(i, 3) = -g'*x;
    stats.risk(i, 3) = x'*H*x;
    stats.cputime(i, 3) = cputime;
end

%% Range-Space
for i=1:1:length(alpha)-1
    H_bi = (1-alpha(i))*H;
    g_bi = -alpha(i)*g;
    tic
    [x,lambda] = EqualityQPSolver(H_bi,g_bi, Aeq',beq,'RangeSpace');
    cputime = toc;
    stats.x_ns(i, :) = x;
    stats.return(i, 4) = g'*x;
    stats.risk(i, 4) = x'*H*x;
    stats.cputime(i, 4) = cputime;
end

%% Interior-Point quadratic
x0 = 0.2*ones(5,1);
y0 = 1;
s0 = ones(10,1);
z0 = ones(10,1);
for i=1:1:length(alpha)-1
    H_bi = (1-alpha(i))*H;
    g_bi = -alpha(i)*g;
    tic
    [x,iter,converged] = PrimalDualInteriorPoint(x0,y0,s0,z0,H_bi,g_bi,Aeq',beq,lb,ub);
    cputime = toc;
    stats.x_pdip(i, :) = x;
    stats.return(i, 5) = -g'*x;
    stats.risk(i, 5) = x'*H*x;
    stats.cputime(i, 5) = cputime;
    stats.iter(i, 1) = iter;
end

%% Quadprog

for i=1:1:length(alpha)-1
    H_bi = (1-alpha(i))*H;
    g_bi = -alpha(i)*g';
    tic
    [x,fval,~,output] = quadprog(H_bi,g_bi,[],[],Aeq,beq,lb,ub,[],options);
    cputime = toc;
    stats.x_qp(i, :) = x;
    stats.return(i, 6) = -g'*x;
    stats.risk(i, 6) = x'*H*x;
    stats.cputime(i, 6) = cputime;
    stats.iter(i, 2) = output.iterations;
    stats.fval(i) = fval;
end

%% CVX
for i = 1:length(alpha)-1
    tic
    cvx_begin
        variable x_c(5)
        minimize((1-alpha(i))*x_c'*H*x_c -alpha(i)*g'*x_c)
        subject to
            Aeq * x_c == beq
            x_c <= ub
            x_c >= lb
    cvx_end
    cputime = toc;
    stats.x_c(i, :) = x_c;
    stats.return(i, 7) = -g'*x_c;
    stats.risk(i, 7) = x_c'*H*x_c;
    stats.cputime(i, 7) = cputime;
end

%% Gurobi
cvx_solver Gurobi
for i = 1:length(alpha)-1
    tic
    cvx_begin
        variable x_c(5)
        minimize((1-alpha(i))*x_c'*H*x_c -alpha(i)*g'*x_c)
        subject to
            Aeq * x_c == beq
            x_c <= ub
            x_c >= lb
    cvx_end
    cputime = toc;
    stats.x_c(i, :) = x_c;
    stats.return(i, 8) = -g'*x_c;
    stats.risk(i, 8) = x_c'*H*x_c;
    stats.cputime(i, 8) = cputime;
end

%% Plotting

% Weight distribution for assets with varying alpha
figure
subplot(2,6,[1,2])
pl1 = plot(alpha(1:100),stats.x_qp(:,1));
xlim([min(alpha) max(alpha)])
ylim([0 1])
title('Asset 1')
xlabel('\alpha')
ylabel('Portfolio % (in decimals)')

subplot(2,6,[3,4])
pl2 = plot(alpha(1:100),stats.x_qp(:,2));
xlim([min(alpha) max(alpha)])
ylim([0 1])
title('Asset 2')
xlabel('\alpha')

subplot(2,6,[5,6])
pl3 = plot(alpha(1:100),stats.x_qp(:,3));
xlim([min(alpha) max(alpha)])
ylim([0 1])
title('Asset 3')
xlabel('\alpha')

subplot(2,6,[8,9])
pl4 = plot(alpha(1:100),stats.x_qp(:,4));
xlim([min(alpha) max(alpha)])
ylim([0 1])
title('Asset 4')
xlabel('\alpha')
ylabel('Portfolio % (in decimals)')

subplot(2,6,[10, 11])
pl5 = plot(alpha(1:100),stats.x_qp(:,5));
xlim([min(alpha) max(alpha)])
ylim([0 1])
title('Asset 5')
xlabel('\alpha')
set([pl1, pl2, pl3, pl4, pl5], 'LineWidth', 1.5)

% Optimal function value 

figure
plot(alpha(1:100), stats.fval, 'LineWidth', 2);
xlim([min(alpha) max(alpha)]);
title('Optimal objective function value')
xlabel('\alpha')
ylabel('F value')
set(gca, 'FontSize', 12)

% CPU times

figure
hold on;
scatter(alpha(1:100), stats.cputime(:, 5), 30, 'filled');
scatter(alpha(1:100), stats.cputime(:, 6), 30, 'filled');
scatter(alpha(1:100), stats.cputime(:, 7), 30, 'filled');
scatter(alpha(1:100), stats.cputime(:, 8), 30, 'filled');
xlim([min(alpha) max(alpha)]);
title('Computation time for quadratic programming')
xlabel('\alpha')
ylabel('Time (s)')
legend({'Primal dual interior point', 'quadprog', 'SDPT3', 'Gurobi'}, 'location', 'best')
set(gca, 'FontSize', 12)


% Number of iterations

figure
hold on;
scatter(alpha(1:100), stats.iter(:, 1), 30, 'filled');
scatter(alpha(1:100), stats.iter(:, 2), 30, 'filled');
xlim([min(alpha) max(alpha)]);
title('Number of iterations for qp methods')
xlabel('\alpha')
ylabel('Iterations')
set(gca, 'FontSize', 12)
hold off
legend({'Primal dual interior point', 'quadprog'}, 'location', 'northeast')

% Return-risk trade-off

figure 
plot(stats.return(:, 5), stats.risk(:,5), 'LineWidth', 2)
title('Risk-return trade-off')
xlabel('Return')
ylabel('Risk')
set(gca, 'FontSize', 12)

% Primal dual computation time

figure
scatter(alpha(1:100), stats.cputime(:, 5), 30, 'filled');
xlim([min(alpha) max(alpha)]);
title('Computation time for primal dual algorithm')
xlabel('\alpha')
ylabel('Time')
set(gca, 'FontSize', 12)

% Primal dual iterations number

figure
scatter(alpha(1:100), stats.iter(:, 1), 30, 'filled');
xlim([min(alpha) max(alpha)]);
title('Number of iterations for primal dual algorithm')
xlabel('\alpha')
ylabel('Iterations')
set(gca, 'FontSize', 12)

%% Risk free asset - 0 covariance with everything
R=12;
% Also rf = 0, so we just add a 0 to the end of our A
H = [2.50 0.93 0.62 0.74 -0.23 0;
    0.93 1.50 0.22 0.56 0.26 0;
    0.62 0.22 1.90 0.78 -0.27 0;
    0.74 0.56 0.78 3.60 -0.56 0;
    -0.23 0.26 -0.27 -0.56 3.90 0
    0 0 0 0 0 0];
A = [16.10; 8.50; 15.70; 10.02; 18.68; 0]';
beq = [1; R];
Aeq = [ones(6,1)'; A];
lb = zeros(6,1);
ub = ones(6,1);
f = zeros(6,1);

[x,fval,exitflag,output,lambda] = quadprog(H,f,[],[],Aeq,beq,lb,ub,[],options);
disp(x)
disp(fval) % Risk

%% Risk free asset - Efficient frontier
Rmin = min(A);
is = zeros((max(A) - Rmin)/0.001, 1);
fs = zeros((max(A) - Rmin)/0.001, 1);
xs = zeros((max(A) - Rmin)/0.001, 6);
cnt = 1;
for i = Rmin:0.001:max(A)
    beq = [1; i];
    [x,fval,exitflag,output,lambda] = quadprog(H,f,[],[],Aeq,beq,lb,ub,[],options);
    is(cnt) = i;
    fs(cnt) = fval;
    xs(cnt, :) = x;
    cnt = cnt + 1;
end
figure
plot(is, fs, "LineWidth", 2)
xlim([Rmin max(A)])
title('Efficient Frontier')
xlabel('Return')
ylabel('Risk')
set(gca, "FontSize", 12);

figure
subplot(2,3,1)
pl1 = plot(is,xs(:,1));
xlim([Rmin max(A)])
ylim([0 1])
title('Asset 1')
xlabel('Return')
ylabel('Portfolio % (in decimals)')

subplot(2,3,2)
pl2 = plot(is,xs(:,2));
xlim([Rmin max(A)])
ylim([0 1])
title('Asset 2')
xlabel('Return')
ylabel('Portfolio % (in decimals)')

subplot(2,3,3)
pl3 = plot(is,xs(:,3));
xlim([Rmin max(A)])
ylim([0 1])
title('Asset 3')
xlabel('Return')
ylabel('Portfolio % (in decimals)')

subplot(2,3,4)
pl4 = plot(is,xs(:,4));
xlim([Rmin max(A)])
ylim([0 1])
title('Asset 4')
xlabel('Return')
ylabel('Portfolio % (in decimals)')

subplot(2,3,5)
pl6 = plot(is,xs(:,5));
xlim([Rmin max(A)])
ylim([0 1])
title('Asset 5')
xlabel('Return')
ylabel('Portfolio % (in decimals)')

subplot(2,3,6)
pl5 = plot(is, xs(:,6));
xlim([Rmin max(A)])
ylim([0 1])
title('Asset 6')
xlabel('Return')
ylabel('Portfolio % (in decimals)')
set([pl1, pl2, pl3, pl4, pl5, pl6], 'LineWidth', 1.5)

%% With R = 14
Rmin = min(A);
is = zeros((max(A) - Rmin)/0.001, 1);
fs = zeros((max(A) - Rmin)/0.001, 1);
xs = zeros((max(A) - Rmin)/0.001, 6);
cnt = 1;
for i = Rmin:0.001:max(A)
    beq = [1; i];
    [x,fval,exitflag,output,lambda] = quadprog(H,f,[],[],Aeq,beq,lb,ub,[],options);
    is(cnt) = i;
    fs(cnt) = fval;
    xs(cnt, :) = x;
    cnt = cnt + 1;
end

%Compute optimal with R = 14
R = 14;
beq = [1; R];
[x,fval,exitflag,output,lambda] = quadprog(H,f,[],[],Aeq,beq,lb,ub,[],options);

% Plots
figure
hold on
plot(is, fs, 'Linewidth', 2)
plot(14, 2*fval, 'r.','LineWidth',2,'MarkerSize',25)
hold off
xlim([Rmin max(A)])
title('Efficient Frontier')
xlabel('Return')
ylabel('Risk')
legend(["Efficient frontier", 'Optimal portfolio for R=14'], "location", "best")

figure
subplot(2,3,1)
hold on
pl1 = plot(is,xs(:,1));
plot(14,x(1),'r.','LineWidth',2,'MarkerSize',25);
hold off
xlim([Rmin max(A)])
ylim([0 1])
title('Asset 1')
xlabel('Return')
ylabel('Portfolio % (in decimals)')

subplot(2,3,2)
hold on
pl2 = plot(is,xs(:,2));
plot(14,x(2),'r.','LineWidth',2,'MarkerSize',25);
hold off
xlim([Rmin max(A)])
ylim([0 1])
title('Asset 2')
xlabel('Return')
ylabel('Portfolio % (in decimals)')

subplot(2,3,3)
hold on
pl3 = plot(is,xs(:,3));
plot(14,x(3),'r.','LineWidth',2,'MarkerSize',25);
hold off
xlim([Rmin max(A)])
ylim([0 1])
title('Asset 3')
xlabel('Return')
ylabel('Portfolio % (in decimals)')

subplot(2,3,4)
hold on
pl4 = plot(is,xs(:,4));
plot(14,x(4),'r.','LineWidth',2,'MarkerSize',25);
hold off
xlim([Rmin max(A)])
ylim([0 1])
title('Asset 4')
xlabel('Return')
ylabel('Portfolio % (in decimals)')

subplot(2,3,5)
hold on
pl5 = plot(is,xs(:,5));
plot(14,x(5),'r.','LineWidth',2,'MarkerSize',25);
hold off
xlim([Rmin max(A)])
ylim([0 1])
title('Asset 5')
xlabel('Return')
ylabel('Portfolio % (in decimals)')

subplot(2,3,6)
hold on
pl6 = plot(is, xs(:,6));
plot(14,x(6),'r.','LineWidth',2,'MarkerSize',25);
hold off
xlim([Rmin max(A)])
ylim([0 1])
title('Asset 6')
xlabel('Return')
ylabel('Portfolio % (in decimals)')
set([pl1, pl2, pl3, pl4, pl5, pl6], 'LineWidth', 1.5)