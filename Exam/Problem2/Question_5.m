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
%qp(2).A = [0; 0];
qp(2).b = 0;

qp(2).l = [-1; -1];
qp(2).u = [1; 1];

% Problem 3
H_temp = rand(10, 10);
qp(3).H = tril(H_temp)+tril(H_temp)';
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
qp(1).y0 = [0.5];
qp(1).s0 = [0.5; 0.5; 0.5; 0.5; 0.5; 0.5];
qp(1).z0 = [0.5; 0.5; 0.5; 0.5; 0.5; 0.5];

% For problem 2:
qp(2).x0 = [0.2 0.5 0.1
            0.2 0.8 -0.4];
qp(2).y0 = [0.5];
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
