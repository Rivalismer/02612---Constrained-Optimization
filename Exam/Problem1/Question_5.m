clc, clear, close all
%{
    This program is analysing size of quadratic equality constrained optimization
    problem impact on all implemented solvers. Problems will be generated 
    randomly. Firstly dense problems will be investigated, than sparse will
    be compared
%}
% Not showing quadprog displays
options_qp = optimoptions('quadprog', 'Display', 'off');
% Methods for used solvers
solvers = {'LUDense', 'LUSparse', 'LDLDense', 'LDLSparse', 'RangeSpace', 'NullSpace'};
% Range of sizes of investigated quadratic problems
sizes = [10:20:1000];

%% Creating problems for optimization 
% define random quadratic dense problems
for i = 1:1:length(sizes)
    H_temp = rand(sizes(i), sizes(i));
    prob(i).H = tril(H_temp)+tril(H_temp)';
    prob(i).g = round(randn(sizes(i), 1), 2);
    prob(i).A = round(randn(sizes(i), 10), 2);
    prob(i).b = round(randn(10, 1), 2);
end

% define random quadratic sparse problems
for i = 1:1:length(sizes)
    H_temp = rand(sizes(i), sizes(i));
    entries = randi([0 , 1], sizes(i), sizes(i));
    H_ent = H_temp.*entries;
    probs(i).H = tril(H_ent)+tril(H_ent)';
    
    entries = randi([0, 1], sizes(i), 1);
    values = round(randn(sizes(i), 1), 2);
    probs(i).g = entries.*values;
    
    entries = randi([0, 1], sizes(i), 10);
    values = round(randn(sizes(i), 10), 2);
    probs(i).A = entries.*values;
    
    entries = randi([0, 1], 10, 1);
    values = round(randn(10, 1), 2);
    probs(i).b = entries.*values;
end

%###############################################
%            DENSE INVESTIGATION
%###############################################

%% Finding solutions with use of quadprog - dense

for i = 1:1:length(sizes)
   try
       x_opt = quadprog(prob(i).H, prob(i).g, [], [], prob(i).A', prob(i).b, [], [], [], options_qp);
   catch
       x_opt = zeros(sizes(i), 1);
   end
   statsqp(i).x_opt = x_opt;
end

%% Trying all implemented algorithms and counting computation time - dense

for i = 1:1:length(solvers)
    for j = 1:1:length(sizes) 
       try 
           tic
            [x_opt, ~] = EqualityQPSolver(prob(j).H, prob(j).g, prob(j).A, prob(j).b, solvers{i});
           cputime = toc;
       catch
           x_opt = zeros(sizes(i), 1);
           cputime = 0;
       end
       % Saving results
       stats(i).xopt(j).x_opt = x_opt;
       stats(i).cputime(j) = cputime;
    end
end

%% Detecting where minimizer wasn't found - dense

for i = 1:1:length(sizes)
   for j = 1:1:length(solvers)
      if ~isempty(stats(j).xopt(i).x_opt) && ~isempty(statsqp(i).x_opt)
         if norm(stats(j).xopt(i).x_opt - statsqp(i).x_opt, inf) > 0.002 
            fprintf('Method %s shows weak performance for size %i \n', solvers{j}, sizes(i)); 
         end
      else
         fprintf('Method %s did not find any solution for size %i \n', solvers{j}, sizes(i)); 
      end
   end
end

%% Plotting results - dense

% Cpu time dense problems
figure(1)
hold on;
p1 = plot(sizes, stats(1).cputime);
p2 = plot(sizes, stats(2).cputime);
p3 = plot(sizes, stats(3).cputime);
p4 = plot(sizes, stats(4).cputime);
p5 = plot(sizes, stats(5).cputime);
p6 = plot(sizes, stats(6).cputime);
title("Computation time - dense analysis");
xlabel("Size of QP problem");
ylabel('Time (s)');
legend({'LUDense', 'LUSparse','LDLDense','LDLSparse','RangeSpace','NullSpace'}, 'location', 'best');
set(gca, "LineWidth", 2, 'FontSize', 12);
set([p1, p2, p3, p4, p5, p6], "LineWidth", 2);

%% 
%###############################################
%            SPARSE INVESTIGATION
%###############################################

%% Finding solutions with use of quadprog - sparse

for i = 1:1:length(sizes)
   try
       x_opt = quadprog(probs(i).H, probs(i).g, [], [], probs(i).A', probs(i).b, [], [], [], options_qp);
   catch
       x_opt = zeros(sizes(i), 1);
   end
   statsqps(i).x_opt = x_opt;
end

%% Trying all implemented algorithms and counting computation time - sparse

for i = 1:1:length(solvers)
    for j = 1:1:length(sizes) 
       try 
           tic
            [x_opt, ~] = EqualityQPSolver(probs(j).H, probs(j).g, probs(j).A, probs(j).b, solvers{i});
           cputime = toc;
       catch
           x_opt = zeros(sizes(i), 1);
           cputime = 0;
       end
       % Saving results
       statss(i).xopt(j).x_opt = x_opt;
       statss(i).cputime(j) = cputime;
    end
end

%% Detecting where minimizer wasn't found - sparse

for i = 1:1:length(sizes)
   for j = 1:1:length(solvers)
      if ~isempty(statss(j).xopt(i).x_opt) && ~isempty(statsqps(i).x_opt)
         if norm(statss(j).xopt(i).x_opt - statsqps(i).x_opt, inf) > 0.002 
            fprintf('Method %s shows weak performance for size %i \n', solvers{j}, sizes(i)); 
         end
      else
         fprintf('Method %s did not find any solution for size %i \n', solvers{j}, sizes(i)); 
      end
   end
end

%% Plotting results - sparse

% Cpu time dense problems
figure(2)
hold on;
p1 = plot(sizes, statss(1).cputime);
p2 = plot(sizes, statss(2).cputime);
p3 = plot(sizes, statss(3).cputime);
p4 = plot(sizes, statss(4).cputime);
p5 = plot(sizes, statss(5).cputime);
p6 = plot(sizes, statss(6).cputime);
title("Computation time - sparse analysis");
xlabel("Size of QP problem");
ylabel('Time (s)');
legend({'LUDense', 'LUSparse','LDLDense','LDLSparse','RangeSpace','NullSpace'}, 'location', 'best');
set(gca, "LineWidth", 2, 'FontSize', 12);
set([p1, p2, p3, p4, p5, p6], "LineWidth", 2);