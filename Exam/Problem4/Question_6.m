clear all; close all; clc;
% This program is investigating optimization of Himmelblau problem with use 
% of damped BFGS Lagrangian approximation method.

% initial guesses
x01 = [1;2];
x02 = [-3;0.75];
x03 = [5;5];
lambda0 = 0.2*ones(6, 1);

% Calculating minimizers
[x_s1, stats_s1] = SQP_BFGS(@himmelblau, @himmel_const, x01, lambda0);
[x_s2, stats_s2] = SQP_BFGS(@himmelblau, @himmel_const, x02, lambda0);
[x_s3, stats_s3] = SQP_BFGS(@himmelblau, @himmel_const, x03, lambda0);

x=-5:0.05:5;
[X,Y]=meshgrid(x);
Z=(X.^2+Y-11).^2+(X+Y.^2-7).^2;
v=[0:2:10 10:10:100 100:20:200];

% Plotting iteration sequence
[c, h]=contour(X, Y, Z, v, 'linewidth', 2);
colorbar, hold on;
yc1 = (x+2).^2;
yc2 = (4*x)/10;
hold on
    title('Himmelblau problem')
    fill(x,yc1,[0.7 0.7 0.7],'facealpha',0.4);
    fill([x x(end) x(1)],[yc2 -5 -5],[0.7 0.7 0.7],'facealpha',0.4);
    xlim([-5,5]);
    ylim([-5,5]);
    pl1 = plot(stats_s1.x(1, :), stats_s1.x(2, :), 'r', 'Linewidth', 2);
    pl2 = plot(stats_s2.x(1, :), stats_s2.x(2, :), 'm', 'Linewidth', 2);
    pl3 = plot(stats_s3.x(1, :), stats_s3.x(2, :), 'b', 'Linewidth', 2);
    legend({'', '', '', 'x_0 = [1  2]''', 'x_0 = [-3  0.75]''', 'x_0 = [5  5]'''}, 'location', 'south');
