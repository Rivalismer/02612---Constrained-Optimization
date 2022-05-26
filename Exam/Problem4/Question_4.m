clear all, close all, clc
%{
    This program is creating contour plot of Himmelblau problem, with 
    saddle points and minimia location
%}
% Creating grid
x=-5:0.05:5;
[X,Y]=meshgrid(x);
Z=(X.^2+Y-11).^2+(X+Y.^2-7).^2;
v=[0:2:10 10:10:100 100:20:200];

% Creating contour plot
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
    scatter([-3.65, -3.55, -0.3, 3], [2.72, -1.42, 2.89, 2], 50, 'filled', 'red');
    scatter([-3.1, -1.4, -0.55, 0.0867, 3.25], [-0.1, 0.36, -0.22, 2.8843, 1.3], 50, 'filled', 'black');