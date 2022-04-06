%% 1.1 Contour plot
x = -4:0.01:7;
y = -4:0.01:7;

[X,Y] = meshgrid(x,y);
F = (X - 1).^2 + (Y - 2.5).^2;
v = [0:2:10 10:10:100 100:20:200];

close all
[c,h] = contour(X,Y,F,v,'linewidth',2);
colorbar

yc1 = 2*y - 1;
yc2 = -2*y + 6;
yc3 = 0*y;
xc2 = 0*x;
xc1 = 0.5*x - 1;

% Choose either -4 or 7 depending on if you want the fill above or below
% the line
hold on
fill([yc1 -4 -4],[y y(end) y(1)], [0.7 0.7 0.7], 'facealpha', 0.5)
fill([yc2 7 7],[y y(end) y(1)], [0.7 0.7 0.7], 'facealpha', 0.5)
fill([yc3 -4 -4],[y y(end) y(1)], [0.7 0.7 0.7], 'facealpha', 0.5)
fill([x x(end) x(1)],[xc1 -4 -4], [0.7 0.7 0.7], 'facealpha', 0.5)
fill([x x(end) x(1)],[xc2 -4 -4], [0.7 0.7 0.7], 'facealpha', 0.5)
xlim([-4 7]);
ylim([-4 7]);

%% 1.2 
A = [1 -2; -1 -2; -1 2; 1 0; 0 1]';
b = [-2; -6; -2; 0; 0];

% Find a feasible initial point - doesn't work yet
Aeq = [];
beq = [];
lb = [];
ub = [];
f = [-2 -5];
rho = 7.25;
options = optimoptions('linprog','Algorithm','interior-point','Display','iter');

[x, fval, exitflag, output, lambda] = linprog(f, A, b, Aeq, beq, lb, ub, options);

%%
% Primal acive set
H = [2 0; 0 2];
g = [-2;-5];
A = [1 -2; -1 -2; -1 2; 1 0; 0 1]';
b = [-2; -6; -2; 0; 0];
x0 = [2; 0];

[x,lambda,W] = PrimalActiveSet(H,g,A,b,x0);
disp(x);

%% Dual active set
