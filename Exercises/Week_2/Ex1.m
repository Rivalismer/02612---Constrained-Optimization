%% Exercise 1
A = [6 2 1 -1 0; 2 5 2 0 -1; 1 2 4 -1 -1; -1 0 -1 0 0; 0 -1 -1 0 0];
b = [8 3 3 -3 0]';

x = A\b;
disp(x);

[L,D,p] = ldl(A,'lower','vector');
x(p) = L' \(D\(L\b(p)));

disp(x);
%%

f = [6 2 1; 2 5 2; 1 2 4] * [2; -1; 1] + [-8; -3; -3] - [1 0 1; 0 1 1]'*[3; -2];
disp(f);

%%
f2 = eig([6 2 1; 2 5 2; 1 2 4]);
disp(f2);

%% Exercise 2
g = [1 -2];
A = [1 0 1 1 -5; 0 1 -1 -5 1];
b = [0 0 -2 -20 -15];

x = 2:0.05:5;
y = 3:0.05:7;

[X,Y] = meshgrid(x,y);
F = X - 2*Y;
v = [-20:0.5:20];

close all
[c, h] = contour(X,Y,F,v,'linewidth', 2);
colorbar

yc1 = x;
yc2 = y;
yc3 = x + 2;
yc4 = x/5 + 4;
yc5 = 5*x - 15;
hold on
    %plot(x,yc1);
    %plot(x,yc2);
    plot(x,yc3);
    plot(x,yc4);
    plot(x,yc5);
    xlim([2 5]);
    ylim([3 7]);
    
%% Active set method

A_1 = [0 0 -1 5; 0 0 1 -1; -1 1 0 0; 5 -1 0 0];
b_1 = [-1;2;0;0];

x_1 = A_1\b_1;
disp(x_1);

%% Exercise 4
x = -4:0.05:4;
y = -4:0.05:4;

[X,Y] = meshgrid(x,y);
F = X.^2 + Y.^2 + 3*Y;
v = [-1:0.5:5];

close all
[c, h] = contour(X,Y,F,v,'linewidth', 2);
colorbar

cy1 = -1 + sqrt(x.^2 - 1);
%cy1 = X.^2 + (Y + 1).^2 -1;

hold on
    plot(x,cy1);
    
%% Exercise 5
x = -5:0.05:5;
y = -5:0.05:5;

[X,Y] = meshgrid(x,y);
F = (X.^2 + Y -11).^2 + (X + Y.^2 -7).^2;
v = [0:2:10 10:10:100 100:20:200];

close all
[c, h] = contour(X,Y,F,v,'linewidth', 2);
colorbar

yc1 = (x+2).^2;
yc2 = 4/10 * x;

hold on
plot(x, yc1);
plot(x, yc2);
ylim([-5 5]);
xlim([-5 5]);

%% 5.6
x0 = [3.5,1.5,0];

x = fsolve(@Butterfly, x0);

plot(x(1),x(2),'d');

%% 5.7
syms h_1 h_2
h = [h_1;h_2];

dc = [10, -4; -1, 10];
b = [0;0];

h = dc'\b;
disp(h)