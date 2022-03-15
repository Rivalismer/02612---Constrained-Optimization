%% Contour plot
x = -10:0.05:10;
y = -10:0.05:10;
[X,Y] = meshgrid(x,y);
F = X.^2 - 2*X + 3*X*Y + 4*Y.^3;

v = [25000:500:600000];

[c,h] = contour(X,Y,F,v,'linewidth',2);
colorbar

hold on

%% Gradient and Hessian
syms X Y
F = X.^2 - 2*X + 3*X*Y + 4*Y.^3;
param = [X; Y];
x = [2; 3];
[y,df,d2f] = GradHess(@Exercise1, x, F, param);

disp([y,df,d2f]);