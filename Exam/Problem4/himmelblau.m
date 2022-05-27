function [f, df, d2f]=himmelblau(x)
% This function is calculating function value, gradiend and Hessian of
% Himmelblau problem, given x

f = (x(1)^2+x(2)-11)^2+(x(1)+x(2)^2-7)^2;
df = [2*(2*x(1)*(x(1)^2+x(2)-11)+x(1)+x(2)^2-7); 2*(x(1)^2+2*x(2)*(x(1)+x(2)^2-7)+x(2)-11)];
d2f = [12*x(1)^2+4*x(2)-42, 4*x(1)+4*x(2); 4*x(1)+4*x(2), 4*x(1)+12*x(2)^2-26];