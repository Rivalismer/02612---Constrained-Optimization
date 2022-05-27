function [f, df, d2f]=custom(x)
% This function is calculating function value, gradiend and Hessian of
% custom problem, given x

f = x(1)*x(2)+x(2)*x(3);
df = [x(2); x(1)+x(3); x(2)];
d2f = [0];