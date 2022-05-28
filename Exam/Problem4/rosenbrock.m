function [f, df, d2f]=rosenbrock(x)
% This function is calculating function value, gradiend and Hessian of
% rosenbrock problem, given x

f = (1-x(1))^2+100*(x(2)-x(1)^2)^2;
df = [2*(200*x(1)^3-200*x(1)*x(2)+x(1)-1); 200*(x(2)-x(1)^2)];
% not used
d2f = [0];