function [f,g] = HimmelWithGrad(x)

f = (x(1)^2 + x(2) - 11)^2 + (x(1) + x(2)^2 - 7)^2;

if nargout > 1
    g(1) = 4*x(1)*(x(1)^2 + x(2) - 11) + 2*(x(1)+x(2)^2 - 7);
    g(2) = 2*(x(1)^2 + x(2) - 11) + 4*x(2)*(x(1) + x(2)^2 - 7);
end