function [f,g,H] = HimmelWithHess(x)

f = (x(1)^2 + x(2) - 11)^2 + (x(1) + x(2)^2 - 7)^2;

if nargout > 1
    g(1) = 4*x(1)*(x(1)^2 + x(2) - 11) + 2*(x(1)+x(2)^2 - 7);
    g(2) = 2*(x(1)^2 + x(2) - 11) + 4*x(2)*(x(1) + x(2)^2 - 7);
    
    if nargout > 2
        H(1,1) = 12*x(1)^2 + 4*(x(2) - 11) + 2;
        H(1,2) = 4*(x(1) + x(2));
        H(2,1) = 4*(x(1) + x(2));
        H(2,2) = 12*x(2)^2 + 4*(x(1) - 7) + 2;
    end
end