function J = HimmelJac(x)

J = zeros(2,1);
J(1) = 2*(x(1) + 2);
J(2) = -1;