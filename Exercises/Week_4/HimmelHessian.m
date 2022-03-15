function H = HimmelHessian(x, lambda)

H(1,1) = 12*x(1)^2 + 4*(x(2) - 11) + 2 + 2*lambda.ineqnonlin(1);
H(1,2) = 4*(x(1) + x(2));
H(2,1) = 4*(x(1) + x(2));
H(2,2) = 12*x(2)^2 + 4*(x(1) - 7) + 2;