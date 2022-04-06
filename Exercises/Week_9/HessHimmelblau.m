function H = HessHimmelblau(x)

H = zeros(2,2);

H(1,1) = 8*x(1) + 2 + 4*(x(1).^2 + x(2) - 11);
H(1,2) = 4*(x(1) + x(2));
H(2,1) = 4*(x(1) + x(2));
H(2,2) = 8*x(2) + 2 + 4*(x(1) + x(2).^2 - 7);