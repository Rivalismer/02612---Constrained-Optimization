function H = ConstraintHess(x)

H = zeros(2,2);
H(1,1) = 2;
H(1,2) = 0;
H(2,1) = 0;
H(2,2) = 0;