function [x, lambda] = RangeSpaceSolver(H,g,A,b,n)

x = zeros(n+1,1);
lambda = zeros(n,1);

[L, info] = chol(H,'lower');

if info == 0
    v = L'\L\(g);
    Ha = A'*(L*L')*A;
    lambda = Ha\(b + A'*v);
    x = (L*L')\(A*lambda - g);
    
    lambda = -lambda;
else
    x = [];
    lambda = [];
end