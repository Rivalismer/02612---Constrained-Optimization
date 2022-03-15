function [H,g,A,b] = CreateQP(n,ubar,d0)

H = eye(n+1);
g = zeros(n+1,1);
b = zeros(n,1);
A = zeros(n,n+1);

b(1) = -d0;
A(1,1) = -1;
A(1,n) = 1;
A(n,n-1) = 1;
A(n,n) = -1;
A(n,n+1) = -1;

for i = 2:(n-1)
    A(i,i-1) = 1;
    A(i,i) = -1;
end

for j = 1:n+1
    g(j) = ubar;
end