%% Variables
d0 = 1;
ubar = 0.2;
n = 1000;

%% Task 1.3
[H,g,A,b] = CreateQP(n,ubar,d0);

%% Task 1.4
KKT = CreateKKT(n,ubar,d0);

%% Using LU
[H,g,A,b] = CreateQP(n,ubar,d0);
KKT = CreateKKT(n,ubar,d0);
[x_lu, lambda_lu] = LUSolver(KKT,g,b,n);

%% Using LDL
[H,g,A,b] = CreateQP(n,ubar,d0);
KKT = CreateKKT(n,ubar,d0);

[x, lambda] = LDLSolver(KKT,g,b,n);

%% Using NULL-Space
[H,g,A,b] = CreateQP(n,ubar,d0);

[x_ns, lambda_ns] = NullSpaceSolver(H,A',g,b,n);

%% Using RANGE-Space
[H,g,A,b] = CreateQP(n,ubar,d0);

[x_rs, lambda_rs] = RangeSpaceSolver(H,g,A',b,n);

%% Sparsity
KKT = CreateKKT(n,ubar,d0);
spy(KKT);

%% LU with sparse
for i = 1:10
    tic
    [H,g,A,b] = CreateQP(n,ubar,d0);
    KKT = CreateKKT(n,ubar,d0);
    SLU = sparse(KKT);
    [x_lu, lambda_lu] = LUSolver(SLU,g,b,n);
    toc
end

%% LDL with sparse

for i = 1:10
    tic
    [H,g,A,b] = CreateQP(n,ubar,d0);
    KKT = CreateKKT(n,ubar,d0);
    SLDL = sparse(KKT);
    [x_lu, lambda_lu] = LDLSolver(SLDL,g,b,n);
    toc
end

%% Task 2.1

x = -4:0.01:7;
y = -4:0.01:7;

[X,Y] = meshgrid(x,y);
F = (X - 1).^2 + (Y - 2.5).^2;
v = [0:2:10 10:10:100 100:20:200];

close all
[c,h] = contour(X,Y,F,v,'linewidth',2);
colorbar

yc1 = 0.5*x + 1;
yc2 = -0.5*x + 3;
yc3 = 0*x;
xc1 = 2*y + 2;
xc2 = 0*y;

hold on
fill([x x(end) x(1)],[yc1 -4 -4], [0.7 0.7 0.7], 'facealpha', 0.2)
fill([x x(end) x(1)],[yc2 -4 -4], [0.7 0.7 0.7], 'facealpha', 0.2)
fill([x x(end) x(1)],[yc3 -4 -4], [0.7 0.7 0.7], 'facealpha', 0.2)
fill([xc1 -4 -4],[y y(end) y(1)], [0.7 0.7 0.7], 'facealpha', 0.2)
fill([xc2 -4 -4],[y y(end) y(1)], [0.7 0.7 0.7], 'facealpha', 0.2)
xlim([-4 7]);
ylim([-4 7]);
hold off