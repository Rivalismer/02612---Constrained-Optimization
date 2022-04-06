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
xc1 = 0.5*y - 1;
xc2 = 0*y;

hold on
fill([x x(end) x(1)],[yc1 -4 -4], [0.3 0.3 0.3], 'facealpha', 0.2)
fill([x x(end) x(1)],[yc2 -4 -4], [0.3 0.3 0.3], 'facealpha', 0.2)
fill([x x(end) x(1)],[yc3 -4 -4], [0.3 0.3 0.3], 'facealpha', 0.2)
fill([xc1 -4 -4],[y y(end) y(1)], [0.3 0.3 0.3], 'facealpha', 0.2)
fill([xc2 -4 -4],[y y(end) y(1)], [0.3 0.3 0.3], 'facealpha', 0.2)
xlim([-4 7]);
ylim([-4 7]);

%% 2.5 Doesn't work sadly
H = [2 0; 0 2];
g = [-2;-5];
A = [1 -2; -1 -2; -1 2; 1 0; 0 1]';
b = [-2; -6; -2; 0; 0];
x0 = [2; 0];

[x, lambda] = ActiveSetTheoretical(H,g,-A,b,x0);

%% 3.1 Markowitz portfolio
R = 10;

H = [2.30 0.93 0.62 0.74 -0.23;
    0.93 1.40 0.22 0.56 0.26;
    0.62 0.22 1.80 0.78 -0.27;
    0.74 0.56 0.78 3.40 -0.56;
    -0.23 0.26 -0.27 -0.56 2.60];


A = [15.10; 12.50; 14.70; 9.02; 17.68]';
beq = [1; R];
Aeq = [ones(5,1)'; A];
lb = zeros(5,1);
ub = ones(5,1);
f = zeros(5,1);

options = optimoptions('quadprog','Algorithm','interior-point-convex',...
    'Display', 'iter', 'TolFun', 1e-10);

[x,fval,exitflag,output,lambda] = quadprog(H,f,[],[],Aeq,beq,lb,ub,[],options);
