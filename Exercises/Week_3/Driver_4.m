%% Driver for week 3

H = [6 2 1; 2 5 2; 1 2 4];
A = [1 0 1; 0 1 1];
g = [-8; -3; -3];
b = [3; 0];

[x,lambda] = EqualityQPSolver(H,g,A',b);

disp(x);
disp(lambda);

%%
[H,g,A,b] = CreateRandomQP(3,2,-20,20);

[x,lambda] = EqualityQPSolver(H,g,A',b);

disp(x);
disp(lambda);