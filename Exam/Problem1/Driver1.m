%% Problem 1.4 - Test with given values
H = [5, 1.86, 1.24, 1.48, -0.46
    1.86, 3, 0.44, 1.12, 0.52
    1.24, 0.44, 3.8, 1.56, -0.54
    1.48, 1.12, 1.56, 7.2, -1.12
    -0.46, 0.52, -0.54, -1.12, 0.78];

g = [-16.1
    -8.5
    -15.7
    -10.02
    -18.68];

A = [16.1, 1
    8.5, 1
    15.7, 1
    10.02, 1
    18.68, 1];

b = [15
    1];

[x, lambda] = EqualityQPSolver(H,g,A,b,'NullSpace');