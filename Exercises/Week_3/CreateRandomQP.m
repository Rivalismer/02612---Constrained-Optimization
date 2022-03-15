function [H,g,A,b] = CreateRandomQP(n_var, n_const, min_range, max_range)

H = randi([min_range max_range], n_var, n_var);
H = tril(H) + tril(H,-1);

g = randi([min_range max_range], n_var,1);
A = randi([min_range max_range], n_const, n_var);
b = randi([min_range max_range], n_const,1);