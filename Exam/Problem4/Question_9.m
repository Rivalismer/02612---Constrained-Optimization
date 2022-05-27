close all, clear all, clc
% This program is investigaring all implemented methods along with
% commercial solvers to solve 3 nonlinear optimization problems, each for 3
% different starting points

%% HIMMELBLAU PROBLEM

% Initial points
x_init(1).x0 = [1; 2];
x_init(2).x0 = [-3; 0.75];
x_init(3).x0 = [5; 5];
lambda0 = [0.2; 0.2; 0.2; 0.2; 0.2; 0.2];

% loop over all methods

%% CUSTOM PROBLEM

% Initial points
x_init(1).x0 = [0.2; 0.2; 0.2];
x_init(2).x0 = [-0.7; -0.7; -0.7];
x_init(3).x0 = [1; -0.2; -0.2];
lambda0 = [0.2; 0.2; 0.2; 0.2; 0.2; 0.2; 0.2; 0.2];

%% CONSTRAINED ROSENBROCK

% Initial points
x_init(1).x0 = [0; 0];
x_init(2).x0 = [1; -0.8];
x_init(3).x0 = [-1; 2];
lambda0 = [0.2; 0.2; 0.2; 0.2; 0.2; 0.2; 0.2; 0.2];