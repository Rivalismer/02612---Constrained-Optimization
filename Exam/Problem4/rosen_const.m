function [c, dc, d2c]=rosen_const(x)
% This function is returning constraints from rosenbrock problem, as well
% as their gradient and hessians

c = [(x(1)+2)^2 - x(2)+1;
     -(x(1)+2)^2 + x(2)+2;
     -4*x(1) + 10*x(2)+1;
     4*x(1) - 10*x(2)+2;
     x(1)+1;
     x(2)+1;
     -x(1)+2;
     -x(2)+2]; 
    
dc = [2*x(1)+4, -2*x(1)-4, -4, 4, 1, 0, -1, 0;
      -1, 1, 10, -10, 0, 1, 0, -1];

% Not used
d2c = 0;