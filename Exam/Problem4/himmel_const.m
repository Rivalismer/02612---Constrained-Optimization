function [c, dc, d2c1, d2c2]=himmel_const(x)
% This function is returning constraints from himmelblau problem, as well
% as their gradient and hessians

c = [(x(1)+2)^2 - x(2); 
     -4*x(1) + 10*x(2);
     x(1)+5;
     x(2)+5;
     -x(1)+5;
     -x(2)+5]; 
    
dc = [2*x(1)+4, -4, 1, 0, -1, 0;
       -1, 10, 0, 1, 0, -1];

% Not used
d2c1 = [2  0;
        0  0];

d2c2 = [0  0;
        0  0];