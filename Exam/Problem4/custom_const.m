function [c, dc, d2c]=custom_const(x)
% This function is returning constraints from custom problem, as well
% as their gradient and hessians

c = [-x(1)^2+x(2)^2-x(3)^2+2; 
     -x(1)^2-x(2)^2-x(3)^2+10;
     x(1);
     x(2)+1;
     x(3);
     -x(1)+3;
     -x(2)+3;
     -x(3)+3]; 
    
dc = [-2*x(1), -2*x(1), 1, 0, 0, -1, 0, 0;
      2*x(2), -2*x(2), 0, 1, 0, 0, -1, 0;
      -2*x(3), -2*x(3), 0, 0, 1, 0, 0, -1];

% Not used
d2c = 0;