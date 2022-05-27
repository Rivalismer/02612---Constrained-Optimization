function [c, ceq]=RosenCons(x)

c(1) = -(x(1)+2)^2+x(2)-1;
c(2) = (x(1)+2)^2-x(2)-2;
ceq=[];