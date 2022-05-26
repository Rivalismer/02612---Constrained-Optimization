function [c, ceq]=CustomCons(x)

c(1) = x(1)^2-x(2)^2+x(3)^2-2;
c(2) = x(1)^2+x(2)^2+x(3)^2-10;
ceq=[];