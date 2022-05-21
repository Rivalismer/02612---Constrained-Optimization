function [c, ceq] = Markowitz(lambda)

c = cumsum(lambda);
ceq = 1;