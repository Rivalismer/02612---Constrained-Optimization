function [f,df,d2f] = GradHess(fun,x,fullfun, params)

f = 0;
df = zeros(2,1);
d2f = zeros(2,2);
f = fun(x);

df(1) = diff(fullfun,params(1));
df(2) = diff(fullfun,params(2));

d2f(1,1) = diff(df(1),params(1));
d2f(1,2) = diff(df(1),params(2));
d2f(2,1) = diff(df(2),params(1));
d2f(2,2) = diff(df(2),params(2));

