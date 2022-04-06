function [x, lambda] = DualActiveSet(H,g,A,b)

x = H\g;
lambda = 0;
W = zeros(size(b),1);
WS = A(W == 1);
NWS = A(W == 0);
bws = b(W == 1);
bnws = b(W == 0);
done = false;

while ~done
    i = (WS*x - bws > 0);
    
    
end