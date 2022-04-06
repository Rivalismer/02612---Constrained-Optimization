function [x,lambda, W] = PrimalActiveSet(H,g,A,b,x0)

%Find working set
W = (A'*x0 == 0);

b_ws = W;

bws = b(b_ws == 1);
bnws = b(b_ws == 0);

Ws = A(:,W == 1);
NWS = A(:,W == 0);

% Initialisations
m = size(A,1);
lambda = zeros(m,1);
x = x0;
done = false;

while ~done
    gk = H*x + g;
    
    if isempty(Ws)
        [p, lambda] = EqSolver(H,gk,[],[]);
    else
        [p, lambda] = EqSolver(H,gk,Ws, zeros(length(bws),1));
    end
    
    if norm(p,2) < 1e-4
        if sum(lambda < 0) == 0
            done = true;
        else
            [~,index] = min(lambda);
            temp = ones(length(lambda),1);
            temp(index) = 0;
            
            W(W == 1) = temp;
            b_ws(b_ws == 1) = temp;
            
            bws = b(b_ws == 1);
            bnws = b(b_ws == 0);
            Ws = A(:,W == 1);
            NWS = A(:,W == 0);
        end
    
    else
        if isempty(NWS)
            alpha = 1;
        else
            i = NWS'*p < 0;
            alph = ones(length(NWS),1);
            alph(i) = (bnws(i) - NWS(:,i)'*x)./(NWS(:,i)'*p);
            [alpha, index] = min(alph);
        end
        
        if alpha < 1
            x = x + alpha*p;
            
            temp = zeros(length(W),1);
            temp(index) = 1;
            
            W(W == 0) = temp;
            b_ws(b_ws == 0) = temp;
            
            bws = b(b_ws == 1);
            bnws = b(b_ws == 0);
            Ws = A(:,W == 1);
            NWS = A(:,W == 0);
        else
            x = x + p;
        end
    end
    disp(p);
    disp(lambda);
    disp(W);
end
