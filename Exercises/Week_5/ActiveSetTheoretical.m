function [x, lambda, W] = ActiveSetTheoretical(H,g,A,b,x0)

% Theoretical active set solver, uses iteration to guess working sets until
% it eventually finds one in which the optimality conditions hold, at which
% point it terminates

% Start with constraint 3 and 5, as per the example
x = x0;
W = [A(:,3) A(:,5)];
done = false;
b_ws = b;
cnt = 1;
alpha = 0;
p = 0;
x = x + alpha*p;
NW = [A(:,1) A(:,2) A(:,4)];
b_nws = [b(1); b(2); b(4)];

% Main loop
while done == false
    % Try a new solution
    gk = H*x + g;
    [p, lambda] = ConvexEqSolver(H,gk,W,b_ws);
    [~,m] = size(W);
    if norm(p,2) == 0
        
        % Check optimality
        done = true;
        for i = 1:m
            if lambda(i) < 0
                done = false;
            end
        end
        
        if done == false
            % Remove something from the working set
            [~, index] = min(lambda);
            
            NW = [NW W(:,index)];
            b_nws = [b_nws; b(index)];
            
            W1 = W(:,1:(index - 1));
            W2 = W(:,(index + 1):m);
            W = [W1 W2];

            b1 = b(1:(index - 1),1);
            b2 = b((index + 1):m ,1);
            b_ws = [b1; b2];
            
        end
    else
        if size(NW) == 0
            alpha = 1;
        else
            [alpha,index] = min(-(NW'*x + b_nws)/(NW'*p));
            alpha = min(alpha);
        end
        
        if alpha < 1
            % Add to the working set
            [~,m] = size(NW);
            x = x + alpha*p;
            W = [W NW(:,index(1))];
            b_ws = [b_ws; b_nws(index(1))];
            
            NW1 = NW(:,1:(index(1) - 1));
            NW2 = NW(:,(index(1) + 1):m);
            NW = [NW1 NW2];

            b1 = b_nws(1:(index(1) - 1),1);
            b2 = b_nws((index(1) + 1):m ,1);
            b_nws = [b1; b2];
        else
            x = x + p;
        end
    end
    
    fprintf("Iteration: %d\n", cnt);
    disp(p);
    disp(lambda);
    disp(W);
    plot(x(1),x(2));
    
    cnt = cnt + 1;
end