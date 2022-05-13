%% Markowitz portfolio as Quadratic Program
R = 12; % Our desired return

H = [2.50 0.93 0.62 0.74 -0.23;
    0.93 1.50 0.22 0.56 0.26;
    0.62 0.22 1.90 0.78 -0.27;
    0.74 0.56 0.78 3.60 -0.56;
    -0.23 0.26 -0.27 -0.56 3.90];


A = [16.10; 8.50; 15.70; 10.02; 18.68]';
beq = [1; R];
Aeq = [ones(5,1)'; A];
lb = zeros(5,1);
ub = ones(5,1);
f = zeros(5,1);

options = optimoptions('quadprog','Algorithm','interior-point-convex',...
    'Display', 'off', 'TolFun', 1e-10);

[x,fval,exitflag,output,lambda] = quadprog(H,f,[],[],Aeq,beq,lb,ub,[],options);
disp(x)
disp(fval) % Risk

%% Efficient Frontier
Rmin = min(A);
is = zeros((max(A) - Rmin)/0.001, 1);
fs = zeros((max(A) - Rmin)/0.001, 1);
xs = zeros((max(A) - Rmin)/0.001, 5);
cnt = 1;
for i = Rmin:0.001:max(A)
    beq = [1; i];
    [x,fval,exitflag,output,lambda] = quadprog(H,f,[],[],Aeq,beq,lb,ub,[],options);
    is(cnt) = i;
    fs(cnt) = fval;
    xs(cnt, :) = x;
    cnt = cnt + 1;
end
figure
plot(is, fs)
xlim([Rmin max(A)])
title('Efficient Frontier')
xlabel('Return')
ylabel('Risk')

figure
subplot(2,3,1)
plot(is,xs(:,1));
xlim([Rmin max(A)])
title('Asset 1')
xlabel('Return')
ylabel('Portfolio % (in decimals)')

subplot(2,3,2)
plot(is,xs(:,2));
xlim([Rmin max(A)])
title('Asset 2')
xlabel('Return')
ylabel('Portfolio % (in decimals)')

subplot(2,3,3)
plot(is,xs(:,3));
xlim([Rmin max(A)])
title('Asset 3')
xlabel('Return')
ylabel('Portfolio % (in decimals)')

subplot(2,3,4)
plot(is,xs(:,4));
xlim([Rmin max(A)])
title('Asset 4')
xlabel('Return')
ylabel('Portfolio % (in decimals)')

subplot(2,3,5)
plot(is,xs(:,5));
xlim([Rmin max(A)])
title('Asset 5')
xlabel('Return')
ylabel('Portfolio % (in decimals)')
