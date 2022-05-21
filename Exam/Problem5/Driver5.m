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
% Min return is 0, but it can't compute that as it makes the problem
% infeasible
% Max return is max(A)
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

%% Bi-criterion - minimise risk while optimising return
% Prioritisation of risk (H) versus return (f)
alpha = linspace(0,0.99,100);
beq = 1;
Aeq = ones(5,1)';
xsol = zeros(length(alpha),5);

options = optimoptions('quadprog','Algorithm','interior-point-convex',...
    'Display', 'off', 'TolFun', 1e-10);

for i=1:length(alpha)
    H_bi = (1-alpha(i))*H;
    g = -alpha(i)*A';
    [x,fval,exitflag,output,lambda] = quadprog(H_bi,g,[],[],Aeq,beq,lb,ub,[],options);
    xsol(i,:) = x;
end

figure
subplot(2,3,1)
plot(alpha,xsol(:,1));
xlim([min(alpha) max(alpha)])
title('Asset 1')
xlabel('\alpha')
ylabel('Portfolio % (in decimals)')

subplot(2,3,2)
plot(alpha,xsol(:,2));
xlim([min(alpha) max(alpha)])
title('Asset 2')
xlabel('\alpha')
ylabel('Portfolio % (in decimals)')

subplot(2,3,3)
plot(alpha,xsol(:,3));
xlim([min(alpha) max(alpha)])
title('Asset 3')
xlabel('\alpha')
ylabel('Portfolio % (in decimals)')

subplot(2,3,4)
plot(alpha,xsol(:,4));
xlim([min(alpha) max(alpha)])
title('Asset 4')
xlabel('\alpha')
ylabel('Portfolio % (in decimals)')

subplot(2,3,5)
plot(alpha,xsol(:,5));
xlim([min(alpha) max(alpha)])
title('Asset 5')
xlabel('\alpha')
ylabel('Portfolio % (in decimals)')

%% Risk free asset - 0 covariance with everything
% Also rf = 0, so we just add a 0 to the end of our A
H = [2.50 0.93 0.62 0.74 -0.23 0;
    0.93 1.50 0.22 0.56 0.26 0;
    0.62 0.22 1.90 0.78 -0.27 0;
    0.74 0.56 0.78 3.60 -0.56 0;
    -0.23 0.26 -0.27 -0.56 3.90 0
    0 0 0 0 0 0];
A = [16.10; 8.50; 15.70; 10.02; 18.68; 0]';
beq = [1; R];
Aeq = [ones(6,1)'; A];
lb = zeros(6,1);
ub = ones(6,1);
f = zeros(6,1);

[x,fval,exitflag,output,lambda] = quadprog(H,f,[],[],Aeq,beq,lb,ub,[],options);
disp(x)
disp(fval) % Risk

%% Risk free asset - Efficient frontier
Rmin = min(A);
is = zeros((max(A) - Rmin)/0.001, 1);
fs = zeros((max(A) - Rmin)/0.001, 1);
xs = zeros((max(A) - Rmin)/0.001, 6);
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

subplot(2,3,6)
plot(is, xs(:,6));
xlim([Rmin max(A)])
title('Asset 6')
xlabel('Return')
ylabel('Portfolio % (in decimals)')

%% With R = 14
Rmin = min(A);
is = zeros((max(A) - Rmin)/0.001, 1);
fs = zeros((max(A) - Rmin)/0.001, 1);
xs = zeros((max(A) - Rmin)/0.001, 6);
cnt = 1;
for i = Rmin:0.001:max(A)
    beq = [1; i];
    [x,fval,exitflag,output,lambda] = quadprog(H,f,[],[],Aeq,beq,lb,ub,[],options);
    is(cnt) = i;
    fs(cnt) = fval;
    xs(cnt, :) = x;
    cnt = cnt + 1;
end

%Compute optimal with R = 14
R = 14;
beq = [1; R];
[x,fval,exitflag,output,lambda] = quadprog(H,f,[],[],Aeq,beq,lb,ub,[],options);

% Plots
figure
hold on
plot(is, fs)
plot(14, fval, 'r.','LineWidth',2,'MarkerSize',25)
hold off
xlim([Rmin max(A)])
title('Efficient Frontier')
xlabel('Return')
ylabel('Risk')

figure
subplot(2,3,1)
hold on
plot(is,xs(:,1));
plot(14,x(1),'r.','LineWidth',2,'MarkerSize',25);
hold off
xlim([Rmin max(A)])
title('Asset 1')
xlabel('Return')
ylabel('Portfolio % (in decimals)')

subplot(2,3,2)
hold on
plot(is,xs(:,2));
plot(14,x(2),'r.','LineWidth',2,'MarkerSize',25);
hold off
xlim([Rmin max(A)])
title('Asset 2')
xlabel('Return')
ylabel('Portfolio % (in decimals)')

subplot(2,3,3)
hold on
plot(is,xs(:,3));
plot(14,x(3),'r.','LineWidth',2,'MarkerSize',25);
hold off
xlim([Rmin max(A)])
title('Asset 3')
xlabel('Return')
ylabel('Portfolio % (in decimals)')

subplot(2,3,4)
hold on
plot(is,xs(:,4));
plot(14,x(4),'r.','LineWidth',2,'MarkerSize',25);
hold off
xlim([Rmin max(A)])
title('Asset 4')
xlabel('Return')
ylabel('Portfolio % (in decimals)')

subplot(2,3,5)
hold on
plot(is,xs(:,5));
plot(14,x(5),'r.','LineWidth',2,'MarkerSize',25);
hold off
xlim([Rmin max(A)])
title('Asset 5')
xlabel('Return')
ylabel('Portfolio % (in decimals)')

subplot(2,3,6)
hold on
plot(is, xs(:,6));
plot(14,x(6),'r.','LineWidth',2,'MarkerSize',25);
hold off
xlim([Rmin max(A)])
title('Asset 6')
xlabel('Return')
ylabel('Portfolio % (in decimals)')