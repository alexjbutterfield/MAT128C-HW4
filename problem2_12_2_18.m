F = @(x) 2;
f = @(x) sin(pi*x) + x.*(1-x);
alpha = 1;
l = 1;
T = 0.25;
m = 10;
N = T/0.01;

[xbd,wbd] = HeatBD_Nonhomogeneous(F,f,alpha,l,T,m,N);
[xcn,wcn] = HeatCN_Nonhomogeneous(F,f,alpha,l,T,m,N);

sol = @(x,t) exp(-pi^2 * t).*sin(pi*x) + x.*(1-x);

data12_2_18 = table(xcn', wbd(:,N+1), wcn(:,N+1),sol(xcn',T),abs(wbd(:,N+1) - sol(xcn',T))./abs(sol(xcn',T)), abs(wcn(:,N+1) - sol(xcn',T))./abs(sol(xcn',T)));
writetable(data12_2_18, '12.2.18.csv');