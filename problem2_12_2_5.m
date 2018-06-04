
alpha = 1;
T = 0.5;

%% Part (a)
fa = @(x) sin(2*pi*x);

la = 2;
ma = la/0.4;

sola = @(x,t) exp(-4*pi^2 * t).*sin(2*pi*x);

N1 = T/0.1;

[x1,w1] = HeatForwardDifference(fa,alpha,la,T,ma,N1);

N2 = T/0.05;

[x2,w2] = HeatForwardDifference(fa,alpha,la,T,ma,N2);

%% Part (b)

fb = @(x) sin(x);
lb = pi;
mb = lb/(pi/10);

solb = @(x,t) exp(-t).*sin(x);

N3 = T/0.05;

[x3,w3] = HeatForwardDifference(fb,alpha,lb,T,mb,N3);

%% Export Data to CSVs

data12_2_5a = table(x1', w1(:,N1+1), w2(:,N2+1), sola(x2,T)', abs(w1(:,N1+1) - sola(x2,T)')./abs(sola(x2,T)'), abs(w2(:,N2+1) - sola(x2,T)')./abs(sola(x2,T)'));
writetable(data12_2_5a,'12.2.5a.csv');

data12_2_5b = table(x3', w3(:,N3+1), solb(x3,T)', abs(w3(:,N3+1) - solb(x3,T)')./abs(solb(x3,T)'));
writetable(data12_2_5b,'12.2.5b.csv');