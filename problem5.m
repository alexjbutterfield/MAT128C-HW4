%%% Computer project 1 for HW 4

% Define problem parameters
F = @(x) 0;
f = @(x) 2*x.*(x<=0.5) + 2*(1-x).*(x>0.5);
alpha = 1;
l = 1;
T = 0.5;
m = l/(1/10);

% define each N for each part
Na = T/(1/1000);
Nb = T/(5/1000);
Nc = T/(1/100);

% consolodate N's into a vector for easy repeated accesses
Ns = [Na Nb Nc];

% define Kth sum element for true solution
sol = @(x,t,K) 8*(sin(K .*pi./2).*(sin(K .*pi.*x)).*exp(-K.^2 .* pi.^2 .* t)./(K.^2))./pi^2;

% solve problems using our codes
for i = 1:3
    [~,w1] = HeatForwardDifference(f,alpha,l,T,m,Ns(i));
    [~,w2] = HeatBD_Nonhomogeneous(F,f,alpha,l,T,m,Ns(i));
    [x,w3] = HeatCN_Nonhomogeneous(F,f,alpha,l,T,m,Ns(i));
    
    wFD(:,i) = w1(:,Ns(i)+1);
    wBD(:,i) = w2(:,Ns(i)+1);
    wCN(:,i) = w3(:,Ns(i)+1);
end

% transpose the x's for table dispaly
x = x';

% generate 100 sum approximation of true solution u(x,t)
u100 = zeros(size(x));
for i = 1:100
    u100 = sol(x,T,i) + u100;
end

% consolodate data to table and export as csv file
data5 = table(x,wFD,wBD,wCN,u100,abs(wFD - repmat(u100,1,3))./repmat(u100,1,3),abs(wBD - repmat(u100,1,3))./repmat(u100,1,3),abs(wCN - repmat(u100,1,3))./repmat(u100,1,3));
writetable(data5, 'problem5.csv');