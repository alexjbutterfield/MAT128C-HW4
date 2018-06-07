
%% Part (a)

fa=@(x,y) -sin(pi*x.*y); % define input function data 
ra=@(x,y) -(x.^2 + y.^2);
xl = 0;
xr = 1;
yb = 0;
yt = 1;
M = 10;
N = 10;
figure(1)
wa = dirchlet_poissonfem(fa, ra, xl,xr,yb,yt,M,N);

fb = @(x,y) exp(2*x.*y);
rb = @(x,y) sin(pi*x.*y);

figure(2)
wb = dirchlet_poissonfem(fb, rb, xl,xr,yb,yt,M,N);