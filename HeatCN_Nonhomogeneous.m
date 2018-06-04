
%%% Nonhomogeneous Heat Equation Crank-Nicholson

function [x,w] = HeatCN_Nonhomogeneous(F,f,alpha,l,T,m,N)
    h = l/m;
    k = T/N;
    lambda = (alpha.^2).*k./h.^2;
    
    x = h:h:(l-h);
    w = zeros(m-1,N+1);
    
    w(:,1) = f(x);
    
    mua = 1 + lambda;
    A = diag(mua*ones(m-1,1), 0) - diag(lambda*ones(m-2,1)./2, 1) - diag(lambda*ones(m-2,1)./2, -1);
    
    mub = 1 - lambda;
    B = diag(mub*ones(m-1,1), 0) + diag(lambda*ones(m-2,1)./2, 1) + diag(lambda*ones(m-2,1)./2, -1);
    b = k*F(x);
    
    for j=1:N
        w(:,j+1) = crout(A,B*w(:,j)+b);
    end
    
end