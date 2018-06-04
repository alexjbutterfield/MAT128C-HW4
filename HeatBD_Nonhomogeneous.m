%%% Nonhomogeneous Heat Equation Backwards-Difference

function [x,w] = HeatBD_Nonhomogeneous(F,f,alpha,l,T,m,N)
    h = l/m;
    k = T/N;
    lambda = (alpha.^2).*k./h.^2;
    
    x = h:h:(l-h);
    w = zeros(m-1,N+1);
    
    w(:,1) = f(x);
    
    mu = 1 + 2*lambda;
    A = diag(mu*ones(m-1,1), 0) - diag(lambda*ones(m-2,1), 1) - diag(lambda*ones(m-2,1), -1);
    b = k*F(x);
    
    for j=1:N
        w(:,j+1) = crout(A,w(:,j)+b);
    end
    
end