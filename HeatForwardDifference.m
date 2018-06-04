%%% Heat Equation Forward-Difference

function [x,w] = HeatForwardDifference(f,alpha,l,T,m,N)
    h = l/m;
    k = T/N;
    lambda = (alpha.^2).*k./h.^2;
    
    x = h:h:(l-h);
    w = zeros(m-1,N+1);
    
    w(:,1) = f(x);
    
    mu = 1 - 2*lambda;
    A = diag(mu*ones(m-1,1), 0) + diag(lambda*ones(m-2,1), 1) + diag(lambda*ones(m-2,1), -1);

    for j=1:N
        w(:,j+1) = A*w(:,j);
    end
    
end