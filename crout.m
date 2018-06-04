

function x = crout(A,b)
    n = size(A,1);
    
    l = zeros(n);
    u = zeros(n);
    z = zeros(n,1);
    x = zeros(n,1);
    
    l(1,1) = A(1,1);
    u(1,2) = A(1,2)/l(1,1);
    z(1) = b(1)/l(1,1);
    
    for i=2:(n-1)
        l(i,i-1) = A(i,i-1);
        l(i,i) = A(i,i) - l(i,i-1)*u(i-1,i);
        u(i,i+1) = A(i,i+1)/l(i,i);
        z(i) = (b(i) - l(i,i-1)*z(i-1))/l(i,i);
    end
    
    l(n,n-1) = A(n,n-1);
    l(n,n) = A(n,n) - l(n,n-1)*u(n-1,n);
    z(n) = (b(n) - l(n,n-1)*z(n-1))/l(n,n);
    
    x(n) = z(n);
    
    for i = (n-1):-1:1
        x(i) = z(i) - u(i,i+1)*x(i+1);
    end
    
end