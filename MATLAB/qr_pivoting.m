function piv = qr_pivoting(A)

% Householder QR with column pivoting
[m, n] = size(A);
Q = zeros(m,n);
R = zeros(n,n);


for j = 1 : n
    c(j) = A(1:m,j)' * A(1:m, j);
end

r = 0;
tau = max(c);

while tau > 0 && r < n
    r = r + 1;
    
    % Find the smallest k with r<=k<=n so c(k) = tau!
    k_min = n;
    for g = r : n
        if(c(g) - tau == 0 && g <= k_min)
            k_min = g;
        end
    end
    k = k_min;

    
    piv(r) = k;
    A(1:m,r) = A(1:m,k);
    c(r) = c(k);
    [v,beta] = house(A(r:m,r));
    A(r:m, r:n) = (eye(m-r+1) - beta * v * v') * A(r:m,r:n) ;
    A(r+1:m,r) = v(2:m-r+1);
    for i = r + 1:n
        c(i) = c(i) - A(r,i)^2;
    end
    tau = max(c(r+1:n));
    
    disp(tau)
    disp(r)
    disp('-----');
    
end

end