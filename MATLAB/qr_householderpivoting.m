% Householder QR factorization with pivoting (Businger and Golub)
% Copyright (c) 2016 by Pranay Seshadri
function [piv, A] = qr_householderpivoting(A)
[m,n] = size(A);

% Computation of column norms
c = zeros(n);
for j = 1 : n
    c(j) = A(1:m,j)' * A(1:m,j);
end

r = 0;
tau = max(c);

while( tau > 0 && r < n)
   r = r + 1;
   k = n;
   % Find the smallest k with r<= k <= n so c(k) = tau
   for i = r : n
       if(c(i) - tau == 0 && i < k)
           k = i;
       end
   end
   piv(r) = k;
   temp = A(1:m,r);
   A(1:m,r) = A(1:m,k);
   A(1:m,k) = temp;

   temp = c(r);
   c(r) = c(k);
   c(k) = temp;
   
   [v,b] = house(A(r:m,r));
   A(r:m,r:n) = (eye(m - r + 1) - b * (v * v')) * A(1:r:m, r:n);
   A(r+1:m,r) = v(2:m - r + 1);
   for i = r + 1:n
       c(i) = c(i) - A(r,i).^2;
   end
   tau = max(c(r+1:n));
end 
  
end