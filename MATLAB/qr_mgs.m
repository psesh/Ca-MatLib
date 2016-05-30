% Thin QR factorization with Modified Gram-Schmidt
% Copyright (c) 2016 by Pranay Seshadri
function [Q,R] = qr_mgs(A)
[m,n] = size(A);

for k = 1 : n
    R(k,k) =  norm(A(1:m,k),2);
    Q(1:m,k) = A(1:m,k)/R(k,k);
    for j = k + 1 :n
        R(k,j) = Q(1:m,k)' * A(1:m,j);
        A(1:m,j) = A(1:m,j) - Q(1:m,k)*R(k,j);
    end
end

end