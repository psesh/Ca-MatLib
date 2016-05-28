% Code for computing QR factorization using Householder vectors
% Copyright (c) 2016 by Pranay Seshadri
function [Q,R] = qr_householder(A)

% Size of A
[m,n] = size(A);

for j = 1 : n
    [v,betav] = house(A(j:m,j));
    A(j:m,j:n) = (eye(m - j) - betav * (v * v') ) * A(j:m,j:n);
    if j < m
        A(j+1:m,j) = v(2:m - j + 1);
    end
end
R = triu(A); % R is the upper triangular matrix of the "new" A

% Computation of Q using backward accumulation
I = eye(m);
Q = I(:, 1:k);

end