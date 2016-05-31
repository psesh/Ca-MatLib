function [permutation_vector] = qr_householderpivoting(A)

% NEED TO FIGURE OUT WHAT TO DO IF A IS "FAT"
%R = A;
[m,n] = size(A);
k = min(m,n);

% Compute the column norms!
colnorms = zeros(n,1);
for i = 1 : n
    colnorms(i) = A(:,i)' * A(:,i);
end

% Initialize the permutation vector!
permutation_vector = 1 : 1 : n;

% Swaping loop
for current_index = 1 : k-1
    
    [c2,p2] = swampcolumns(colnorms(current_index:n), permutation_vector(current_index:n) );
    colnorms = [colnorms(1:current_index-1); c2];
    permutation_vector = [permutation_vector(1: current_index-1), p2];
    
    % Now compute the householder vector of A(:,1)
%     R(current_index:m,current_index);
%     [v,b] = house(R(current_index:m,current_index));
%     
%     % Apply the transformation to R from the left
%     R(current_index:m,current_index:n) = (eye(m-current_index+1) - b * (v * v')) * R(current_index:m,current_index:n);
%     
    [v,betav] = house(A(current_index:m,current_index));
    A(current_index:m,current_index:n) = (eye(m-current_index+1) - betav * (v * v') ) * A(current_index:m,current_index:n);
    if current_index < m
        A(current_index+1:m,current_index) = v(2:m - current_index + 1);
    end
    
    for i = current_index + 1 : n
        colnorms(i) = colnorms(i) - A(current_index+1,i)^2;
    end
    %permutation_vector
end

% % Computation of Q using backward accumulation
% k = min(m,n);
% Q = eye(m,m);
% for j = k : -1 : 1
%     v = [1; R(j+1:m,j)];
%     betav = 2/(1 + norm(R(j+1:m,j), 2)^2); % We get the beta's from the stored Householder vectors!
%     Q(j:m,j:m) = Q(j:m,j:m) - (betav * (v*v') * Q(j:m,j:m));
% end


end

