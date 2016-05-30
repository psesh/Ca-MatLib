function
R = A; 
[m,n] = size(A);
colnorms = zeros(n,1);
for i = 1 : n
    colnorms(i) = R(:,i)' * R(:,i);
end
permutation_vector = 1 : 1 : n;

for current_index = 1 : n
    
    [c2,p2] = swampcolumns(colnorms(current_index:n), permutation_vector(current_index:n) );
    colnorms = [colnorms(1:current_index-1); c2];
    permutation_vector = [permutation_vector(1: current_index-1), p2];
    
    % Now compute the householder vector of A(:,1)
    [v,b] = house(R(current_index:m,current_index));
    
    % Apply the transformation to R from the left
    R(current_index:m,current_index:n) = (eye(m-current_index+1) - b * (v * v')) * R(current_index:m,current_index:n);
    
    for i = current_index + 1 : n
        colnorms(i) = colnorms(i) - R(current_index,i)^2;
    end
    
end

