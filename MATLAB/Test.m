clear; clc;

% Testing QR Householder with Pivoting!
A = rand(6,4); R = A; Q = eye(6);
[m,n] = size(A);
colnorms = zeros(n,1);
for i = 1 : n
    colnorms(i) = R(:,i)' * R(:,i);
end
permutation_vector = 1 : 1 : n;

%---------------
% Current index 1
%---------------
current_index = 1;
[colnorms, permutation_vector ] = swampcolumns(colnorms, permutation_vector);

% Now compute the householder vector of A(:,1)
[v,b] = house(R(current_index:m,current_index));

% Apply the transformation to R from the left
R(1:m,1:n) = (eye(m-current_index+1) - b * v* v') * R(1:m,1:n)


%---------------
% Current index 2
%---------------
current_index  = current_index + 1; % Here its "2"
for i = current_index : n
    colnorms(i) = colnorms(i) - R(1,i)^2;
end

% Now once again we have to swap. Exchange the second column with the max
% norm location!
[c2, p2] = swampcolumns(colnorms(current_index:n), permutation_vector(current_index:n));
colnorms = [colnorms(1:current_index-1); c2]; 
permutation_vector = [permutation_vector(1: current_index-1), p2];

% Compute the householder vector of A(:,2)
[v,b] = house(R(current_index:m,current_index));
R(current_index:m,current_index:n) = (eye(m-current_index+1) - b * (v * v')) * R(current_index:m,current_index:n)

%---------------
% Current index 3
%---------------
current_index = current_index + 1;
for i = current_index : n
    colnorms(i) = colnorms(i) - R(1,i)^2;
end

[c3,p3] = swampcolumns(colnorms(current_index:n), permutation_vector(current_index:n) );
colnorms = [colnorms(1:current_index - 1); c3];
permutation_vector = [permutation_vector(1: current_index - 1), p3];
[v,b] = house(R(current_index:m, current_index) );
R(current_index:m,current_index:n) = (eye(m-current_index+1) - b * (v * v')) * R(current_index:m,current_index:n)


%---------------
% Current index 4
%---------------
current_index = current_index + 1;
for i = current_index : n
    colnorms(i) = colnorms(i) - R(1,i)^2;
end

[c4,p4] = swampcolumns(colnorms(current_index:n), permutation_vector(current_index:n) );
colnorms = [colnorms(1:current_index - 1); c4];
permutation_vector = [permutation_vector(1: current_index - 1), p4];
[v,b] = house(R(current_index:m, current_index) );
R(current_index:m,current_index:n) = (eye(m-current_index+1) - b * (v * v')) * R(current_index:m,current_index:n)

