function [colnorms, permutation_vector ] = swampcolumns(colnorms, permutation_vector)

% Always swapping with the first index; implies that the size of colnorms
% will reduce in loop
first = 1;

% Figure out where the maximum norm is:
[~, max_norm_location] = max(colnorms);

% Swap the first and "max_norm_location" columns:
if( max_norm_location ~= first)
    temp = colnorms(max_norm_location);
    colnorms(max_norm_location) = colnorms(first);
    colnorms(first) = temp;
    clear temp;
    
    % Apply same to the permutation vector
    temp = permutation_vector(max_norm_location);
    permutation_vector(max_norm_location) = permutation_vector(first);
    permutation_vector(first) = temp;
end

end
