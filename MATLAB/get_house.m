function v  = get_house(x,i,j)

n = length(x);
v = zeros(n,1);

v(i:j) = x(i:j);
v(i) = v(i) - norm(x(i:j));

if( ( v' * v) > 0 )
    v = v * sqrt(2/(v' * v));
end


end