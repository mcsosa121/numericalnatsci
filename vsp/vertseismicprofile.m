function [G, s] = vertseismicprofile( depth, m, n )
    % model spacing and vector (m,1)
    dz = depth / m;
    z = dz .* (1:m)';
    
    % calculate slowness vector
    s = arrayfun( @(x) slowness(x), z );
    
    if ( m == n )
        % matrix is lower triangular
        G = dz .* tril( ones(n) );
    else
        % over/under/mixed determined
        % gamma vector
        dy = (m / n) * dz;
        y = dy .* (1:n)';
        
        % initializing
        G = zeros(n,m);
        for i = 1:n
           for j = 1:m
               G(i,j) = hside( y(i), z(j)) * dz;
           end
        end
    end   
end



% heaviside helper function
function[h] = hside( yi, zj )
    if ( yi - zj >= 0)
        h = 1;
    else
        h = 0;
    end
end