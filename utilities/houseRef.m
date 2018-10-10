function [H] = houseRef( v, i )
    % Householder reflection 
    % Create matrix H such that H*v has 
    % entries 0 below index i
    n = length(v);
    temp = v(i:n);
    tnorm = -norm(temp);
    if ( temp(1) < 0 )
        tnorm = -1 * tnorm;
    end
    
    if ( tnorm == 0 )
        H = eye(n);
    else
       m = length(temp);
       x = zeros(m, 1);
       
       x(1,1) = sqrt( (1/2) * (1 - ( temp(1) / tnorm ) ) );
       p = -tnorm * x(1);
       x(2:m) = temp(2:m) / (2*p);
       w = [zeros(i-1,1);x];
       H = eye(n) - 2*(w*w');
    end
end

