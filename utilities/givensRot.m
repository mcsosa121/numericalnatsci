function [c, s, r] = givensRot( a, b )
    % Givens rotation 
    % Rotates [a;b] to [r;0] via [c -s; s c]*[a;b]
    if ( b == 0 )
        c = sign(a);
        s = 0.0;
        r = abs(a);
        return;
    end
    
    if ( a == 0 )
        c = 0.0;
        s = -sign(b);
        r = abs(b);
        return;
    end
    
    if ( abs(a) > abs(b) )
        x = b / a;
        y = sign(a) * abs( sqrt( 1 + x^2 ) );
        c = 1 / y;
        s = -c * x;
        r = a * y;
    else
       x = a / b;
       y = sign(b) * abs( sqrt( 1 + x^2 ) );
       s = -1 / y;
       c = - s * x;
       r = b * y;
    end
end

