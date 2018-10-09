function [s] = slowness( x )
    if ( x < 15 || x > 25 )
       s = 7.5e-5 * x + 2e-3;
    else
       s = 0.002;
    end
    
end