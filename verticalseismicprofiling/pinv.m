% The psuedoinverse. If the value is 0, then still 0. Otherwise 1/sigma_p
function [sigmai] = pinv(s)
    if s > 1e-7
       sigmai = 1 / s;
    else 
       sigmai = 0; 
    end
end