% The psuedoinverse. If the value is 0, then still 0. Otherwise 1/sigma_p
function [sigmai] = pinv(s)
    if s == 0
       sigmai = 0;
    else 
       sigmai = 1 / s; 
    end
end