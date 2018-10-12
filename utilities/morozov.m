function [lambda] = morozov( func, G, d, expec, lmin, lmax )
    % Solves for optimal lambda given morozov's discrepancy principle
    % p(l) = ||Gm(l) - d|| =(approx) ||e|| where e is the error
    % func is a function computing ||Gm(1) - d||
    if ( nargin < 5 )
        lmin = -1;
        lmax = 1;
    end
    
    fun = @(x) resid(func, x, G, d) - expec;
    bounds = [lmin, lmax];
    lambda = fzero(fun, bounds);
end

