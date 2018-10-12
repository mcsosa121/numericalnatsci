function [r] = resid( func, lm, G, d)
    % Compute residual for given problem
    % func - function computing inverse solution
    % lm either model or lambda 
    
    % using lambda 
    if ( length(lm) == 1 )
        m = func(G,d,l);
    else
        m = lm;
    end
    
    r = norm( G*m - d );
end

