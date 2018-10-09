function [lambda] = lambdacurve(G, Cd, Wm, de, d, lmin, lmax)
    fun = @(x) norm(G * wdls(G, Cd, Wm, de, x) - d, 'fro') - 0;
    
    x = linspace(lmin, lmax);
    res = zeros(length(x),1);
    
    lambda = 0;
    err = 100;
    for k=1:length(x)
       res(k) = fun(x(k));
       if ( res(k) < err )
          lambda = x(k);
          err = res(k);
       end
    end
    figure
    plot(x, res );
    title('L-curve');
    xlabel('Lambda');
    ylabel('Error');
end

