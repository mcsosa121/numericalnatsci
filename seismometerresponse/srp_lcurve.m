function [rnorm, xnorm, kappa, deta] = srp_lcurve(G, d,lambda)
    % Lcurve for srp problem 
    % Rnorm is ||Gx - b||
    % xnorm is ||x||
    % kappa, deta 
    mtik = srp_tik2(G,d,lambda);
    resid = G*mtik - d;
    z = srp_tik2(G,resid,lambda);
    
    % Note: eta = ||mtik||_2
    %       rho = ||G*mtik - d||_2 = ||resid||_2
    eta = mtik'*mtik;
    rho = resid'*resid;
    
    % derivative of eta with respect to lambda
    deta = 4/lambda * (d' * z);

    % curvature of the L-curve
    c = - 2 * (eta * rho ) / deta;
    numer = lambda^2*deta*rho + 2*lambda*eta*rho + lambda^4*eta*deta;
    denom = lambda^2 + eta^2 + rho^2;
    kappa = c*numer / denom^(1.5);

    xnorm = eta;
    rnorm = rho;
end

