function [rnorm, xnorm, kappa, deta] = srp_lcurve(G, d,lambda)
    % Lcurve for srp problem 
    mtik = srp_tik(G,d,lambda);
    resid = G*mtik - d;
    z = srp_tik(G,resid,lambda);
    
    eta = mtik'*mtik;
    rho = resid'*resid;
    deta = 4/lambda * (d' * z);

    kappa = -2*eta*rho/deta * ...
            (lambda^2*deta*rho + 2*lambda*eta*rho + lambda^4*eta*deta) / ...
            (lambda^4*eta^2 + rho^2)^(1.5);

    xnorm = sqrt(eta);
    rnorm = sqrt(rho);
end

