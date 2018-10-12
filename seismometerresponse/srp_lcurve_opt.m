function [lambda] = srp_lcurve_opt(G, d, lmin, lmax)
  if nargin < 3
    lmin = 1e-4;
    lmax = 1;
  end
  function [invkappa] = srp_invkappa(lambda)
    [~, ~, kappa, ~] = srp_lcurve(G, d, lambda);
    invkappa = 1 / kappa;
  end
  lambda = fminbnd(@(lambda) srp_invkappa(lambda), lmin, lmax);
end

