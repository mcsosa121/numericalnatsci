function [rhos, etas, lambdas, kappas] = srp_lcurve_plot(G,d, lmin, lmax)

  if nargin < 3
    lmin = 1e-4;
    lmax = 1;
  end

  lambdas = exp(linspace(log(lmin), log(lmax)));
  rhos = 0*lambdas;
  etas = 0*lambdas;
  kappas = 0*lambdas;
  for k = 1:length(lambdas)
    [rhos(k), etas(k), kappas(k), ~] = srp_lcurve(G,d, lambdas(k));
  end
end
