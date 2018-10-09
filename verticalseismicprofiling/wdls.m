function [mest] = wdls(G, Cd, Wm, d, eps, m0)
    % no initial guess
    if nargin < 6
        [~,m] = size(G);
        m0 = zeros(m,1);
    end
    
    cg = Cd \ G;
    A = G' * cg + eps^2 .* Wm;

    % B = G' * Cd^{-1} * (d - G*m0)
    dgm = d - G * m0;
    cdgm = Cd \ dgm;
    B = G' * cdgm;

    % m0 + A^{-1}*B
    mest = m0 + ( A \ B );
end

