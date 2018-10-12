function [U, Sigma, V, m] = srp_svd( G, d )
    % Solves Seismometer (instrument) response problem using SVD
    % Input:
    %   G - Model
    %   d - data
    % Output:
    %   U, Sigma, V' as computed in SVD
    %   Solution m*
    [U, Sigma, V] = svd(G);
    sings = diag(Sigma);
    sinv = arrayfun( @(x) pinv(x), sings );
    m = U' * d;
    m = diag(sinv) * m;
    m = V * m;
end

