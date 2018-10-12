function [mtik] = srp_tik2(G, d, lambda)
    % To compare
    [U,S,V] = svd(G);
    [n,~] = size(S);

    spli = S'*S + lambda^2 * eye(n);
    sud = S' * U' * d;
    sudp = V * sud;
    mtik = spli \ sudp;
end

