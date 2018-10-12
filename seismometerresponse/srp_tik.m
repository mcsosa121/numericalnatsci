function [mtik] = srp_tik(G, d, lambda)
    % Solve tikhonov regularization problem
    [U,S,V] = svd(G);
    sings = diag(S);
    
    [n,~] = size(S);
    
    % filter factors
    ff = sings.^2 ./ (sings.^2 + lambda^2);
    
    mtik = 0.0;
    for k=1:n
       temp = U(:,k)' * d;
       mtik = mtik + ff(k) * temp * V(:,k); 
    end
end

