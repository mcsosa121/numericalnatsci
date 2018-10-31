function [sol, l2prede, l2mc] = cgls(G, d, n, m0)
    %CGLS Conjugate Gradient Least Squares
    % Input: G  - Theory matrix 
    %        d  - Data vector 
    %        n  - number of iterations to perform
    %        m0 - inital guess
    % Output: sol    - predicted solution
    %         l2pred - l2 prediction error ||r||^2 = ||d-Gm||^2
    %         l2mc   - l2 model complexity length ||m||^2
    
    % setting up vectors
    pe = zeros(n,1);
    mc = zeros(n,1);
    alphas = zeros(n,1); 
    betas = zeros(n,1);
    
    % setting initial values
    s = d - G*m0;
    r = G'*s(1);
    p = r(1);
    q = G*p(1);
    
    for k = 0:n
        alphas(k+1) = (r'*r) / (q'*q);
        m = m + alphas(k+1) * p;
        
        % storing model complexity length
        mc(k+1) = norm(m)^2;
        
        % computing residual vector
        s = s - alphas(k+1) * q;
        rtemp = r;
        r = G' * s;
        
        % storing prediction error
        pe(k+1) = norm(r)^2;
        
        % computing new conjugate vectors
        betas(k+1) = (r' * r) / (rtemp' * rtemp);
        p = r + betas(k+1) * p;
        q = G * p;
    end

    sol = m;
    l2prede = pe;
    l2mc = mc;
end

