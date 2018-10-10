function [ V, U, Z, S, C ] = csdec( Q1, Q2 )
    % Cosine-Sine decomposition 
    % Originally implemented by David Hysell in Python
    m,p = size(Q1);
    n, pb = size(Q2);
    
    if ( pb ~= p )
        error('Matrix dimensions must match');
    end
    
    if ( m < n )
        % Recurse
       [ V, U, Z, S, C ] = csdec( Q2, Q1 );
       % flip columns L-R to R-L
       C = fliplr(C);
       S = fliplr(S);
       Z = fliplr(Z);
       
       m = min(m,p);
       % flip some more
       i = (m:-1:2);
       C(1:m,:) = C(i,:);
       U(:,1:m) = U(:,i);
       
       % and some more
       n = min(n,p);
       i = (n:-1:1);
       S(1:n,:) = S(i,:);
       V(:,1:n) = V(:,i);
       
       return;
    end
    
    [U,C,Z] = svd(Q1);
    C = add_zeros(C,Q1);
        


end

