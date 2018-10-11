function [ U, V, Z, C, S ] = csdec( Q1, Q2 )
    % Cosine-Sine decomposition 
    % Originally implemented by David Hysell in Python
    [m,p] = size(Q1);
    [n, pb] = size(Q2);

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

    [U,C,Z] = svd(Q1, 'econ');
    % change size of C to match Q1
    C = addZeros(C,Q1);
    
    q = min(m,p);
    % flip only the diagonal
    C(1:q, 1:q) = C((q:-1:1), (q:-1:1));

    % more flips
    U(:,1:q) = U(:,(q:-1:1));
    Z(:,1:q) = Z(:,(q:-1:1));

    Z = rot90(Z,2)';
    S = Q2*Z;

    % Setting k value
    if q == 1
        k = 1;
    elseif ( m < p )
        k = n;    
    else
        % first index not less than 1/sqrt(2)
        [~, fi] = max( diag(C) <= ( 1 / sqrt(2) ) );
        k = max(1, fi );
    end
    
    % As noted, works when k=0. But test when k > 0
    [V, ~] = qr( S(:, 1:k+1) );
    
    
    S = V'*S;
    r = min(k,m);
    
    [val,~]=diagSpec(S(:,1:r));
    S(:,1:r) = val;
    
    if ( m == 1 && p > 1)
        S(0,0) = 0;
    end
    
    if ( k < min(n,p) )
        r = min(n,p);
        % get SVD
        [UT, ST, VT] = svd(S(k:n,k:r), 'econ');
        
        if ( k > 0 )
            S(1:k,k:r) = 0;
        end
        
        S(k:p, k:p) = ST;

        C(:,k:r) = C(:,k:r)*VT;
        V(:,k:n) = V(:,k:n)*UT;
        Z(:,k:r) = Z(:,k:n)*VT;
        
        i = k:q;
        t = 1:length(i);
        
        [Q,R] = qr(C(k:q,k:r));
        [val,~] = diagSpec(R);
        C(k:q,k:r) = val;
        U(:,t) = U(:,i)*Q;
    end
    
    if ( m < p )
       % Diagonalize final block and permute
       [temp1,~] = diagSpec(C,0);
       [temp2,~] = diagSpec(S,0);
       q = min( nnz( sparse( abs( temp1 > 10 * m * eps ) ) ), ...
                nnz( sparse( abs( temp2 > 10 * n * eps ) ) ) );
            
       i = q:n;
       % j = m:p;
       
       % Now S(i,j) has othogonal cols and elements of 
       % S(:,q+1:p) outside of S(i,j) are negligible
       [Q,R] = qr(S(q:n,m:p));
       S(:,q:p) = 0;
       [val,~] = diagSpec(R);
       S(q:n,m:p) = val;
       V(:,i) = V(:,i)*Q;
       
       if ( n > 1 )
          i = [q:q+p-m; 1:q; q+p-m:n];
       else
          i = 1;
       end
       
       j = [m:p; 1:m];
       t = 1:length(j);
       
       C(:,t) = C(:,j);
       S = S(q:n, m:p);
       Z = Z(:, j);
       V = V(:, i);
    end
    
    if ( n < p )
       S(:, n:p) = 0; 
    end
    
    % Make sure C, S are real and positive
    [U,C] = diagSpec( U, max(1,p-m), C );
    C = real(C);
    [V,S] = diagSpec( V, 1, S );
    S = real(S);
end

