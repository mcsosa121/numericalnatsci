function [ U, V, X, C, S ] = gsvdHysell( A, B, ~ )
    % GSVD 
    % As originally defined by David Hysell in python
    [m,p] = size(A);
    [n,pb] = size(B);
    
    if ( p ~= pb )
       error('Matrix dimensions mismatch'); 
    end
    
    % economy sized
    if ( nargin > 2 )
        if ( m > p )
            [QA, A] = qr(A);
            [QA, A] = diagSpec( QA, A, 0 );
            m = p;
        end
        
        if ( n > p )
            [QB, B] = qr(B);
            [QB, B] = diagSpec( QB, B, 0 );
            n = p;
        end
    end
    
    [Q,R] = qr([A;B]);
    
    % Trim the cols
    Q = Q(:,1:p);
    % Trim the rows
    R = R(1:p,:);

    [U, V, Z, C, S] = csdec( Q(1:m,:), Q(m+1:m+n,:) );
    X = dot(R', Z);
    if ( exist('QA', 'var') )
        U = dot(QA, U);
    end
    
    if ( exist('QB', 'var') )
        V = dot(QB, V);
    end
end

