function [ Y, X ] = diagSpec( X, k, Y )
    if ( nargin == 1 )
        % Diag force
        type = 'F';
    elseif ( nargin == 2 )
        % Kth matrix Diag
        type = 'K';
    else 
        % Diagonal Positive
        type = 'P';
    end

    if ( type == 'K' )
        % Kth matrix diagonal
        % The k-th diagonal of X, even if X a vector
        [m,n] = size(X);
        if ( min(m,n) > 1 )
           Y = diag(X, k);
        elseif ( 0 < k && 1 + k < n )
           Y = X(1+k); 
        elseif ( k < 1 && 2 - k < m )
           Y = X(1-k);
        else
           Y = []; 
        end
    elseif ( type == 'F' )
        % Diagonal Force
        % Zeros all the elements off the main diag of X
        Y = triu(tril(X));
    elseif ( type == 'P' )
        % Diagonal positive
        % Scales the columns of Y and Rows of X by unimodular factors
        % to make kth diagonal of X real and positive
        [D,~] = diagSpec(Y,k);
        [d1, d2] = size(D);
        
        for i=1:d2
            if( ~isreal(D(i,:)) )
                D = diag( conj(D(i,:)), abs(D(i,:)) );
                Y(:,i) = dot( Y(:,i), D' );
                X(i,:) = dot( D, X(i,:) );
            end
        end 
    end
end

