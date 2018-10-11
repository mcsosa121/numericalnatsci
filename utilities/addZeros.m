function [ C ] = addZeros( C, Q )
    % Returns C padded with zeros to be same size as Q
    % Values of C along its diagonal are included
    [m,n] = size(Q);
    [x,y] = size(C);
    i = min(m,n);

    temp = zeros(m, n);
    
    if ( y == 1 )
        temp(1:i,1:i) = diag(C);
    else
        temp(1:i,1:i) = C;
    end
    C = temp;
end

