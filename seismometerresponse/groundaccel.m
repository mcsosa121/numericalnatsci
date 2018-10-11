function [g] = groundaccel(peaks, consts, width )
    % True model for ground acceleration
    % peaks are [p1,p2,...] how many peaks
    % consts are [c1, c2,...] how high (percentage of max response)
    % width w, how wide each peak is (this is not a vector)
    m = length(peaks);
    n = length(consts);
    
    if ( m ~= n )
       error('Array lengths must match'); 
    end
    
    g = @(t) sum( arrayfun( @(x,y) x*exp(-(t - y)^2 / (2 * width^2 ) ), consts, peaks ) );
end

