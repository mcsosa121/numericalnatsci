function [ G, m, t ] = seismometerResponseProb( t0, t1, dt, gc, tc )
    % A Seismometer measures displacement of weights experiencing mechanical 
    % forces. However response times can be skewed and results in broadened
    % signals from mechanical perturbations.  
    % So accurate estimations require the instrument function to be removed 
    % the observed signal via deconvolution
    %
    % Inputs:
    %   t0 - Initial time
    %   t1 - Final time
    %   dt - Time intervals (ie 0.5 s )
    %   gc - gain constant
    %   tc - time constant 
    % Outputs:
    %   G - System matrix
    %   m - the true ground acceleration model
    %   t - vector of times
    
    % Signals received at times from t0 to t1 sec in dt sec increments
    t = ( t0 : dt : t1 );
    n = length(t);
    % Constants (time and gain)
    peaks = [8,20];
    consts = [1,0.38];
    width = 2;

    % true ground acceleration function
    ga = groundaccel( peaks, consts, width );
    % true model
    m = arrayfun( @(x) ga(x), t )';

    % Discretizing problem
    G = zeros(n,n);
    for i=1:n
        for j=1:n
           diff = t(i) - t(j);
           if ( diff >= 0 )
             G(i,j) = diff*exp(-diff / tc ) * dt;
           end
        end
    end
end

