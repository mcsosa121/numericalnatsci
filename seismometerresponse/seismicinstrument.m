% Solving the Even Determined Sesimic Instrument Response Inverse Problem

% A Seismometer measures displacement of weights experiencing mechanical 
% forces. However response times can be skewed and results in broadened
% signals from mechanical perturbations.  
% times from -5 to 99.5 in 0.5 increments. So accurate estimations 
% require the instrument function to be removed from the observed signal
% via deconvolution. 

% Signals received at times from -5.0 sec to 99.5 sec in 0.5 sec increments
t = [ -5.0 : 0.5 : 99.5 ];
