% Solving the Even Determined Sesimic Instrument Response Inverse Problem

% A Seismometer measures displacement of weights experiencing mechanical 
% forces. However response times can be skewed and results in broadened
% signals from mechanical perturbations.  
% times from -5 to 99.5 in 0.5 increments. So accurate estimations 
% require the instrument function to be removed from the observed signal
% via deconvolution. 

% Signals received at times from -5.0 sec to 99.5 sec in 0.5 sec increments
t = ( -5.0 : 0.5 : 99.5 );
% Constants (time and gain)
tc = 10;
gc = tc / exp(1);
peaks = [8,20];
consts = [1,0.38];
width = 2;
% true ground acceleration function
ga = groundaccel( peaks, consts, width );
% true model
mt = arrayfun( @(x) ga(x), t )';

% figure
%     plot((1:length(mt)),mt);
%     title('True Ground Acceleration');
%     xlabel('Time');
%     ylabel('Response');

% Discretizing problem
G = zeros(210,210);
for i=1:length(t)
    for j=1:length(t)
       if ( i >= j )
         diff = t(i) - t(j);
         G(i,j) = diff*exp(-diff / tc ) * 0.5;
       end
    end
end

% noise free data
d = G*mt;
% figure
%     plot((1:length(t)),d);
%     title('Noise free data');
%     xlabel('Time');
%     ylabel('Voltage');

% noise is gaussian with mean 0, var = 0.01
d = d + normrnd( 0, 0.0001, [210,1] );

figure
    plot((-5:0.5:99.5),d);
    title('Noisy data');
    xlabel('Time');
    ylabel('Voltage');
