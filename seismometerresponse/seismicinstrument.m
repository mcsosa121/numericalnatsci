% Solving the Even Determined Sesimic Instrument Response Inverse Problem

% Constants (time and gain)
tc = 12;
gc = tc / exp(1);

% Model
[G, mtrue, t] = seismometerResponseProb( -5.0, 99.5, 0.5, gc, tc );

% True Signal Figure
% figure
%     plot(t,mtrue);
%     xlim([-5 100]);
%     ylim([0 1]);
%     title('True Ground Acceleration');
%     xlabel('Time');
%     ylabel('Response');

% noise free data
d = G*mtrue;

% figure
%     plot(t,d);
%     xlim([-5 100]);
%     ylim([0 max(d)]);
%     title('Noise free data');
%     xlabel('Time');
%     ylabel('Voltage');

% noise is gaussian with mean 0, var = 1% of peak output
v = 0.01 * max(d);
dnoisey = d + normrnd( 0, v, [210,1] );

figure
    plot(t,dnoisey);
    xlim([-5 100]);
    ylim([0 max(dnoisey)]);
    title('Noisy data');
    xlabel('Time');
    ylabel('Voltage');
