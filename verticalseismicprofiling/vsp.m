% Even determined. m = n = 80
[G, s, z, y] = vertseismicprofile( 40, 80, 80 );
% predicted data vector
t1 = G * s;
% add noise ( convert 0.1 ms to 0.0001 s )
t1 = t1 + normrnd( 0, 0.0001, [80, 1] );

% true generated slowness data 
figure
    plot(z, s );
    title('True slowness');
    xlabel('depth (m)');
    ylabel('slowness (s/m)');

% data recovered by direct differentiation
% add temp 0, and x(80) to make midpoint easier
ztemp = [0;z;z(80)];
t1temp = [0;t1;t1(80)];
sdd = zeros(80);

for i = 2:81
    sdd(i-1) = ( t1temp(i+1) - t1temp(i-1) ) / ( ztemp(i+1) - ztemp(i-1) );
end

figure
    plot(z, sdd );
    title('Slowness: Direct Diff');
    xlabel('depth (m)');
    ylabel('slowness (s/m)');

% data recovered by min length inverse 
smli = (G'*G) \ (G' * t1);
% 
figure
    plot(z, smli );
    title('Slowness: Min Length Inverse');
    xlabel('depth (m)');
    ylabel('slowness (s/m)');


% ------------------------------------------------------------------------

% Overdetermined. m = 80, n = 160
[Go, so, zo, yo] = vertseismicprofile( 40, 80, 160 );
% predicted data vector and adding noise
t2 = Go * so;
t2 = t2 + normrnd( 0, 0.0001, [160,1] );

% true generated slowness data 
figure
    plot(zo, so );
    title('Overdetermined : True slowness');
    xlabel('depth (m)');
    ylabel('slowness (s/m)');

% data recovered by direct differentiation
ztemp = [yo;yo(160)];
t1temp = [t2;t2(160)];
sdd2 = zeros(80,1);

for i = 2:2:160
    sdd2(i/2) = ( t1temp(i+1) - t1temp(i-1) ) / ( ztemp(i+1) - ztemp(i-1) );
end

figure
    plot(z, sdd2 );
    title('Overdetermined : Slowness: Direct Diff');
    xlabel('depth (m)');
    ylabel('slowness (s/m)');

% data recovered by min length inverse 
smli = (G'*G) \ (G' * t1);

figure
    plot(z, smli );
    title('Overdetermined : Slowness: Min Length Inverse');
    xlabel('depth (m)');
    ylabel('slowness (s/m)');

% comparison1
resdd = norm( sdd - s );
resmli = norm( smli - s );
fprintf('Direct diff error %d \n', resdd);
fprintf('Min Length inverse error %d \n', resmli);
% comparison2
of_resdd = norm( sdd2 - so );
of_resmli = norm( smli - so );
fprintf('Overfitting, Direct diff error %d \n', of_resdd);
fprintf('Overfitting, Min Length inverse error %d \n', of_resmli);