% Even determined. m = n = 80
[G,s] = vertseismicprofile( 40, 80, 80 );
% predicted data vector
t1 = G * s;
% add noise
t1 = t1 + normrnd( 0, 0.1, [80, 1] );



% Overdetermined. m = 80, n = 160
[Go, so] = vertseismicprofile( 40, 80, 160 );
% predicted data vector and adding noise
t2 = Go * so;
t2 = t2 + normrnd( 0, 0.1, [160,1] );