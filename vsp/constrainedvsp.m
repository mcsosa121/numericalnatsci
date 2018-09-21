% Overdetermined. m = 80, n = 160
[G, s, z, y] = vertseismicprofile( 40, 80, 160 );
% predicted data vector and adding noise
t2 = G * s;
t2 = t2 + normrnd( 0, 0.0001, [160,1] );

% true generated slowness data 
figure
    plot(z, s );
    title('Overdetermined : True slowness');
    xlabel('depth (m)');
    ylabel('slowness (s/m)');

    
% constructing constraints
b1 = 15 / (depth / m);
b2 = 25 / (depth / m);
f = eye(m,m);
f = f(:, (b1:b2));

aug = [G'*G, f; f', zeros(21,21)];

% constructing new data mat
cons = 0.002 .* ones(21,1);
t2 = [G'*t2; cons];

smli = aug \ t2;
z = [z;zeros(21,1)];

figure
    plot(z(1:m), smli(1:m));
    title('Overdetermined : Slowness: Constrained Inverse');
    xlabel('depth (m)');
    ylabel('slowness (s/m)');