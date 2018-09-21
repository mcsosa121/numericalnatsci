% Solve the Vertical Seismic problem using damped least squares
% aka G*= m0 + (G'*Cd^{-1}*G + eps^2*Wm)^{-1}*G'*Cd^{-1}*(d - G*m0)
% Wm = L'L  where L is second derivative damping matrix

% m = n = 80
% Cd = I
depth = 40;
m = 80;
n = 10;
dz = depth / n;
Cd = eye(n);

% creating the weight matrix (2nd derivative)
c = zeros(m,1);
c(1) = -2;
c(2) = 1;
c(end)=1;
L = (1/dz^2) .* toeplitz(c, c');
Wm = L'*L;

% setting up model
[G, s, z, y] = vertseismicprofile( depth, m, n);

% no initial guess/apriori estimate so eq is now
% s* = (G'*Cd^-1*G + eps^2*Wm)^-1 * G'*Cd^-1 * d
% predicted data vector
t1 = G * s;
% add noise ( convert 0.1 ms to 0.0001 s )
t1 = t1 + normrnd( 0, 0.0001, [n, 1] );


aopt = morozov( G, Cd, Wm, t1, G*s, 0, 1);
mest = wdls( G, Cd, Wm, t1, aopt);

figure
    plot(z, mest );
    title('Slowness: WDLS, optimal lambda');
    xlabel('depth (m)');
    ylabel('slowness (s/m)');
