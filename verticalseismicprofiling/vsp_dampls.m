% Solve the Vertical Seismic problem using damped least squares
% aka G*= m0 + (G'*Cd^{-1}*G + eps^2*Wm)^{-1}*G'*Cd^{-1}*(d - G*m0)
% Wm = L'L  where L is second derivative damping matrix

% m = n = 80
% Cd = I
depth = 40;
m = 80;
n = 80;
dz = depth / n;
Cd = eye(n);

% creating the weight matrix (2nd derivative)
c = zeros(m,1);
c(1) = -2;
c(2) = 1;
c(end)=1;
L = (1/dz^2) .* toeplitz(c, c');

%for weighting
b1 = round(15 / (depth / m));
b2 = round(25 / (depth / m));
wn = zeros(1,m);
wn(b1:b2) = 0.5;
L = L + repmat(wn, m, 1);
Wm = L'*L;

% setting up model
[G, s, z, y] = vertseismicprofile( depth, m, n);

% no initial guess/apriori estimate so eq is now
% s* = (G'*Cd^-1*G + eps^2*Wm)^-1 * G'*Cd^-1 * d
% predicted data vector
t1 = G * s;
% add noise ( convert 0.1 ms to 0.0001 s )
t1 = t1 + normrnd( 0, 0.0001, [80, 1] );


mest = wdls( G, Cd, Wm, t1, 0.001);

figure
    plot(z, mest );
    title('Slowness: WDLS, Small lambda');
    xlabel('depth (m)');
    ylabel('slowness (s/m)');



% with initial guess of all 2's
