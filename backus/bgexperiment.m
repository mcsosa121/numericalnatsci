% Backus and Gilbert
Re = 6.3708*10e6;
Ie = 8.02*10e37;

g1 = @(r) 4.0 * pi * r^2;
g2 = @(r) (8.0/3.0) * pi * r^4;

h11 = @(r,rh) g1(r)*g1(r)*(r-rh)^2;
h12 = @(r,rh) g1(r)*g2(r)*(r-rh)^2;
h22 = @(r,rh) g2(r)*g2(r)*(r-rh)^2;

