% Underdetermined. m = 400, n = 100
[G, s, z, y] = vertseismicprofile( 40, 400, 100 );
% predicted data vector and adding noise
t2 = G * s;
t2 = t2 + normrnd( 0, 0.0001, [100,1] );

% taking the SVD of G
[U, Sigma, Vt] = svd(G, 'econ');
sings = diag(Sigma);
sinv = arrayfun( @(x) pinv(x), sings );
figure
    plot((1:length(sings)), sings);
    title('Singular Values of G : Underdetermined System');
    xlabel('index');
    ylabel('Singular value');

% model where all singular values retained
solreg = Vt * diag(sinv) * U' * t2;

% model where all the singular values are not retained
%half
solhlf = Vt(:,(1:50)) * diag(sinv(1:50)) * U(:,(1:50))' * t2;
% quarter
solquart = Vt(:,(1:25)) * diag(sinv(1:25)) * U(:,(1:25))' * t2;
% tenth
solten = Vt(:,(1:10)) * diag(sinv(1:10)) * U(:,(1:10))' * t2;


figure
    plot((1:length(solten)), solten);
    title('Solution using SVD psuedoinverse');
    xlabel('depth');
    ylabel('slowness');