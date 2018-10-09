% setting up G mat
sq2 = sqrt(2);
G = [ 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0; ...
      0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1; ...
      1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0; 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0; ...
      0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0; 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 1; ...
      sq2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; sq2 0 0 0 0 sq2 0 0 0 0 sq2 0 0 ...
      0 0 sq2 ];
  
[U,S,V] = svd(G);
sings = diag(S);
sinv = arrayfun( @(x) pinv(x), sings );
Sp = diag(sinv);
% G tilde
Gt = V(:,1:9) * Sp(1:9,1:9) * U(:,1:9)';

% model res Gt * G
Rm = Gt * G;

% data res 
Rd = G * Gt;

% noiseless spike test
nst = zeros(16,1);
% unity slowness and no noise
nst([6,7,10,11]) = 1;

% generate time
t2 = G*nst;

% resolve for time
sol = Gt * t2;

figure
    plot((1:length(nst)), nst);
    hold on 
    plot((1:length(sol)), sol);
    title('Original Solution verses PsuedoInverse');
    xlabel('time');
    ylabel('slowness');

% plot model and data res
% imagesc( Rm );
% ploting model nullspace
% imagesc( reshape( V(:,10) , [4,4] ) );
% imagesc( reshape( V(:,11) , [4,4] ) );
% imagesc( reshape( V(:,12) , [4,4] ) );
% imagesc( reshape( V(:,13) , [4,4] ) );
% imagesc( reshape( V(:,14) , [4,4] ) );
% imagesc( reshape( V(:,15) , [4,4] ) );
% imagesc( reshape( V(:,16) , [4,4] ) );

% ploting model rowspace
% imagesc( reshape( V(:,1) , [4,4] ) );
% imagesc( reshape( V(:,2) , [4,4] ) );
% imagesc( reshape( V(:,3) , [4,4] ) );
% imagesc( reshape( V(:,4) , [4,4] ) );
% imagesc( reshape( V(:,5) , [4,4] ) );
% imagesc( reshape( V(:,6) , [4,4] ) );
% imagesc( reshape( V(:,7) , [4,4] ) );
% imagesc( reshape( V(:,8) , [4,4] ) );
% imagesc( reshape( V(:,9) , [4,4] ) );