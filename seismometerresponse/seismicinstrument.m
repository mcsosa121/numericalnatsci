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
dnoisy = d + normrnd( 0, v, [210,1] );

% figure
%     plot(t,dnoisy);
%     xlim([-5 100]);
%     ylim([0 max(dnoisey)]);
%     title('Noisy data');
%     xlabel('Time');
%     ylabel('Voltage');

% solving using singular value decomp
[U, Sigma, V, m] = srp_svd(G, d);
% noisy solution
[Un, Sigman, Vn, mn] = srp_svd(G,dnoisy);

% figure
%     plot(t,m);
%     xlim([-5 100]);
%     ylim([0 1]);
%     title('Noise free SVD model');
%     xlabel('Time');
%     ylabel('Response');
% 
% figure
%     plot(t,mn);
%     xlim([-5 100]);
%     ylim([min(mn) max(mn)]);
%     title('Noisy SVD model');
%     xlabel('Time');
%     ylabel('Response');

sings = diag(Sigma);
ind = find(sings<10e-5);
sings = sings(1:ind-1);
U = Un(:,1:length(sings));
V = Vn(:,1:length(sings));
% figure
%     clf
%     semilogy(sings);
%     axis tight;
%     title('Singular values');

condK = max(sings) / min(sings);
disp(condK);
disp('Yikes that is big');

% Using truncated SVD, 25 values
sings = sings(1:25);
Ut = U(:,1:25);
Vt = V(:,1:25);
Gp = Ut * diag(sings) * Vt';
sinv = arrayfun( @(x) pinv(x), sings );
sol_trunc = Vt * diag(sinv) * Ut' * dnoisy;

% figure
%     plot(t,sol_trunc);
%     xlim([-5 100]);
%     ylim([min(sol_trunc) max(sol_trunc)]);
%     title('TSVD model');
%     xlabel('Time');
%     ylabel('Response');

% Model resolution matrix Rm = Gtilde*G = Vt*Vt'
Rm = Vt*Vt';

% figure
%     imagesc(Rm);
%     title('Truncated Model Resolution Matrix');

% Seeing how well the model approximates the identity
cut = Rm(:,100);

% figure
%     plot((1:210), cut);
%     xlabel('Time');
%     ylabel('Val');
%     xlim([0 210]);
%     ylim([min(cut) max(cut)]);
%     title('Cut through Rm');

% Now plotting the L-curve
[rhos, etas, lambdas, kappas] = srp_lcurve_plot(G,d);
rhos = rhos';
etas = flip(etas)';
rho_corner = srp_lcurve_opt(G,d, 1e-4, 1);
optlambda = srp_lcurve_opt(G,d,1e-4,1);
[ro, eo, ~, ~] = srp_lcurve(G,d,optlambda);
figure(1)
    clf
    plot( lambdas, rhos, 'o'); 
    hold on;
    axis square;
    %plot(rhos,lambdas, etas, lambdas);
%     loglog(ro,eo);
    title('Lcurve');
    xlabel('||Gm-d||');
    ylabel('||m||');
    
% figure(2)
% clf
% bookfonts;
% plot(rhos,etas,'k-');
% xlabel('Residual Norm ||Gm - d||_2')
% ylabel('Solution Norm ||m||_2')
% axis tight
% hold on
% H=plot(ro,eo,'ko');
% set(H,'markersize',16)
% hold off
% print -deps2 c4flcurve0
% display('Displaying the L-curve (fig. 1)')
