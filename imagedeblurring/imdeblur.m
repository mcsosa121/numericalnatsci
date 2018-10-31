% load scary skeletons
im = imread('spook.jpg');
% colormap(gray(256));
% by converting to greyscale no longer need 3 dim
% https://www.mathworks.com/help/matlab/ref/rgb2gray.html
im = .2989*im(:,:,1) + .5870*im(:,:,2) + .1140*im(:,:,3);
im = im2double(im);
[m, n, ~] = size(im);
npts = m*n;

% point-spread function
sig=1.5;
wid=3;
% Add noise at 1% level
scale=0.01;
a= [exp(-((0:wid-1).^2)/(2*sig^2)), zeros(1,m-wid)];
g= toeplitz(a);
g= sparse(g);
g= (1/(2*pi*sig^2))*kron(g,g);

% generate blurred image
imr = reshape(im, npts, 1);
imb = g * imr + scale * randn(npts, 1);

% imblur = reshape(imb, m,n);
% imshow(imblur);
m0 = zeros(npts,1);

[mopt, m2, e2] = cgls(g, imb, 100, m0);
figure
    plot(e2, m2 );
    title('L-curve');
    xlabel('L2 error');
    ylabel('L2 length');