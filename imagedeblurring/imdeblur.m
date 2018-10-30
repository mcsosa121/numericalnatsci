% load scary skeletons
im = imread('spook.jpg');
% image(im);
[m, n, ~] = size(im);
npts = m*n;

% point-spread function
sig=1.5;
wid=3;
% Add noise at 1% level
scale=0.01;
tset = zeros(1,m-wid);
a= [exp(-((0:wid-1).^2)/(2*sig^2)), zeros(1,m-wid)];
g= toeplitz(a);
g= sparse(g);
g= (1/(2*pi*sig^2))*kron(g,g);

% generate blurred image
imblur = im;
for i=1:3
    imr = double(reshape(im(:,:,i), npts, 1));
    imb = g * temp + scale * randn(npts, 1);
    imblur(:,:,i) = reshape(imb, m,n);
end

vblur = imshow(imblur);