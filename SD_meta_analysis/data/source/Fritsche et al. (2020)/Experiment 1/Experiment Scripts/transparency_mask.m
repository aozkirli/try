function [image] = transparency_mask(sigma,background,Gparams)

imSize = 14*sigma;

imSize_px = round(deg2pix(imSize,Gparams));
sigma_px = round(deg2pix(sigma,Gparams));



X = 1:imSize_px;            % X is a vector from 1 to imageSize
X0 = (X / imSize_px) - .5;  % rescale X -> -.5 to .5

[Xm Ym] = meshgrid(X0, X0);     % 2D matrices



% Gaussian window
s = sigma_px / imSize_px;

gauss = exp( -(((Xm.^2)+(Ym.^2)) ./ (2* s.^2))); % formula for 2D gaussian

% gauss(gauss < trim) = 0;

r = background(1) * ones(imSize_px);
g = background(2) * ones(imSize_px);
b = background(3) * ones(imSize_px);


a = abs(gauss-1)*255;


image = cat(3,r,g,b,a);