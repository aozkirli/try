function [im, log] = transparency_mask(sigma,imSize,background,Gparams)
% TRANSPARENCY_MASK Generate a rgba matrix to mask stimulus images with a gaussian mask.
%
% Input
% sigma (double): standard devation of gaussian mask
% imSize (double): image size in visual degrees
% background (double): greyscale value of background
% Gparams (struct): structure with information about screen (size/resolution) and distance between screen and observer
%
% Output
% im (double matrix):   matrix with rgba values of mask image
% log (struct):         structure with information about the generated mask
%
% Matthias Fritsche - 19/05/2016

imSize_px = round(deg2pix(imSize,Gparams));
sigma_px = round(deg2pix(sigma,Gparams));


X = 1:imSize_px;            % X is a vector from 1 to imageSize
X0 = (X / imSize_px) - .5;  % rescale X -> -.5 to .5

[Xm Ym] = meshgrid(X0, X0);     % 2D matrices


% Gaussian window
s = sigma_px / imSize_px;

gauss = exp( -(((Xm.^2)+(Ym.^2)) ./ (2* s.^2))); % formula for 2D gaussian

r = background * ones(imSize_px);
g = background * ones(imSize_px);
b = background * ones(imSize_px);


a = abs(gauss-1)*255;


im = cat(3,r,g,b,a);

log                 = [];
log.sigma           = sigma;
log.imSize          = imSize;
log.background      = background;