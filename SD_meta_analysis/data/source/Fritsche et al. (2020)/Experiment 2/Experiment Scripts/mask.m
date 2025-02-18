function [im,log] = mask(sd_smooth,contrast,imSize,background,Gparams)
% MASK Generate a smoothed noisepatch that can be displayed with Psychtoolbox
%
% Input
% sd_smooth (double): s.d. of smoothing kernel for smoothing white noise
% contrast (double):  Michelson contrast from 0 to 100%
% imSize (double):    image size in visual degrees
% background (double): greyscale value of background
% Gparams (struct): structure with information about screen (size/resolution) and distance between screen and observer
%
% Output
% im (double matrix):   matrix with grayscale values of smoothed noisepatch
% log (struct):         structure with information about the generated noisepatch
%
% Matthias Fritsche - 19/05/2016

% convert from degree to pixel
imSize_px = round(deg2pix(imSize,Gparams));
sd_smooth_px  = round(deg2pix(sd_smooth,Gparams));

noisepatch = rand(imSize_px);

X = 1:imSize_px;
X0 = (X / imSize_px) - .5;  % rescale X -> -.5 to .5
[x_kernel, y_kernel] = meshgrid(X0, X0);     % 2D matrices

s = sd_smooth_px /size(X,2);

kernel = exp(-(x_kernel.^2 + y_kernel.^2)/(2*s^2))/(2*pi*s^2);

ftkernel = fft2(kernel);
ftnoise = fft2(noisepatch);
ftfilterednoise = ftkernel .* ftnoise;

filterednoise = real(ifft2(ftfilterednoise));
filterednoise = (filterednoise - min(filterednoise(:))) / range(filterednoise(:)) * 2 - 1;

% Apply contrast
filterednoise = filterednoise*contrast;
im = filterednoise * 127.5 + background;

% Log
log                 = [];
log.sd_smooth       = sd_smooth;
log.contrast        = contrast;
log.imSize          = imSize;
log.background      = background;
