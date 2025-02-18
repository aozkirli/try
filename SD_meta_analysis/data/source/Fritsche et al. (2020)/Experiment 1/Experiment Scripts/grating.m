function [image, log] = grating(theta, f_spat, phase, sigma, color1, color2, contrast, background, Gparams)

% Input:
% theta: orientation in degrees
% f_spat: spatial frequency in cycles/visual degrees
% phase: phase of the grating (0->1)
% sigma: SD of the gaussian window in visual degrees
% color1: first color of grating in [r g b]
% color2: second color of grating in [r g b]
% background: background color in [r g b]
% Gparams: parameters of the setup
% 
% Output:
% image: gabor image
% log: log structure containing parameters of the gabor

imSize = 8*sigma;

cycles = f_spat * imSize;

% Convert visual degree into px
imSize_px = round(deg2pix(imSize,Gparams));
sigma_px = round(deg2pix(sigma,Gparams));
f_spat_px = round(deg2pix(f_spat,Gparams));%



X = 1:imSize_px;            % X is a vector from 1 to imageSize
X0 = (X / size(X,2)) - .5;  % rescale X -> -.5 to .5


phaseRad = (phase * 2* pi);     % convert to radians: 0 -> 2*pi
[Xm Ym] = meshgrid(X0, X0);     % 2D matrices
thetaRad = (theta / 360) * 2*pi;% convert theta (orientation) to radians

Xt = Xm * cos(thetaRad);        % compute proportion of Xm for given orientation
Yt = Ym * sin(thetaRad);        % compute proportion of Ym for given orientation
XYt = [ Xt + Yt ];              % sum X and Y components
XYf = XYt * cycles * 2*pi;      % convert to radians and scale by frequency


im = sin( XYf + phaseRad);      % make 2D sinewave

im = im * contrast;

% Gaussian window
s = sigma_px / size(X,2);
gauss = exp( -(((Xm.^2)+(Ym.^2)) ./ (2* s.^2)) ); % formula for 2D gaussian


im = (im+1)./2; % rescale between 0 and 1


r = color1(1) * im + color2(1) * (1-im);
g = color1(2) * im + color2(2) * (1-im);
b = color1(3) * im + color2(3) * (1-im);


r = round(r.* gauss + background(1)*(1-gauss));
g = round(g.* gauss + background(2)*(1-gauss));
b = round(b.* gauss + background(3)*(1-gauss));


image = cat(3,r,g,b);


log = [];
log.imSize = imSize;
log.imSize_px = imSize_px;
log.theta = theta;
log.f_spat = f_spat;
log.f_spat_px = f_spat_px;
log.phase = phase;
log.sigma = sigma;
log.sigma_px = sigma_px;
log.color1 = color1;
log.color2 = color2;
log.background = background;