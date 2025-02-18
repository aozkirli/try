function [im, log] = grating(theta, f_spat, phase, contrast, imSize, background, Gparams)
% GRATING Generate a sinusoidal grating that can be displayed with Psychtoolbox.
%
% Input
% theta (double):     grating angle in (deg) - zero angle: north; rotation direction: cw
% f_spat (double):    grating spatial frequency (cycles/visual degree)
% phase (double):     phase from 0 to 1
% contrast (double):  Michelson contrast from 0 to 100%
% imSize (double):    image size in visual degrees
% background (double): greyscale value of background
% Gparams (struct): structure with information about screen (size/resolution) and distance between screen and observer
%
% Output
% im (double matrix):   matrix with grayscale values of grating image
% log (struct):         structure with information about the generated grating
%
% Matthias Fritsche - 19/05/2016


%% Create 2D sinewave image matrix
cycles = f_spat * imSize; % compute number of cycles for sinusoidal grating

% Convert visual degree into px
imSize_px           = round(deg2pix(imSize,Gparams));

X = 1:imSize_px;                    % X is a vector from 1 to imageSize
X0 = (X / size(X,2)) - .5;          % rescale X -> -.5 to .5
[Xm, Ym] = meshgrid(X0, X0);        % create image matrix

phaseRad = (phase * 2* pi);         % convert phase to radians: 0 -> 2*pi
thetaRad = (theta / 360) * 2*pi;    % convert theta (orientation) to radians


Xt  = Xm * cos(thetaRad);           % compute proportion of Xm for given orientation
Yt  = Ym * sin(thetaRad);           % compute proportion of Ym for given orientation
XYt = Xt + Yt;                      % sum X and Y components
XYf = XYt * cycles * 2*pi;          % convert to radians and scale by frequency
 
im = sin( XYf + phaseRad);          % make 2D sinewave


%% Apply contrast
im = im*contrast;
im = im * 127.5 + background;


%% Log 
log                 = [];
log.theta           = theta;
log.f_spat          = f_spat;
log.phase           = phase;
log.contrast        = contrast;
log.imSize          = imSize;
log.background      = background;
