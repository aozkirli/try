clear all;

%% Script for generating stimulus figures

Gparams = [];

Gparams.BGint = 128;
Gparams.DistToScreen = .57;
Gparams.ScreenWidth = .48;

allScreens = Screen('Screens');
Gparams.pWindow = Screen('OpenWindow', allScreens(1), Gparams.BGint);
Screen('BlendFunction', Gparams.pWindow, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
[Gparams.ScreenResX, Gparams.ScreenResY] = Screen('WindowSize', Gparams.pWindow);

Screen('CloseAll');

Gparams.stim.f_spat = 0.4;

Gparams.stim.sigma = 1.5; % sd of gaussian envelope in degrees
Gparams.stim_size = 10; % in degrees
Gparams.stim.mask_sd_smooth = .5;  % in degrees
Gparams.stim.eccentricity = 6.5; % in degrees
Gparams.stim.michelson = 0.25; % Michelson contrast of the stimulus

Gparams.background = [128 128 128]; % rgb color for grating background

% Generate stimulus 1
[im, log.stim] = grating(...
    Gparams.stim_size,...
    -35,...
    Gparams.stim.f_spat,...
    0,...
    Gparams.stim.sigma,...
    [255 255 255],...
    [0 0 0],...
    Gparams.stim.michelson,...
    Gparams.background,...
    Gparams);

[mask_im] = mask(...
    Gparams.stim_size,...
    Gparams.stim.sigma,...
    Gparams.stim.mask_sd_smooth,...
    [255 255 255],...
    [0 0 0],...
    Gparams.background,...
    Gparams...
);

im = im./255;

figure
image(im)
axis off
axis image

mask_im = mask_im./255;
figure
image(mask_im)
axis off
axis image


