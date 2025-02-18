clear all;

%% Test Script for testing spatial frequencies
Gparams = [];

Gparams.BGint = 128;
Gparams.DistToScreen = .57;
Gparams.ScreenWidth = .48;

allScreens = Screen('Screens');
Gparams.pWindow = Screen('OpenWindow', allScreens(1), Gparams.BGint);
Screen('BlendFunction', Gparams.pWindow, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
[Gparams.ScreenResX, Gparams.ScreenResY] = Screen('WindowSize', Gparams.pWindow);
Gparams.flipinterval = Screen('GetFlipInterval',Gparams.pWindow);
Screen('TextSize',Gparams.pWindow, 20);

Gparams.flip_lag = .015;

% Spatial frequencies of the gratings
Gparams.stim.f_spat1 = 0.5; % in cycles per degrees
Gparams.stim.f_spat2 = 0.6;

% Gparams.stim.sigma = 1.5; % sd of gaussian envelope in degrees
Gparams.stim.sigma = 1.5;
Gparams.stim_size = 10; % in degrees
Gparams.stim.mask_sd_smooth = .5;  % in degrees
Gparams.stim.eccentricity = 11; % in degrees
Gparams.stim.michelson = 0.25; % Michelson contrast of the stimulus

Gparams.background = [128 128 128]; % rgb color for grating background

Gparams.FixDotSize = deg2pix(.5, Gparams);

%% Initialize
KbName('UnifyKeyNames');
HideCursor;
ListenChar(2);
priorityLevel = MaxPriority(Gparams.pWindow);
Priority(priorityLevel);

fullscreenloop(Gparams);

% Generate stimulus 1
[im1, log.stim1] = grating(...
    45,...
    Gparams.stim.f_spat1,...
    0.25,...
    Gparams.stim.sigma,...
    [255 255 255],...
    [0 0 0],...
    Gparams.stim.michelson,...
    Gparams.background,...
    Gparams);



% Generate stimulus 2
[im2, log.stim2] = grating2(...
    -30,...
    Gparams.stim.f_spat2,...
    0.25,...
    Gparams.stim.sigma,...
    [255 255 255],...
    [0 0 0],...
    Gparams.stim.michelson,...
    Gparams.background,...
    Gparams);

% % Generate mask
% [mask_im] = mask(...
%     Gparams.stim.sigma,...
%     Gparams.stim.mask_sd_smooth,...
%     [255 255 255],...
%     [0 0 0],...
%     Gparams.background,...
%     Gparams...
% );


% Stimulus position 1 
yPos2 = (Gparams.ScreenResY-1)/2;
xPos2 = (Gparams.ScreenResX - 1)/2 - deg2pix(Gparams.stim.eccentricity,Gparams);
[s1, s2, s3] = size(im1);
baseRect = [0 0 s1 s2];
dstRect1 =  CenterRectOnPointd(baseRect, xPos2, yPos2);

% Stimulus position 2
yPos2 = (Gparams.ScreenResY-1)/2;
xPos2 = (Gparams.ScreenResX - 1)/2 + deg2pix(Gparams.stim.eccentricity,Gparams);
[s1, s2, s3] = size(im2);
baseRect = [0 0 s1 s2];
dstRect2 =  CenterRectOnPointd(baseRect, xPos2, yPos2);

% % Stimulus position 2 
% yPos2 = (Gparams.ScreenResY-1)/2;
% xPos2 = (Gparams.ScreenResX - 1)/2;
% [s1, s2, s3] = size(im2);
% baseRect = [0 0 s1 s2];
% dstRect2 =  CenterRectOnPointd(baseRect, xPos2, yPos2);

%% Present stimuli

% Display stimulus
pTexture1 = Screen('MakeTexture', Gparams.pWindow, im1);
pTexture2 = Screen('MakeTexture', Gparams.pWindow, im2);
Screen('DrawTexture', Gparams.pWindow, pTexture1,[],dstRect1);
Screen('DrawTexture', Gparams.pWindow, pTexture2,[],dstRect2);
Screen('DrawDots', Gparams.pWindow, [(Gparams.ScreenResX - 1)/2, (Gparams.ScreenResY-1)/2], Gparams.FixDotSize, 255, [], 1);


Screen('Flip', Gparams.pWindow);
KbReleaseWait;
KbWait;

Screen('CloseAll');
ShowCursor;
ListenChar(0);
Priority(0);
