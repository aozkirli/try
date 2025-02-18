clear all;

%% Script for generating stimulus figures

Gparams = []; 

Gparams.BGint = 128;
Gparams.DistToScreen = .57;
Gparams.ScreenWidth = .528;

allScreens = Screen('Screens');
Gparams.pWindow = Screen('OpenWindow', allScreens(1), Gparams.BGint);
Screen('BlendFunction', Gparams.pWindow, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
[Gparams.ScreenResX, Gparams.ScreenResY] = Screen('WindowSize', Gparams.pWindow);


Gparams.stim.horizontal_eccentricity = 9; % in degrees
Gparams.stim.vertical_eccentricity = 5; % in degrees
   

[gauss_mask] = transparency_mask(1.5,9,128,Gparams);
 
[im1, ~] = grating(...
    -5,...
    0.33,...
    0,...
    0.25,...
    9 ,...
    128,...
    Gparams);

[im2, ~] = grating(...
    -15 ,...
    0.33,...
    0,...
    0.25,...
    9,...
    128,...
    Gparams);   
   

% Stimulus position
yPos = (Gparams.ScreenResY-1)/2 + deg2pix(Gparams.stim.vertical_eccentricity,Gparams);
xPos = (Gparams.ScreenResX - 1)/2 - deg2pix(Gparams.stim.horizontal_eccentricity,Gparams);
[s1, s2, s3] = size(im1);
baseRect = [0 0 s1 s2];
dstRect1 =  CenterRectOnPointd(baseRect, xPos, yPos);

yPos = (Gparams.ScreenResY-1)/2 + deg2pix(Gparams.stim.vertical_eccentricity,Gparams);
xPos = (Gparams.ScreenResX - 1)/2 - deg2pix(Gparams.stim.horizontal_eccentricity,Gparams);
[s1, s2, s3] = size(gauss_mask);
baseRect = [0 0 s1 s2];
dstRect1_1 =  CenterRectOnPointd(baseRect, xPos, yPos);

yPos = (Gparams.ScreenResY-1)/2 + deg2pix(Gparams.stim.vertical_eccentricity,Gparams);
xPos = (Gparams.ScreenResX - 1)/2 + deg2pix(Gparams.stim.horizontal_eccentricity,Gparams);
[s1, s2, s3] = size(im2);
baseRect = [0 0 s1 s2];
dstRect2 =  CenterRectOnPointd(baseRect, xPos, yPos);

yPos = (Gparams.ScreenResY-1)/2  + deg2pix(Gparams.stim.vertical_eccentricity,Gparams);
xPos = (Gparams.ScreenResX - 1)/2 +  deg2pix(Gparams.stim.horizontal_eccentricity,Gparams);
[s1, s2, s3] = size(gauss_mask);
baseRect = [0 0 s1 s2];
dstRect2_2=  CenterRectOnPointd(baseRect, xPos, yPos);


% Display stimulus
maskTexture = Screen('MakeTexture', Gparams.pWindow, gauss_mask);

stimTexture = Screen('MakeTexture', Gparams.pWindow, im1);
Screen('DrawTexture', Gparams.pWindow, stimTexture, [], dstRect1);
Screen('DrawTexture', Gparams.pWindow, maskTexture, [], dstRect1_1);

stimTexture = Screen('MakeTexture', Gparams.pWindow, im2);
Screen('DrawTexture', Gparams.pWindow, stimTexture, [], dstRect2);
Screen('DrawTexture', Gparams.pWindow, maskTexture, [], dstRect2_2);


Screen('DrawDots', Gparams.pWindow, [(Gparams.ScreenResX - 1)/2, (Gparams.ScreenResY-1)/2], 10, 255, [], 1);

Screen('Flip', Gparams.pWindow);

imageArray = Screen('GetImage', Gparams.pWindow); 

KbWait();

Screen('CloseAll');

imwrite(imageArray,'afc_example3 .png');
   