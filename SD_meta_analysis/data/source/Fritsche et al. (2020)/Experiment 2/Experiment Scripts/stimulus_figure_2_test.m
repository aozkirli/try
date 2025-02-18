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
[gauss_mask2] = transparency_mask(1.3,9,128,Gparams);
 
[im1, ~] = grating(...
    20,...
    0.33,...
    0,...
    0.25,...
    9 ,...
    128,...
    Gparams);

[im2, ~] = grating(...
    20,...
    0.33,...
    0,...
    0.25,...
    9,...
    128,...
    Gparams);   

[im3, ~] = mask(...
    0.4,...
    0.8 ,...
    9,...
    128,...
    Gparams);
   

imSize = round(deg2pix(9,Gparams));

%create bar texture
bar = ones(imSize) * 128;
bar(:,round(imSize/2 + 1 - deg2pix(0.6, Gparams)/2) : round(imSize/2 + 1 + deg2pix(0.6, Gparams)/2)) = 255;
bartex = Screen('MakeTexture', Gparams.pWindow, bar);



% Stimulus position
yPos = (Gparams.ScreenResY-1)/2 - deg2pix(Gparams.stim.vertical_eccentricity,Gparams);
xPos = (Gparams.ScreenResX - 1)/2 - deg2pix(Gparams.stim.horizontal_eccentricity,Gparams);
[s1, s2, s3] = size(im1);
baseRect = [0 0 s1 s2];
dstRect1 =  CenterRectOnPointd(baseRect, xPos, yPos);

yPos = (Gparams.ScreenResY-1)/2 - deg2pix(Gparams.stim.vertical_eccentricity,Gparams);
xPos = (Gparams.ScreenResX - 1)/2 - deg2pix(Gparams.stim.horizontal_eccentricity,Gparams);
[s1, s2, s3] = size(gauss_mask);
baseRect = [0 0 s1 s2];
dstRect1_1 =  CenterRectOnPointd(baseRect, xPos, yPos);

yPos = (Gparams.ScreenResY-1)/2 + deg2pix(Gparams.stim.vertical_eccentricity,Gparams);
xPos = (Gparams.ScreenResX - 1)/2 - deg2pix(Gparams.stim.horizontal_eccentricity,Gparams);
[s1, s2, s3] = size(im2);
baseRect = [0 0 s1 s2];
dstRect2 =  CenterRectOnPointd(baseRect, xPos, yPos);

yPos = (Gparams.ScreenResY-1)/2 + deg2pix(Gparams.stim.vertical_eccentricity,Gparams);
xPos = (Gparams.ScreenResX - 1)/2 - deg2pix(Gparams.stim.horizontal_eccentricity,Gparams);
[s1, s2, s3] = size(gauss_mask);
baseRect = [0 0 s1 s2];
dstRect2_2=  CenterRectOnPointd(baseRect, xPos, yPos);

yPos = (Gparams.ScreenResY-1)/2 + deg2pix(Gparams.stim.vertical_eccentricity,Gparams);
xPos = (Gparams.ScreenResX - 1)/2 + deg2pix(Gparams.stim.horizontal_eccentricity,Gparams);
[s1, s2, s3] = size(im3);
baseRect = [0 0 s1 s2];
dstRect3 =  CenterRectOnPointd(baseRect, xPos, yPos);

yPos = (Gparams.ScreenResY-1)/2 + deg2pix(Gparams.stim.vertical_eccentricity,Gparams);
xPos = (Gparams.ScreenResX - 1)/2 + deg2pix(Gparams.stim.horizontal_eccentricity,Gparams);
[s1, s2, s3] = size(gauss_mask);
baseRect = [0 0 s1 s2];
dstRect3_3 =  CenterRectOnPointd(baseRect, xPos, yPos);

yPos = (Gparams.ScreenResY-1)/2 - deg2pix(Gparams.stim.vertical_eccentricity,Gparams);
xPos = (Gparams.ScreenResX - 1)/2 + deg2pix(Gparams.stim.horizontal_eccentricity,Gparams);
[s1, s2, s3] = size(bar);
baseRect = [0 0 s1 s2];
dstRect4 =  CenterRectOnPointd(baseRect, xPos, yPos);

yPos = (Gparams.ScreenResY-1)/2 - deg2pix(Gparams.stim.vertical_eccentricity,Gparams);
xPos = (Gparams.ScreenResX - 1)/2 + deg2pix(Gparams.stim.horizontal_eccentricity,Gparams);
[s1, s2, s3] = size(gauss_mask);
baseRect = [0 0 s1 s2];
dstRect4_4 =  CenterRectOnPointd(baseRect, xPos, yPos);

% Display stimulus
maskTexture = Screen('MakeTexture', Gparams.pWindow, gauss_mask);
maskTexture2 = Screen('MakeTexture', Gparams.pWindow, gauss_mask2);

stimTexture = Screen('MakeTexture', Gparams.pWindow, im1);
Screen('DrawTexture', Gparams.pWindow, stimTexture, [], dstRect1);
Screen('DrawTexture', Gparams.pWindow, maskTexture, [], dstRect1_1);

stimTexture = Screen('MakeTexture', Gparams.pWindow, im2);
Screen('DrawTexture', Gparams.pWindow, stimTexture, [], dstRect2);
Screen('DrawTexture', Gparams.pWindow, maskTexture, [], dstRect2_2);

stimTexture = Screen('MakeTexture', Gparams.pWindow, im3);
Screen('DrawTexture', Gparams.pWindow, stimTexture, [], dstRect3);
Screen('DrawTexture', Gparams.pWindow, maskTexture, [], dstRect3_3);
 
Screen('DrawTexture', Gparams.pWindow, bartex, [], dstRect4,20);
Screen('DrawTexture', Gparams.pWindow, maskTexture2, [], dstRect4_4,20);

Screen('DrawDots', Gparams.pWindow, [(Gparams.ScreenResX - 1)/2, (Gparams.ScreenResY-1)/2], 10, 255, [], 1);

Screen('Flip', Gparams.pWindow);

imageArray = Screen('GetImage', Gparams.pWindow); 

KbWait();

Screen('CloseAll');

%imwrite(imageArray,'layout.png');
   