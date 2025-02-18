function log = trial(funcfg, Gparams)

log = [];
log.funcfg = funcfg;
log.onsets = [];

% Pre-stim fixation
Screen('DrawDots', Gparams.pWindow, [(Gparams.ScreenResX - 1)/2, (Gparams.ScreenResY-1)/2], Gparams.FixDotSize, 255, [], 1);
time = Screen('Flip', Gparams.pWindow);
log.onsets = [log.onsets, time];

if funcfg.trialfreq == 1
    
    f_spat = Gparams.stim.f_spat1;

elseif funcfg.trialfreq == 2
    
    f_spat = Gparams.stim.f_spat2;
    
end


% Generate stimulus
[im, log.stim] = grating(...
    funcfg.orientation,...
    f_spat,...
    0.25,...
    Gparams.stim.sigma,...
    Gparams.stim.color,...
    [0 0 0],...
    Gparams.stim.michelson,...
    Gparams.background,...
    Gparams...
);

% Generate mask
[mask_im] = mask(...
    Gparams.stim.sigma,...
    Gparams.stim.mask_sd_smooth,...
    [255 255 255],...
    [0 0 0],...
    Gparams.background,...
    Gparams...
);

% Stimulus position
yPos = (Gparams.ScreenResY-1)/2;
xPos = (Gparams.ScreenResX - 1)/2 + funcfg.position;
[s1, s2, s3] = size(im);
baseRect = [0 0 s1 s2];
dstRect =  CenterRectOnPointd(baseRect, xPos, yPos);

% Display stimulus
stimTexture = Screen('MakeTexture', Gparams.pWindow, im);
Screen('DrawTexture', Gparams.pWindow, stimTexture, [], dstRect);
Screen('DrawDots', Gparams.pWindow, [(Gparams.ScreenResX - 1)/2, (Gparams.ScreenResY-1)/2], Gparams.FixDotSize, 255, [], 1);
time = Screen('Flip', Gparams.pWindow, time + Gparams.trial.prefixdur - Gparams.flipinterval); % substract 1 frame here that is spent outside of the function
log.onsets = [log.onsets, time];

% Display mask
maskTexture = Screen('MakeTexture', Gparams.pWindow, mask_im);
Screen('DrawTexture', Gparams.pWindow, maskTexture, [], dstRect);
Screen('DrawDots', Gparams.pWindow, [(Gparams.ScreenResX - 1)/2, (Gparams.ScreenResY-1)/2], Gparams.FixDotSize, 255, [], 1);
time = Screen('Flip', Gparams.pWindow, time + Gparams.trial.stimdur);
log.onsets = [log.onsets, time];

Screen('DrawDots', Gparams.pWindow, [(Gparams.ScreenResX - 1)/2, (Gparams.ScreenResY-1)/2], Gparams.FixDotSize, 255, [], 1);
time = Screen('Flip', Gparams.pWindow, time + Gparams.trial.maskdur);
log.onsets = [log.onsets, time];

WaitSecs('UntilTime', time + Gparams.trial.responsedelay);



% Reponse dial
log.dial = responsedial(funcfg, Gparams);

% ITI
Screen('DrawDots', Gparams.pWindow, [(Gparams.ScreenResX - 1)/2, (Gparams.ScreenResY-1)/2], Gparams.FixDotSize, 255, [], 1);
time = Screen('Flip', Gparams.pWindow);
log.onsets = [log.onsets, time];

% WaitSecs('UntilTime', time + Gparams.trial.postfixdur);

% Close textures
Screen('Close', stimTexture);
Screen('Close', maskTexture);

end
