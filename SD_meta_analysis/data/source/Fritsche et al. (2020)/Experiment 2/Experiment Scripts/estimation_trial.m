function log = estimation_trial(funcfg, Gparams)

log = [];
log.funcfg = funcfg;
log.onsets = [];

% Presentation position
h_pos = funcfg.trial_hpos * deg2pix(Gparams.stim.h_eccentricity,Gparams);     
v_pos = funcfg.trial_vpos * deg2pix(Gparams.stim.v_eccentricity,Gparams); 

% cue mask
[cue_mask,~] = transparency_mask(...
    0.7,...
    Gparams.stim.imSize,...
    Gparams.bg,...
    Gparams);

% cue position
yPos = (Gparams.ScreenResY-1)/2 + v_pos; 
xPos = (Gparams.ScreenResX - 1)/2 + h_pos;
[s1, s2, ~] = size(cue_mask);
baseRect = [0 0 s1 s2];
dstRect =  CenterRectOnPointd(baseRect, xPos, yPos);

% present cue
Screen('DrawDots', Gparams.pWindow, [(Gparams.ScreenResX - 1)/2, (Gparams.ScreenResY-1)/2], Gparams.FixDotSize, 255, [], 1);
Screen('FillRect', Gparams.pWindow, 170, dstRect);
cuemaskTexture = Screen('MakeTexture', Gparams.pWindow, cue_mask);
Screen('DrawTexture', Gparams.pWindow, cuemaskTexture , [], dstRect);
time = Screen('Flip', Gparams.pWindow);
log.onsets = [log.onsets, time];

% pre-stimulus fixation
Screen('DrawDots', Gparams.pWindow, [(Gparams.ScreenResX - 1)/2, (Gparams.ScreenResY-1)/2], Gparams.FixDotSize, 255, [], 1);
time = Screen('Flip', Gparams.pWindow, time + Gparams.timing.cuedur);
log.onsets = [log.onsets, time];

% generate the stimulus
[stim, log.stim] = grating(...
    funcfg.orientation,...
    Gparams.stim.f_spat,...
    0,...
    Gparams.stim.contrast,...
    Gparams.stim.imSize,...
    Gparams.bg,...
    Gparams);

% generate noise mask
[noise_mask,~] = mask(...
    Gparams.mask.sd_smooth,...
    Gparams.mask.contrast,...
    Gparams.stim.imSize,...
    Gparams.bg,...
    Gparams);

% generate transparency mask
[gaussian_window,~] = transparency_mask(...
    Gparams.stim.window_sd,...
    Gparams.stim.imSize,...
    Gparams.bg,...
    Gparams);


% display stimulus
stimTexture     = Screen('MakeTexture', Gparams.pWindow, stim);
gaussTexture    = Screen('MakeTexture', Gparams.pWindow, gaussian_window);

Screen('DrawTexture', Gparams.pWindow, stimTexture, [], dstRect);
Screen('DrawTexture', Gparams.pWindow, gaussTexture, [], dstRect);

Screen('DrawDots', Gparams.pWindow, [(Gparams.ScreenResX - 1)/2, (Gparams.ScreenResY-1)/2], Gparams.FixDotSize, 255, [], 1);
time = Screen('Flip', Gparams.pWindow, time + Gparams.timing.inducer_prestim);
log.onsets = [log.onsets, time];


% Display mask
maskTexture = Screen('MakeTexture', Gparams.pWindow, noise_mask);

Screen('DrawTexture', Gparams.pWindow, maskTexture, [], dstRect);
Screen('DrawTexture', Gparams.pWindow, gaussTexture, [], dstRect);

Screen('DrawDots', Gparams.pWindow, [(Gparams.ScreenResX - 1)/2, (Gparams.ScreenResY-1)/2], Gparams.FixDotSize, 255, [], 1);
time = Screen('Flip', Gparams.pWindow, time + Gparams.timing.inducer_stimdur);
log.onsets = [log.onsets, time];

% response delay
Screen('DrawDots', Gparams.pWindow, [(Gparams.ScreenResX - 1)/2, (Gparams.ScreenResY-1)/2], Gparams.FixDotSize, 255, [], 1);
time = Screen('Flip', Gparams.pWindow, time + Gparams.timing.inducer_maskdur);
log.onsets = [log.onsets, time];

Screen('DrawDots', Gparams.pWindow, [(Gparams.ScreenResX - 1)/2, (Gparams.ScreenResY-1)/2], Gparams.FixDotSize, 255, [], 1);
time = Screen('Flip', Gparams.pWindow, time + Gparams.timing.pre_response_fixation);
log.onsets = [log.onsets, time];

% response
cfg = [];
cfg.h_pos = h_pos;
cfg.v_pos = v_pos;
log.response = responsedial(cfg, Gparams);

% ITI
Screen('DrawDots', Gparams.pWindow, [(Gparams.ScreenResX - 1)/2, (Gparams.ScreenResY-1)/2], Gparams.FixDotSize, 255, [], 1);
time = Screen('Flip', Gparams.pWindow);
log.onsets = [log.onsets, time];

Screen('DrawDots', Gparams.pWindow, [(Gparams.ScreenResX - 1)/2, (Gparams.ScreenResY-1)/2], Gparams.FixDotSize, 255, [], 1);
time = Screen('Flip', Gparams.pWindow, time + Gparams.timing.estimation_iti);
log.onsets = [log.onsets, time];

% close textures
Screen('Close', cuemaskTexture);
Screen('Close', stimTexture);
Screen('Close', gaussTexture);
Screen('Close', maskTexture);

end
