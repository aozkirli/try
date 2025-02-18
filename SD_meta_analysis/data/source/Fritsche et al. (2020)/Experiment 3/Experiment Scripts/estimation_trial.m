function log = estimation_trial(funcfg, Gparams)

log = [];
log.funcfg = funcfg;
log.onsets = [];

% present cue
Screen('DrawDots', Gparams.pWindow, [(Gparams.ScreenResX - 1)/2, (Gparams.ScreenResY-1)/2], Gparams.FixDotSize, 180, [], 1);
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


% Presentation position
h_pos = funcfg.trial_hpos * deg2pix(Gparams.stim.h_eccentricity,Gparams);

yPos = (Gparams.ScreenResY-1)/2; 
xPos = (Gparams.ScreenResX - 1)/2 + h_pos;
[s1, s2, ~] = size(stim);
baseRect = [0 0 s1 s2];
dstRect =  CenterRectOnPointd(baseRect, xPos, yPos);


% display stimulus
stimTexture     = Screen('MakeTexture', Gparams.pWindow, stim);
gaussTexture    = Screen('MakeTexture', Gparams.pWindow, gaussian_window);

Screen('DrawTexture', Gparams.pWindow, stimTexture, [], dstRect);
Screen('DrawTexture', Gparams.pWindow, gaussTexture, [], dstRect);

Screen('DrawDots', Gparams.pWindow, [(Gparams.ScreenResX - 1)/2, (Gparams.ScreenResY-1)/2], Gparams.FixDotSize, 255, [], 1);
time = Screen('Flip', Gparams.pWindow, time + Gparams.timing.prestim_fixation);
log.onsets = [log.onsets, time];


% Display mask
maskTexture = Screen('MakeTexture', Gparams.pWindow, noise_mask);

Screen('DrawTexture', Gparams.pWindow, maskTexture, [], dstRect);
Screen('DrawTexture', Gparams.pWindow, gaussTexture, [], dstRect);

Screen('DrawDots', Gparams.pWindow, [(Gparams.ScreenResX - 1)/2, (Gparams.ScreenResY-1)/2], Gparams.FixDotSize, 255, [], 1);
time = Screen('Flip', Gparams.pWindow, time + Gparams.timing.stimdur);
log.onsets = [log.onsets, time];

% response delay
Screen('DrawDots', Gparams.pWindow, [(Gparams.ScreenResX - 1)/2, (Gparams.ScreenResY-1)/2], Gparams.FixDotSize, 255, [], 1);
time_start_response_interval = Screen('Flip', Gparams.pWindow, time + Gparams.timing.maskdur);
log.onsets = [log.onsets, time_start_response_interval];

Screen('DrawDots', Gparams.pWindow, [(Gparams.ScreenResX - 1)/2, (Gparams.ScreenResY-1)/2], Gparams.FixDotSize, 255, [], 1);
time = Screen('Flip', Gparams.pWindow, time_start_response_interval + Gparams.timing.response_delay(funcfg.response_delay));
log.onsets = [log.onsets, time];

% response
cfg = [];
cfg.h_pos = h_pos;
cfg.maxtime = Gparams.timing.response_timeout;
cfg.fadeout_time = 0;
log.response = responsedial(cfg, Gparams);

% warning when response timed out
if isnan(log.response.rt)
    
    % feedback - red warning
    Screen('DrawDots', Gparams.pWindow, [(Gparams.ScreenResX - 1)/2, (Gparams.ScreenResY-1)/2], Gparams.FixDotSize, [255 0 0], [], 1);
    time = Screen('Flip', Gparams.pWindow);
    log.onsets = [log.onsets, time];
      
    % ITI fixation
    Screen('DrawDots', Gparams.pWindow, [(Gparams.ScreenResX - 1)/2, (Gparams.ScreenResY-1)/2], Gparams.FixDotSize, 255, [], 1);
    Screen('Flip', Gparams.pWindow, time + Gparams.timing.warning);
    
else
    
    % feedback - response registered
    Screen('DrawDots', Gparams.pWindow, [(Gparams.ScreenResX - 1)/2, (Gparams.ScreenResY-1)/2], Gparams.FixDotSize, 255, [], 1);
    time = Screen('Flip', Gparams.pWindow);
    log.onsets = [log.onsets, time];
    
    % ITI fixation
    Screen('DrawDots', Gparams.pWindow, [(Gparams.ScreenResX - 1)/2, (Gparams.ScreenResY-1)/2], Gparams.FixDotSize, 255, [], 1);
    Screen('Flip', Gparams.pWindow, time + Gparams.timing.warning);
    
end


% Wait until the complete response period + half of ITI is over
Screen('DrawDots', Gparams.pWindow, [(Gparams.ScreenResX - 1)/2, (Gparams.ScreenResY-1)/2], Gparams.FixDotSize, 255, [], 1);
time = Screen('Flip', Gparams.pWindow, time_start_response_interval + Gparams.timing.response_interval);
log.onsets = [log.onsets, time];

% close textures
Screen('Close', stimTexture);
Screen('Close', gaussTexture);
Screen('Close', maskTexture);

end
