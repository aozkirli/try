function log = practice_adjustment_trial(cfg, Gparams)

log = [];
log.cfg = cfg;
log.onsets = [];

%% Adjustment part
% Presentation position
h_pos = cfg.h_location * deg2pix(Gparams.stim.h_eccentricity,Gparams);     
v_pos = cfg.v_location * deg2pix(Gparams.stim.v_eccentricity,Gparams); 

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

% make stimuli
[inducer_stim, ~] = grating(...
    cfg.inducer_orientation,...,...
    Gparams.stim.f_spat,...
    0,...
    Gparams.stim.contrast,...
    Gparams.stim.imSize,...
    Gparams.bg,...
    Gparams);

[distractor_stim, ~] = grating(...
    cfg.inducer_distractor_orientation,...
    Gparams.stim.f_spat,...
    0,...
    Gparams.stim.contrast,...
    Gparams.stim.imSize,...
    Gparams.bg,...
    Gparams);

% generate transparency mask
[gaussian_window,~] = transparency_mask(...
    Gparams.stim.window_sd,...
    Gparams.stim.imSize,...
    Gparams.bg,...
    Gparams);


% determine stimulus positions
[s1, s2, ~] = size(inducer_stim);
baseRect = [0 0 s1 s2];

yPos1 = (Gparams.ScreenResY-1)/2 + v_pos;
xPos1 = (Gparams.ScreenResX - 1)/2 + h_pos;
inducer_dstRect =  CenterRectOnPointd(baseRect, xPos1, yPos1);

yPos2 = (Gparams.ScreenResY-1)/2 + v_pos;
xPos2 = (Gparams.ScreenResX - 1)/2 - h_pos;
distractor_dstRect =  CenterRectOnPointd(baseRect, xPos2, yPos2);

% display stimuli
gaussTexture = Screen('MakeTexture', Gparams.pWindow, gaussian_window);

inducerTexture = Screen('MakeTexture', Gparams.pWindow, inducer_stim);
Screen('DrawTexture', Gparams.pWindow, inducerTexture, [], inducer_dstRect);
Screen('DrawTexture', Gparams.pWindow, gaussTexture, [], inducer_dstRect);

distractorTexture = Screen('MakeTexture', Gparams.pWindow, distractor_stim);
Screen('DrawTexture', Gparams.pWindow, distractorTexture, [], distractor_dstRect);
Screen('DrawTexture', Gparams.pWindow, gaussTexture, [], distractor_dstRect);

Screen('DrawDots', Gparams.pWindow, [(Gparams.ScreenResX - 1)/2, (Gparams.ScreenResY-1)/2], Gparams.FixDotSize, 255, [], 1);
time = Screen('Flip', Gparams.pWindow, time + Gparams.timing.inducer_prestim);
log.onsets = [log.onsets, time];


% generate noise masks
[inducer_noise_mask,~] = mask(...
    Gparams.mask.sd_smooth,...
    Gparams.mask.contrast,...
    Gparams.stim.imSize,...
    Gparams.bg,...
    Gparams);

[distractor_noise_mask,~] = mask(...
    Gparams.mask.sd_smooth,...
    Gparams.mask.contrast,...
    Gparams.stim.imSize,...
    Gparams.bg,...
    Gparams);

% display mask
inducer_maskTexture = Screen('MakeTexture', Gparams.pWindow, inducer_noise_mask);
Screen('DrawTexture', Gparams.pWindow, inducer_maskTexture, [], inducer_dstRect);
Screen('DrawTexture', Gparams.pWindow, gaussTexture, [], inducer_dstRect);

distractor_maskTexture = Screen('MakeTexture', Gparams.pWindow, distractor_noise_mask);
Screen('DrawTexture', Gparams.pWindow, distractor_maskTexture, [], distractor_dstRect);
Screen('DrawTexture', Gparams.pWindow, gaussTexture, [], distractor_dstRect);

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
funcfg = [];
funcfg.h_pos = h_pos;
funcfg.v_pos = v_pos;
log.adjustment_response = responsedial(funcfg, Gparams);

% clean up
Screen('Close', cuemaskTexture);
Screen('Close', gaussTexture);
Screen('Close', inducerTexture);
Screen('Close', distractorTexture);
Screen('Close', inducer_maskTexture);
Screen('Close', distractor_maskTexture);

% ITI period
Screen('DrawDots', Gparams.pWindow, [(Gparams.ScreenResX - 1)/2, (Gparams.ScreenResY-1)/2], Gparams.FixDotSize, 255, [], 1);
time = Screen('Flip', Gparams.pWindow);
log.onsets = [log.onsets, time];

Screen('DrawDots', Gparams.pWindow, [(Gparams.ScreenResX - 1)/2, (Gparams.ScreenResY-1)/2], Gparams.FixDotSize, 255, [], 1);
time = Screen('Flip', Gparams.pWindow, time + Gparams.timing.afc_iti);
log.onsets = [log.onsets, time];