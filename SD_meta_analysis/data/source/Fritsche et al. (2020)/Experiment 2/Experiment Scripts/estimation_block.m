function  [dataPath, logPath] = estimation_block(cfg, Gparams)

log     = [];
log.cfg = cfg;

% Set up data and log paths
dataPath = [Gparams.DataPath, 'data_S', num2str(Gparams.subject_number), '_Estimation_Block_', num2str(cfg.curBlock)];
while(exist([dataPath, '.mat'], 'file'))
    dataPath = [dataPath, '_1'];
end
dataPath = [dataPath, '.mat'];

logPath = [Gparams.LogPath, 'log_S', num2str(Gparams.subject_number),'_Estimation_Block_', num2str(cfg.curBlock)];
while(exist([logPath, '.mat'], 'file'))
    logPath = [logPath, '_1'];
end
logPath = [logPath, '.mat'];


% generate trialsequence
cb_sequence  = carryoverCounterbalance(2,cfg.cbOrder,cfg.reps,0);
tmp = [-1 1];
v_pos_sequence  = tmp(cb_sequence);

log.orientation = nan(size(v_pos_sequence,2), 1);
log.trial_hpos  = nan(size(v_pos_sequence,2), 1);
log.trial_vpos  = nan(size(v_pos_sequence,2), 1);

data            = cell(size(v_pos_sequence,2), 1);
log.trial       = cell(size(v_pos_sequence,2), 1);

% Show fixation screen before first trial
Screen('DrawDots', Gparams.pWindow, [(Gparams.ScreenResX - 1)/2, (Gparams.ScreenResY-1)/2], Gparams.FixDotSize, 255, [], 1);
time = Screen('Flip', Gparams.pWindow);
WaitSecs('UntilTime', time + Gparams.timing.estimation_iti + 2);

% loop over trials
for iTrial = 1:size(v_pos_sequence,2)
    
    funcfg                  = [];
    funcfg.orientation      = rand*180 - 90;
    funcfg.trial_hpos       = cfg.h_position;
    funcfg.trial_vpos       = v_pos_sequence(iTrial);
    
    log.orientation(iTrial) = funcfg.orientation;
    log.trial_hpos(iTrial)  = funcfg.trial_hpos;
    log.trial_vpos(iTrial)  = funcfg.trial_vpos;
    
    
    log.trial{iTrial} = estimation_trial(funcfg, Gparams);
    
    data{iTrial} = [ ...
        Gparams.subject_number, ...
        cfg.curBlock, ...
        iTrial, ...
        log.orientation(iTrial), ...
        log.trial{iTrial}.response.alpha, ...
        log.trial{iTrial}.response.rt, ...
        log.trial_vpos(iTrial)
        ];
    
    save(dataPath, 'data');
    save(logPath, 'log');
    
    % Interrupt trial?
    [~, ~, keyCode] = KbCheck;
    if (ismember(KbName('F8'), find(keyCode)))
        DrawFormattedText(Gparams.pWindow, 'Experiment temporarily interrupted by researcher.\n\nPlease wait...', 'center', 'center', 255);
        Screen('Flip', Gparams.pWindow);
        
        while(KbCheck);
            WaitSecs(0.01);
        end;
        
        keyCode = 0;
        while(~ismember(KbName('F8'), find(keyCode)))
            [~, ~, keyCode] = KbCheck;
            WaitSecs(0.01);
        end
        
        Screen('Flip', Gparams.pWindow);
        
        WaitSecs(2);
    end
    
end

% Feedback
data_mat = cell2mat(data);

errors = mod(data_mat(:,4) - data_mat(:,5) + 90, 180) - 90;
mean_error = sqrt(mean(errors.^2));

msg = ['Performance on this block:\n\n\n\n Mean adjustment error: ' num2str(round(mean_error*100)/100) ' degrees\n\n\n\nPlease press any button...'];

showtext(msg,Gparams)

end