function  [dataPath, logPath] = afc_block(cfg, Gparams)

log = [];
log.cfg = cfg;

% Set up data and log paths
dataPath = [Gparams.DataPath, 'data_S', num2str(Gparams.subject_number),'_AFC_Phase_' num2str(Gparams.phase) '_Block_', num2str(cfg.curBlock)];
while(exist([dataPath, '.mat'], 'file'))
    dataPath = [dataPath, '_1'];
end
dataPath = [dataPath, '.mat'];

logPath = [Gparams.LogPath, 'log_S', num2str(Gparams.subject_number),'_AFC_Phase_' num2str(Gparams.phase) '_Block_', num2str(cfg.curBlock)];
while(exist([logPath, '.mat'], 'file'))
    logPath = [logPath, '_1'];
end
logPath = [logPath, '.mat'];


trialmatrix = cfg.trialmatrix;

log.trialmatrix = trialmatrix;

% trialmatrix:
% 1st row: stimulus level
% 2nd row: location change between inducer and afc stimulus (0 = no; 1 = yes)
% 3rd row: bias direction of inducer stimuli
% 4th row: horizontal screen location
% 5th row: vertical screen location (of inducer)
% 6th row: orientation of sd stimulus
% 7th row: orientation of distractor in the adjustment task

data            = cell(size(trialmatrix,2), 1);
log.trialinfo   = cell(size(trialmatrix,2), 1);
log.trialoutput = cell(size(trialmatrix,2), 1);

save(logPath, 'log');

% Show fixation screen before first trial
Screen('DrawDots', Gparams.pWindow, [(Gparams.ScreenResX - 1)/2, (Gparams.ScreenResY-1)/2], Gparams.FixDotSize, 255, [], 1);
time = Screen('Flip', Gparams.pWindow);
WaitSecs('UntilTime', time + Gparams.timing.afc_iti + 2);


% loop over trials
for iTrial = 1:size(trialmatrix,2)

    funcfg                                  = [];
    funcfg.afc_response                     = cfg.afc_response;
    funcfg.loc_change                       = trialmatrix(2,iTrial);
    funcfg.h_location                       = trialmatrix(4,iTrial);
    funcfg.v_location                       = trialmatrix(5,iTrial);
    funcfg.inducer_orientation              = trialmatrix(6,iTrial) + trialmatrix(3,iTrial) * 20;
    funcfg.inducer_distractor_orientation   = trialmatrix(7,iTrial);
    funcfg.sd_ref_stim_orientation          = trialmatrix(6,iTrial);
    funcfg.sd_probe_stim_orientation        = trialmatrix(6,iTrial) + trialmatrix(1,iTrial);

    % save information about trial in log structure
    log.trialinfo{iTrial} = [ ...
        funcfg.afc_response,...
        funcfg.loc_change,...
        funcfg.h_location,...
        funcfg.v_location,...
        funcfg.inducer_orientation,...
        funcfg.inducer_distractor_orientation,...
        funcfg.sd_ref_stim_orientation,...
        funcfg.sd_probe_stim_orientation];


    log.trialoutput{iTrial} = afc_trial(funcfg, Gparams);

    data{iTrial} = [ ...
        Gparams.subject_number, ...
        cfg.curBlock, ...
        cfg.afc_response,...
        iTrial, ...
        funcfg.loc_change,... % location change between inducer and test?
        funcfg.h_location,... % location of the reference stimulus
        trialmatrix(1,iTrial),... % stimulus level (direction of probe)
        trialmatrix(3,iTrial),... % inducer bias direction
        log.trialoutput{iTrial}.adjustment_response.alpha, ...
        log.trialoutput{iTrial}.adjustment_response.rt, ...
        log.trialoutput{iTrial}.afc_response.probe_dir,... % sd probe stimulus cw?
        log.trialoutput{iTrial}.afc_response.response,... % left/right button
        ];


    % save data
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

%% Feedback
data_mat = cell2mat(data);
trialinfo_mat = cell2mat(log.trialinfo);

errors = mod(trialinfo_mat(:,5) - data_mat(:, 9) + 90, 180) - 90;
mean_error = sqrt(mean(errors.^2));

pCorrect_2AFC = nansum((data_mat(:,7)>0 & data_mat(:,11) == 1) |  (data_mat(:,7)<0 & data_mat(:,11) == 0)) / nansum(data_mat(:,7) ~= 0);

if pCorrect_2AFC < 0.70
    afc_msg = 'Please try to improve your performance in this task!\n\n\n';
elseif pCorrect_2AFC > 0.90
    afc_msg = 'Wow, impressive! Keep up your good performance in this task!\n\n\n';
else
    afc_msg = '\n\n\n';
end


msg = ['Performance on this block:\n\n\n\n Task 1 - Mean adjustment error: ' num2str(round(mean_error*100)/100) '% degrees\n\n' ...
    'Task 2: ' num2str(round(pCorrect_2AFC*10000)/100) '% correct\n\n' afc_msg ...
    '\n\n Please press any button...'];

showtext(msg,Gparams)

end
