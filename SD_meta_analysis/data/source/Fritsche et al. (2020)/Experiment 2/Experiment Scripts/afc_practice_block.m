function  [dataPath, logPath] = afc_practice_block(cfg, Gparams)

log = [];
log.cfg = cfg;

% Set up data and log paths
dataPath = [Gparams.DataPath, 'data_S', num2str(Gparams.subject_number),'_AFC_Practice_Block_', num2str(cfg.curBlock)];
while(exist([dataPath, '.mat'], 'file'))
    dataPath = [dataPath, '_1'];
end
dataPath = [dataPath, '.mat'];

logPath = [Gparams.LogPath, 'log_S', num2str(Gparams.subject_number),'_AFC_Practice_Block_', num2str(cfg.curBlock)];
while(exist([logPath, '.mat'], 'file'))
    logPath = [logPath, '_1'];
end
logPath = [logPath, '.mat'];

% trialmatrix:
% 1st row: stimulus level
% 2nd row: location change between inducer and afc stimulus (0 = no; 1 = yes)
% 3rd row: bias direction of inducer stimuli
% 4th row: horizontal screen location
% 5th row: vertical screen location (of inducer)
% 6th row: orientation of sd stimulus
% 7th row: orientation of distractor in the adjustment task

stim_levels     = repmat(cfg.stim_levels,1,cfg.repetitions);
% on half of the trials inducer and afc stimulus location is the same, on the other half the location is switched
loc_change      = [zeros(1,size(stim_levels,2)/2), ones(1,size(stim_levels,2)/2)];
% half of the trials is biased cw (1), other half ccw (-1)
inducer_bias    = [ones(1,size(stim_levels,2)/4), -1*ones(1,size(stim_levels,2)/4)...
    ones(1,size(stim_levels,2)/4), -1*ones(1,size(stim_levels,2)/4)];
% half of the trials is presented on left side (-1) other half on right side (1)
h_location      = [ones(1,size(stim_levels,2)/8), -1*ones(1,size(stim_levels,2)/8)...
    ones(1,size(stim_levels,2)/8), -1*ones(1,size(stim_levels,2)/8)...
    ones(1,size(stim_levels,2)/8), -1*ones(1,size(stim_levels,2)/8)...
    ones(1,size(stim_levels,2)/8), -1*ones(1,size(stim_levels,2)/8)];

v_location = repmat([-1 1],1,cfg.repetitions);
v_location = v_location(randperm(size(v_location,2)));

% randomly draw orientations for 2AFC reference stimulus -
% for each trial draw from specific interval such that the reference and
% probe stimulus are guaranteed to fall in the pre-specified interal (see defGparams.m)
for i = 1:size(stim_levels,2)
    if stim_levels(i) <= 0
        r1 = Gparams.sd2afc.orientation_range(1) - stim_levels(i);
        r2 = Gparams.sd2afc.orientation_range(2);
        sd2afc_orientations(i) = (r2-r1).*rand + r1;
    elseif stim_levels(i) > 0
        r1 = Gparams.sd2afc.orientation_range(1);
        r2 = Gparams.sd2afc.orientation_range(2) - stim_levels(i);
        sd2afc_orientations(i) = (r2-r1).*rand + r1;
    end
end


% determine orientations of distractors in the 2IFC task
ifc_distractor_ori = (Gparams.inducer2ifc.orientation_range(2)-Gparams.inducer2ifc.orientation_range(1))...
    .*rand(1,size(stim_levels,2)) + Gparams.inducer2ifc.orientation_range(1);


trialmatrix = [stim_levels;...
    loc_change;...
    inducer_bias;...
    h_location;...
    v_location;...
    sd2afc_orientations;...
    ifc_distractor_ori];


% randomize trial order
trialmatrix = trialmatrix(:,randperm(size(trialmatrix,2)));


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

fprintf(['\n\n Mean Adjustment Error: ' num2str(round(mean_error*100)/100) ' degrees...\n']);
fprintf(['\n\n 2AFC Accuracy: ' num2str(round(pCorrect_2AFC*10000)/100) ' percent correct...\n']);

end