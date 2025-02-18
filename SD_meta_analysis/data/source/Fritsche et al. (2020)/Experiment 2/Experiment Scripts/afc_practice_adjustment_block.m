function  [dataPath, logPath] = afc_practice_adjustment_block(cfg, Gparams)

log = [];
log.cfg = cfg;

% Set up data and log paths
dataPath = [Gparams.DataPath, 'data_S', num2str(Gparams.subject_number),'_AFC_Practice_Adjustment_Block_', num2str(cfg.curBlock)];
while(exist([dataPath, '.mat'], 'file'))
    dataPath = [dataPath, '_1'];
end
dataPath = [dataPath, '.mat'];

logPath = [Gparams.LogPath, 'log_S', num2str(Gparams.subject_number),'_AFC_Practice_Adjustment_Block_', num2str(cfg.curBlock)];
while(exist([logPath, '.mat'], 'file'))
    logPath = [logPath, '_1'];
end
logPath = [logPath, '.mat'];


% Set up the trialsequence for this block
repetitions = cfg.repetitions;

h_location = [ones(1,repetitions/2), -1*ones(1,repetitions/2)];
v_location = [ones(1,repetitions/4), -1*ones(1,repetitions/4)...
                ones(1,repetitions/4), -1*ones(1,repetitions/4)];

ifc_ori = (Gparams.inducer2ifc.orientation_range(2)-Gparams.inducer2ifc.orientation_range(1))...
    .*rand(1,repetitions) + Gparams.inducer2ifc.orientation_range(1);

% determine orientations of distractor inducer in the 2IFC task
ifc_distractor_ori = (Gparams.inducer2ifc.orientation_range(2)-Gparams.inducer2ifc.orientation_range(1))...
    .*rand(1,repetitions) + Gparams.inducer2ifc.orientation_range(1);

trialmatrix = [h_location;...
               v_location;...
               ifc_ori;...
               ifc_distractor_ori];

% randomize trial order
trialmatrix = trialmatrix(:,randperm(size(trialmatrix,2)));
                  
log.trialmatrix = trialmatrix;

% trialmatrix:
% 1st row: horizontal location
% 2nd row: vertical location
% 3nd row: orientation of adjustment stimulus
% 4rd row: orientation of distractor in the adjustment task

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
    funcfg.h_location                       = trialmatrix(1,iTrial);
    funcfg.v_location                       = trialmatrix(2,iTrial);
    funcfg.inducer_orientation              = trialmatrix(3,iTrial);
    funcfg.inducer_distractor_orientation   = trialmatrix(4,iTrial);  
    
    % save information about trial in log structure
    log.trialinfo{iTrial} = [ ...
        funcfg.h_location,...
        funcfg.v_location,...
        funcfg.inducer_orientation,...
        funcfg.inducer_distractor_orientation];
    
   
    log.trialoutput{iTrial} = practice_adjustment_trial(funcfg, Gparams);
    
    data{iTrial} = [ ...
        Gparams.subject_number, ...
        cfg.curBlock, ...
        iTrial, ...
        funcfg.h_location,... % location of the reference stimulus
        funcfg.v_location,... % location of the reference stimulus
        trialmatrix(3,iTrial),... % inducer orientation
        log.trialoutput{iTrial}.adjustment_response.alpha, ...
        log.trialoutput{iTrial}.adjustment_response.rt, ...
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

errors = mod(data_mat(:,6) - data_mat(:, 7) + 90, 180) - 90;
mean_error = sqrt(mean(errors.^2));

fprintf(['\n\n Mean Adjustment Error: ' num2str(round(mean_error*100)/100) ' degrees...\n']);

end