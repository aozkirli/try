clear all; close all;

%% Estimation Task - Practice round

%% Set up general experiment parameters
Gparams         = [];
Gparams         = defGParams(Gparams);

nBlocks         = 1;


%% Initialize Screen
HideCursor;
priorityLevel = MaxPriority(Gparams.pWindow);
Priority(priorityLevel);

fullscreenloop(Gparams); % ensures that window is fullscreen

estimation_task_msg    = ['Task: *** Adjust orientation ***\n\n Adjust the response bar to the orientation of the pattern.\n\n'...
    'Use the arrow keys to rotate the response bar. Press the Space bar to submit your response.\n\nYou have 3 seconds for your response. If you were too slow, the fixation dot will turn red and the next trial will begin.\n\n'...
    'Please respond accurately, but within the time limit!  \n\n\n\n\nPlease press any key to continue...'];

showtext(estimation_task_msg,Gparams);

% Side of stimulus presentation (left - right on interleaved runs)
% counterbalanced for subjects and session
h_position = ones(1,nBlocks);

if mod(Gparams.subject_number,2) == 0
    if mod(Gparams.session_number,10) == 1
        h_position(2:2:nBlocks) = -1;
    elseif mod(Gparams.session_number,10) == 2
        h_position(1:2:nBlocks) = -1;
    end        
else 
    if mod(Gparams.session_number,10) == 1
        h_position(1:2:nBlocks) = -1;  
    elseif mod(Gparams.session_number,10) == 2
        h_position(2:2:nBlocks) = -1;
    end        
end

for iBlock = 1:nBlocks
    
    if h_position(iBlock) == 1
        side = 'right';
    else
        side = 'left';
    end
    
    showtext(sprintf(['Press a key to begin block %g/%g... \n\nIn this block the patterns will be presented on the ' side ' side. \n\nPlease remember to maintain fixation.'], iBlock, nBlocks),Gparams);

    cfg             = [];
    cfg.curBlock    = iBlock;
    cfg.h_position  = h_position(iBlock);
    cfg.reps        = 9;
    cfg.cbOrder     = 1;
    cfg.practice    = true;
      
    estimation_block(cfg, Gparams);
    
    if (iBlock ~= nBlocks)
        showtext(sprintf('You just finished block %g/%g.\n\n\nFeel free to take a break.\n\n\nIf you are ready,\npress any key to continue...', iBlock, nBlocks), Gparams);
    end
    
end

msg = 'This is the end of the practice block.\n\n\nPlease report to the experimenter now...'; 
showtext(msg,Gparams);

%% Cleaning up
Screen('CloseAll');
ShowCursor;
ListenChar(0);
Priority(0);


%% Console Stats
for iBlock = 1:nBlocks
    
    tmp = load(['../Data/S' num2str(Gparams.subject_number) '/practice_data_S' num2str(Gparams.subject_number) '_Session_' num2str(Gparams.session_number) '_Block_' num2str(iBlock) '.mat']);
    data_mat = cell2mat(tmp.data);
    valid_trials = not(isnan(data_mat(:,6)));
    errors = mod(data_mat(valid_trials,5) - data_mat(valid_trials,4) + 90, 180) - 90;
    mean_error = sqrt(mean(errors.^2));
    missed_trials = sum(isnan(data_mat(:,6)));
    
    disp(['Mean Error: ' num2str(round(mean_error*100)/100) ' deg. Missed trials: ' num2str(missed_trials)]);
    
end
