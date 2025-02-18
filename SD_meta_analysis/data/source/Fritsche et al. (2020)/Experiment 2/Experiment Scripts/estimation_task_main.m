clear all; close all;

%% Estimation Task

%% Set up general experiment parameters
Gparams         = [];
Gparams.task    = 'estimation';
Gparams         = defGParams(Gparams);

nBlocks         = 8;

%% Initialize Screen
HideCursor;
priorityLevel = MaxPriority(Gparams.pWindow);
Priority(priorityLevel);

fullscreenloop(Gparams); % ensures that window is fullscreen

estimation_task_msg    = ['Task: *** Adjust orientation ***\n\n Adjust the response bar to the orientation of the grating.\n\n'...
    'Use the arrow keys to rotate the response bar. Press the Space bar to submit your response\n\n\n\n\nPlease press any key to continue...'];

showtext(estimation_task_msg,Gparams);

% Side of stimulus presentation (left - right on interleaved runs)
% counterbalanced for subjects
% If subject's number is even start on the right side, else start on the
% left side
h_position = ones(1,nBlocks);

if mod(Gparams.subject_number,2) == 0
    h_position(2:2:nBlocks) = -1;
else 
    h_position(1:2:nBlocks) = -1;
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
    cfg.reps        = 25;
    cfg.cbOrder     = 1;
      
    estimation_block(cfg, Gparams);
    
    if (iBlock ~= nBlocks)
        showtext(sprintf('You just finished block %g/%g.\n\n\nFeel free to take a break.\n\n\nIf you are ready,\npress any key to continue...', iBlock, nBlocks), Gparams);
    end
    
end

msg = 'This is the end of the this part.\n\n\nPlease report to the experimenter now...'; 
showtext(msg,Gparams);

%% Cleaning up
Screen('CloseAll');
ShowCursor;
ListenChar(0);
Priority(0);