clear all; close all;

%% Practice adjustment task

%% Set up general experiment parameters
prompt = {'Enter block number:'};
dlg_title = 'Block Information';
num_lines = 1;
answer = inputdlg(prompt,dlg_title,num_lines);

Gparams                 = [];
Gparams.task            = 'estimation';
Gparams.block_number    = str2num(answer{1});

Gparams = defGParams(Gparams);


%% Initialize Screen
HideCursor;
priorityLevel = MaxPriority(Gparams.pWindow);
Priority(priorityLevel);

fullscreenloop(Gparams); % ensures that window is fullscreen

orientation_task_msg    = ['*** Adjust orientation ***\n\n Adjust the response bar to the orientation of the grating.\n\n'...
    'Use the arrow keys to rotate the response bar. Press the Space bar to submit your response\n\n\n\n\n'];


cfg                 = [];
cfg.curBlock        = Gparams.block_number;
cfg.repetitions     = 20;


msg = [orientation_task_msg ...
    '\n\n***Please remember to maintain fixation on the central dot at all times***\n\n\n'...
    'Press any button to start with the block...'];

showtext(msg,Gparams)

afc_practice_adjustment_block(cfg, Gparams);

%% Cleaning up
Screen('CloseAll');
ShowCursor;
ListenChar(0);
Priority(0);