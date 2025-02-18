clear all; close all;

%% Practice AFC Task

%% Set up general experiment parameters
prompt = {'Enter block number:'};
dlg_title = 'Block Information';
num_lines = 1;
answer = inputdlg(prompt,dlg_title,num_lines);

Gparams                 = [];
Gparams.task            = 'afc';
Gparams.block_number    = str2num(answer{1});

Gparams = defGParams(Gparams);

%% Session parameters
if strcmp(Gparams.response_order,'cw')
    afc_response = 1;
elseif strcmp(Gparams.response_order,'ccw')
    afc_response = 2;
else
    
    Screen('CloseAll');
    ShowCursor;
    ListenChar(0);
    Priority(0);
    
    disp('Invalid response specification!')
    
end


%% Initialize Screen
HideCursor;
priorityLevel = MaxPriority(Gparams.pWindow);
Priority(priorityLevel);

fullscreenloop(Gparams); % ensures that window is fullscreen

orientation_task_msg    = ['First task: *** Adjust orientation ***\n\n Adjust the response bar to the orientation of the grating.\n\n'...
    'Use the arrow keys to rotate the response bar. Press the Space bar to submit your response\n\n\n\n\n'];
afc_response_cw_msg     = ['Second task:\n\nWhich grating (left/right) is oriented more ***clockwise***? \n\n' ... 
    'Left arrow key for left grating, right arrow key for right grating\n\n\n\n']; 
afc_response_ccw_msg    = ['Second task:\n\nWhich grating (left/right) is oriented more ***counter-clockwise***? \n\n' ... 
    'Left arrow key for left grating, right arrow key for right grating\n\n\n\n'];


cfg                 = [];
cfg.curBlock        = Gparams.block_number;
cfg.stim_levels     = [-9 9];
cfg.repetitions     = 8;
cfg.afc_response    = afc_response; % 1: which grating is more cw; 2: which grating is more ccw

if cfg.afc_response == 1
    task_msg = [orientation_task_msg,afc_response_cw_msg];
elseif cfg.afc_response == 2
    task_msg = [orientation_task_msg,afc_response_ccw_msg];
end

msg = ['Tasks in this block:\n\n\n ' task_msg...
    '\n\n***Please remember to maintain fixation on the central dot at all times***\n\n\n'...
    'Press any button to start with the block...'];

showtext(msg,Gparams)

afc_practice_block(cfg, Gparams);

%% Cleaning up
Screen('CloseAll');
ShowCursor;
ListenChar(0);
Priority(0);
