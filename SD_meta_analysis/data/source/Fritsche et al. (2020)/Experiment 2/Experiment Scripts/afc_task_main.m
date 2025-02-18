clear all; close all;

%% 2AFC Task

%% Set up general experiment parameters
Gparams         = [];
Gparams.task    = 'afc';
Gparams         = defGParams(Gparams);

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


%% create the trialmatrix
nBlocks         = 4;
repetitions     = 8;

stim_levels     = repmat(Gparams.sd2afc.stim_levels,1,repetitions*nBlocks);
% on half of the trials inducer and afc stimulus location is the same, on the other half the location is switched
loc_change      = [zeros(1,size(stim_levels,2)/2), ones(1,size(stim_levels,2)/2)];
% half of the trials is biased cw (1), other half ccw (-1)
inducer_bias    = [ones(1,size(stim_levels,2)/4), -1*ones(1,size(stim_levels,2)/4)...
                    ones(1,size(stim_levels,2)/4), -1*ones(1,size(stim_levels,2)/4)];
% half of the trials is presented on left side (-1) other half on right side (1)
h_location        = [ones(1,size(stim_levels,2)/8), -1*ones(1,size(stim_levels,2)/8)...
                    ones(1,size(stim_levels,2)/8), -1*ones(1,size(stim_levels,2)/8)...
                    ones(1,size(stim_levels,2)/8), -1*ones(1,size(stim_levels,2)/8)...
                    ones(1,size(stim_levels,2)/8), -1*ones(1,size(stim_levels,2)/8)];
% half of the inducers is presented in the upper visual field (-1) other half in the lower visual field (1)
v_location        = [ones(1,size(stim_levels,2)/16), -1*ones(1,size(stim_levels,2)/16)...
                    ones(1,size(stim_levels,2)/16), -1*ones(1,size(stim_levels,2)/16)...
                    ones(1,size(stim_levels,2)/16), -1*ones(1,size(stim_levels,2)/16)...
                    ones(1,size(stim_levels,2)/16), -1*ones(1,size(stim_levels,2)/16)...
                    ones(1,size(stim_levels,2)/16), -1*ones(1,size(stim_levels,2)/16)...
                    ones(1,size(stim_levels,2)/16), -1*ones(1,size(stim_levels,2)/16)...
                    ones(1,size(stim_levels,2)/16), -1*ones(1,size(stim_levels,2)/16)...
                    ones(1,size(stim_levels,2)/16), -1*ones(1,size(stim_levels,2)/16)];                
           
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

% trialmatrix:
% 1st row: stimulus level
% 2nd row: location change between inducer and afc stimulus (0 = no; 1 = yes)
% 3rd row: bias direction of inducer stimuli
% 4th row: horizontal screen location
% 5th row: vertical screen location (of inducer)
% 6th row: orientation of sd stimulus
% 7th row: orientation of distractor in the adjustment task

% Save the trialmatrix
savePath = [Gparams.DataPath, 'AFC_trialmatrix_S', num2str(Gparams.subject_number),'_Phase' num2str(Gparams.phase)];
while(exist([savePath, '.mat'], 'file'))
    savePath = [savePath, '_1'];
end
savePath = [savePath, '.mat'];

save(savePath,'trialmatrix')
% load(savePath); % Comment this in to load an old trialmatrix in case the script crashed !!! Comment out saving the new trialmatrix !!!



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

for iBlock = 1:nBlocks
    
    trialstart = size(Gparams.sd2afc.stim_levels,2)*repetitions * (iBlock-1) + 1;
    trialend   = size(Gparams.sd2afc.stim_levels,2)*repetitions * iBlock;
    
    cfg                 = [];
    cfg.curBlock        = iBlock;
    cfg.afc_response    = afc_response; % 1: which grating is more cw; 2: which grating is more ccw
    cfg.trialmatrix     = trialmatrix(:,trialstart:trialend);
    
    if cfg.afc_response == 1
        task_msg = [orientation_task_msg,afc_response_cw_msg];
    elseif cfg.afc_response == 2
        task_msg = [orientation_task_msg,afc_response_ccw_msg];
    end
    
    msg = ['This is block ' num2str(iBlock) ' of 4...\n\n\n\n\n' ...
        'Tasks in this block:\n\n\n ' task_msg...
        '\n\n***Please remember to maintain fixation on the central dot at all times***\n\n\n'...
        'Press any button to start with the block...'];
    
    showtext(msg,Gparams)
    
    afc_block(cfg, Gparams);
    
    if iBlock < 4       
        msg = 'You can take a short break now.\n\n\nPress any button when you are ready for the next block...';
    else
        msg = 'This is the end of the this part.\n\n\nPlease report to the experimenter now...'; 
    end
    
    showtext(msg,Gparams);
    
end

%% Cleaning up
Screen('CloseAll');
ShowCursor;
ListenChar(0);
Priority(0);
