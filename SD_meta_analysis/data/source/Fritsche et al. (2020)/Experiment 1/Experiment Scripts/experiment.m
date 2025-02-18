clear all;

% Randomize seed of random number generator
% RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));
rng('shuffle');

%% Main Experiment Script

%% Define Gparams
Gparams = [];
Gparams.curSubj = input('Subject nr: ');

Gparams.BGint = 128;
Gparams.DistToScreen = .57;
Gparams.ScreenWidth = .48;

allScreens = Screen('Screens');
Gparams.pWindow = Screen('OpenWindow', allScreens(1), Gparams.BGint);
Screen('BlendFunction', Gparams.pWindow, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
[Gparams.ScreenResX, Gparams.ScreenResY] = Screen('WindowSize', Gparams.pWindow);
Gparams.flipinterval = Screen('GetFlipInterval',Gparams.pWindow);
Screen('TextSize',Gparams.pWindow, 16);
Screen('TextFont',Gparams.pWindow,'Helvatica');


Gparams.flip_lag = .015;

Gparams.dial.r_line = deg2pix(5, Gparams);
Gparams.dial.line_width = 15;
Gparams.dial.sd = 0.8;


Gparams.trial.stimdur = .5 - Gparams.flip_lag;
Gparams.trial.maskdur = 1 - Gparams.flip_lag;
Gparams.trial.responsedelay = 0.25 - Gparams.flip_lag;
Gparams.trial.prefixdur = 2 - Gparams.flip_lag;
% Gparams.trial.postfixdur = 1 - Gparams.flip_lag;

% Spatial frequencies for the gratings
Gparams.stim.f_spat1 = 0.33; % in cycles per degrees
Gparams.stim.f_spat2 = 0.5;

Gparams.stim.sigma = 1.5; % sd of gaussian envelope in degrees
Gparams.stim.mask_sd_smooth = .5;  % in degrees
Gparams.stim.eccentricity = 6.5; % in degrees
Gparams.stim.michelson = 0.25; % Michelson contrast of the stimulus
Gparams.stim.color = [255 255 255];

Gparams.background = [128 128 128]; % rgb color for grating background

Gparams.FixDotSize = deg2pix(.3, Gparams);

Gparams.datestr = datestr(now, 'HH:MM:SS dd-mmm-yyyy');
Gparams.DataPath = ['Data - Subj', num2str(Gparams.curSubj), '/'];
Gparams.LogPath = [Gparams.DataPath, 'Log/'];

if (~exist(Gparams.DataPath))
    mkdir(Gparams.DataPath);
end
if (~exist(Gparams.LogPath))
    mkdir(Gparams.LogPath);
end


%% Initialize
KbName('UnifyKeyNames');
HideCursor;
ListenChar(2);
priorityLevel = MaxPriority(Gparams.pWindow);
Priority(priorityLevel);

fullscreenloop(Gparams);

%% Instructions
instructions(Gparams);


%% Experiment
logPath = cell(1, 1);
dataPath = cell(1, 1);

Nblocks = 10;


% Side of stimulus presentation (left - right on interleaved runs)
% counterbalanced for subjects
% If subject's number is even start on the right side, else start on the
% left side
position = ones(1,Nblocks);

if mod(Gparams.curSubj,2) == 0
    position(2:2:Nblocks) = -1;
else 
    position(1:2:Nblocks) = -1;
end

position = position * deg2pix(Gparams.stim.eccentricity,Gparams);


for a = 1:Nblocks
    
    if position(a) == deg2pix(Gparams.stim.eccentricity,Gparams)
        side = 'right';
    else
        side = 'left';
    end
    showtext(sprintf(['Press a key to begin block %g/%g... \n\nIn this block the patterns will be presented on the ' side ' side. \n\nPlease remember to maintain fixation.'], a, Nblocks),Gparams);

    Screen('DrawDots', Gparams.pWindow, [(Gparams.ScreenResX - 1)/2, (Gparams.ScreenResY-1)/2], Gparams.FixDotSize, 255, [], 1);
    Screen('Flip', Gparams.pWindow);
    WaitSecs(2);
    
	cfg = [];
	cfg.curBlock = a;
    cfg.reps = 10; % The number of times each condition is preceded by every other nTuple. Determines the amount of trials.
    cfg.cbOrder = 2; % Depth of counterbalancing
    cfg.position = position(a);

	[dataPath{cfg.curBlock}, logPath{cfg.curBlock}] = block(cfg, Gparams);

    if (a ~= Nblocks)
        showtext(sprintf('You just finished block %g/%g.\n\n\nFeel free to take a break.\n\n\nIf you are ready,\npress any key to continue...', a, Nblocks), Gparams);
    end
end



%% Collect data from all blocks
totaldata = [];
totaldata.dataPath = dataPath;
totaldata.logPath = logPath;
totaldata.Gparams = Gparams;

totalPath = [Gparams.DataPath, 'totaldata_Subj', num2str(Gparams.curSubj)];
while(exist([totalPath, '.mat'], 'file'))
    totalPath = [totalPath, '_1'];
end

save([totalPath, '.mat'], 'totaldata');

%% Clean-up
showtext('This is the end of the experiment.\n\n\nThank you very much!', Gparams);

Screen('CloseAll');
ShowCursor;
ListenChar(0);
Priority(0);