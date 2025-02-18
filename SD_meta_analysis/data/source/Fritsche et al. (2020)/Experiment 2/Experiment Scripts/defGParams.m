function Gparams = defGParams(Gparams)

%% Get information about experiment
switch (Gparams.task)
    
    case 'estimation'
        
        prompt = {'Enter subject number:','Enter location:'};
        dlg_title = 'Experiment Information';
        num_lines = 1;
        defaultans = {'','cubicle'};
        answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
        
        Gparams.subject_number      = str2num(answer{1});
        Gparams.location            = answer{2};
        
    case 'afc'
        
        prompt = {'Enter subject number:','Enter location:','Phase','Response'};
        dlg_title = 'Experiment Information';
        num_lines = 1;
        defaultans = {'','cubicle','',''};
        answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
        
        Gparams.subject_number      = str2num(answer{1});
        Gparams.location            = answer{2};
        Gparams.phase               = answer{3};
        Gparams.response_order      = answer{4};
end

Gparams.datetime = datestr(now);

%% Set up Psychtoolbox Screen and parameters
KbName('UnifyKeyNames');

Gparams.bg = 128; % set background color

allScreens = Screen('Screens');

switch Gparams.location
    case 'office'
        rng('shuffle'); % Randomize seed of random number generator
        Gparams.DistToScreen    = .57;		% in meters
        Gparams.ScreenWidth     = .53;		% in meters
        Gparams.pWindow = Screen('OpenWindow', allScreens(1), Gparams.bg);

        Gparams.DataPath    = ['../Data/S', num2str(Gparams.subject_number) '/'];
        Gparams.LogPath     = [Gparams.DataPath, 'Log/'];
    case 'laptop'
        rng('shuffle'); % Randomize seed of random number generator
        Gparams.DistToScreen    = .57;		% in meters
        Gparams.ScreenWidth     = .286;		% in meters
        Gparams.pWindow = Screen('OpenWindow', allScreens(1), Gparams.bg);

        Gparams.DataPath    = ['../Data/S', num2str(Gparams.subject_number) '/'];
        Gparams.LogPath     = [Gparams.DataPath, 'Log/'];
    case 'cubicle'
        rng('shuffle'); % Randomize seed of random number generator
        Gparams.DistToScreen    = .53;		% in meters
        Gparams.ScreenWidth     = .531;		% in meters
        Gparams.pWindow = Screen('OpenWindow', allScreens(1), Gparams.bg);

        Gparams.DataPath    = ['../Data/S', num2str(Gparams.subject_number) '/'];
        Gparams.LogPath     = [Gparams.DataPath, 'Log/'];
    case 'home'
        rng('shuffle'); % Randomize seed of random number generator
        Gparams.DistToScreen    = .53;		% in meters
        Gparams.ScreenWidth     = .407;		% in meters
        Gparams.pWindow = Screen('OpenWindow', allScreens(2), Gparams.bg);

        Gparams.DataPath    = ['../Data/S', num2str(Gparams.subject_number) '/'];
        Gparams.LogPath     = [Gparams.DataPath, 'Log/'];
end

if (~exist(Gparams.DataPath))
    mkdir(Gparams.DataPath);
end
if (~exist(Gparams.LogPath))
    mkdir(Gparams.LogPath);
end


% general screen parameters
[Gparams.ScreenResX, Gparams.ScreenResY]    = Screen('WindowSize', Gparams.pWindow);
Gparams.flipinterval                        = Screen('GetFlipInterval',Gparams.pWindow);
Screen('BlendFunction', Gparams.pWindow, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

Gparams.FixDotSize = deg2pix(.25, Gparams);
Gparams.flip_lag = .015;
Screen('TextSize',Gparams.pWindow, 16);
Screen('TextFont',Gparams.pWindow,'Helvatica');

% general stimulus parameters
Gparams.stim.imSize             = 9;
Gparams.stim.contrast           = .25;
Gparams.stim.window_sd          = 1.5;
Gparams.stim.f_spat             = 0.33;
Gparams.stim.h_eccentricity     = 10;
Gparams.stim.v_eccentricity     = 5;

% general mask parameters
Gparams.mask.window_sd      = 1.5;
Gparams.mask.sd_smooth      = 0.3;
Gparams.mask.contrast       = .5;

% cue parameters
Gparams.cue.size            = deg2pix(0.4, Gparams);
Gparams.cue.color           = [255 255 255];

% response dial parameters
Gparams.dial.line_width = deg2pix(0.6, Gparams);
Gparams.dial.sd = 1.3;

% task parameters
Gparams.sd2afc.orientation_range        = [-14.5, 14.5]; % Draw orientation of sd stimulus from this range
Gparams.sd2afc.stim_levels              = -9:3:9;
Gparams.inducer2ifc.orientation_range   = [-34.5, 34.5]; % Draw orientation of irrelevant inducer stimulus from this range

% timing information
Gparams.timing.cuedur               = 0.35 - Gparams.flip_lag;
Gparams.timing.inducer_prestim      = 0.35 - Gparams.flip_lag;
Gparams.timing.intertask_fixation   = 1.3 - Gparams.flip_lag;
Gparams.timing.sd2afc_prestim       = 0.35 - Gparams.flip_lag;

Gparams.timing.pre_response_fixation = 0.25 - Gparams.flip_lag;

Gparams.timing.afc_iti                  = 1.5 - Gparams.flip_lag - Gparams.flipinterval; % One frame is spent in the block function
Gparams.timing.estimation_iti           = 2 - Gparams.flip_lag - Gparams.flipinterval;

Gparams.timing.inducer_stimdur      = 0.5 - Gparams.flip_lag;
Gparams.timing.inducer_maskdur      = 1 - Gparams.flip_lag;
Gparams.timing.afcstimdur           = 0.5 - Gparams.flip_lag;
Gparams.timing.afcmaskdur           = 1 - Gparams.flip_lag;

% Save the parameters of the experiment
save([Gparams.DataPath, 'Gparams.mat'], 'Gparams');
