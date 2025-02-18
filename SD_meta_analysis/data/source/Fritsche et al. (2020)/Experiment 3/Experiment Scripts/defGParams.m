function Gparams = defGParams(Gparams)

%% Get information about experiment
prompt = {'Enter subject number:','Enter session number','Enter location:'};
dlg_title = 'Experiment Information';
num_lines = 1;
defaultans = {'','','cubicle'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

Gparams.subject_number      = str2num(answer{1});
Gparams.session_number      = str2num(answer{2});
Gparams.location            = answer{3};

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
Gparams.stim.h_eccentricity     = 6.5;
Gparams.stim.v_eccentricity     = 0;

% general mask parameters
Gparams.mask.window_sd      = 1.5;
Gparams.mask.sd_smooth      = 0.3;
Gparams.mask.contrast       = .5;

% response dial parameters
Gparams.dial.line_width = deg2pix(0.6, Gparams);
Gparams.dial.sd = 1.3;

% timing information
Gparams.timing.stimdur      = 0.25 - Gparams.flip_lag;
Gparams.timing.maskdur      = 0.25 - Gparams.flip_lag;

Gparams.timing.cuedur           = 0.25 - Gparams.flip_lag;
Gparams.timing.prestim_fixation = 0.25 - Gparams.flip_lag;

Gparams.timing.response_delay       = [0.05, 3.5] - Gparams.flip_lag;
Gparams.timing.response_timeout     = 3;
Gparams.timing.warning              = 0.25 - Gparams.flip_lag;

Gparams.timing.response_interval = Gparams.timing.response_delay(2) + Gparams.timing.response_timeout + Gparams.timing.warning + 0.25;

Gparams.timing.iti_half    = 0.25 - 2*Gparams.flipinterval;


% Save the parameters of the experiment
save([Gparams.DataPath, 'Gparams_Session_' num2str(Gparams.session_number) '.mat'], 'Gparams');
