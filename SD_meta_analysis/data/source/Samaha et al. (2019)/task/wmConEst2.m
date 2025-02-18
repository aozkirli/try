function myscreen = wmConEst2(thresh)
%addpath(genpath('/Users/localadmin/Documents/MATLAB/ecog_task/tasks/mgl'));
%addpath(genpath('/Users/localadmin/Documents/MATLAB/wmConEst/circStats'));

dbclear all
clear global
commandwindow
if nargin == 0
    error('YO DOG! you forgot to put in a threshold')
end

% init screen and open up the window
myscreen.background = [0.5 0.5 0.5];
myscreen.screenNumber = 1;  % zero for small window (debugging), 1 for main display, 2 for secondary disp
myscreen.keyboard.nums = [50 19 20 21 22]; %left and right arrow, respectively 
myscreen.saveData = 1; % use -1 to have it ask to save, 1 for auto save
myscreen.displaySize = [52 32.5]; %WxH in cm
myscreen.displayDistance = 62; % distance from chin rest to screen (cm)
myscreen.eatKeys=1; %stop response keysfrom printing to command window
myscreen.resolution = mglResolution;
myscreen = initScreen(myscreen);
mglDisplayCursor(0); %use '1' to display curser 

% initalize the task params
task{1}.waitForBacktick = 1;
task{1}.segmin =      [0.3 0.3 0.033 nan inf inf 0.05]; % iti, warning, target, resp, confresp, response delay
task{1}.segmax =      [0.5 0.3 0.033 nan inf inf 0.05];
task{1}.getResponse = [0.0 0.0 0.000 0.0 1.0 1.0 0.00];
task{1}.segdur{4} =   .5;%linspace(.6, 12, 5); %delay period randomly decided from  [0.60 3.45 6.30 9.15 12.00] secs
task{1}.numTrials = 60;
%task{1}.randVars.uniform.orientation = 1:180; %in degrees
task{1}.randVars.uniform.posev = [1 .5]; %1 = 1 posevidence (hi conf), .5 = low conf; decrese pos evidence/noise by factor of 1 (half)
task{1}.randVars.uniform.delay = linspace(.6, 12, 5); %delay period randomly decided from  [0.60 3.45 6.30 9.15 12.00] secs
task{1}.randVars.len_ = task{1}.numTrials;
task{1}.randVars.calculated.resp = nan;
task{1}.randVars.calculated.dist = nan; %negative is ccw
task{1}.randVars.calculated.conf = nan;
task{1}.randVars.calculated.orientation = nan;
task{1}.randVars.calculated.targ = nan;
task{1}.randVars.calculated.delay = nan;


% initialize the task
phaseNum = 1;
[task{phaseNum} myscreen] = initTask(task{phaseNum},...
    myscreen,...
    @startSegmentCallback,...
    @screenUpdateCallback,...
    @responseCallback);

% init the stimulus
global stimulus;
myscreen = initStimulus('stimulus',myscreen);
stimulus.thresh =thresh;

%%%%%%%%%%%%%%%%%%%%%%
% instruction screen %
%%%%%%%%%%%%%%%%%%%%%%
mglClearScreen(0.5);mglFlush
mglClearScreen(0.5);mglFlush
mglStrokeText('Now you have to remember the oriention for a variable delay, up to 12 seconds',-11,3,.3,.3);
mglStrokeText('After giving your response with the mouse, please rate how close you were',-11,1,.3,.3);
mglStrokeText('1 = complete guess',-11,0,.3,.3);
mglStrokeText('4 = very close the the true orientation',-11,-1,.3,.3);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main display loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while (phaseNum <= length(task)) && ~myscreen.userHitEsc
  % update the task
  [task myscreen phaseNum] = updateTask(task,myscreen,phaseNum);
  % flip screen
  myscreen = tickScreen(myscreen,task);
end

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task, myscreen)

global stimulus;
mglClearScreen;
if (task.thistrial.thisseg==1) % make stimulus with threshold during ITI
    %define arperture
    mglStencilCreateBegin(1)
    mglFillOval(0,0,[2 2])
    mglStencilCreateEnd;
    mglStencilSelect(1)
    
%     noisefactor=1/task.thistrial.posev;
%     thresh = stimulus.thresh *noisefactor; 
%     gaussian = mglMakeGaussian(3,3,.5,.5);
%     gabor = 256*(1+  mglMakeGrating(3,3,1,task.thistrial.orientation,0).*gaussian  /thresh)/2;
%     noise = 256*(1+(  -1 + (1+1).*rand(size(gabor))).*gaussian  /noisefactor)/2;
%     targ = mean(cat(3,noise,gabor),3);
%     stimulus.targ = mglCreateTexture(targ);
    
    %no gauss
    task.thistrial.orientation = 180;%randsample(1:180,1);
    task.thistrial.targ = task.thistrial.orientation;
    
    noisefactor=1/task.thistrial.posev;
    thresh = stimulus.thresh *noisefactor; 
    gabor = 256*(1+  mglMakeGrating(3,3,1.5,task.thistrial.orientation,0)  /thresh)/2;
    noise = 256*(1+(  -1 + (1+1).*rand(size(gabor)))  /noisefactor)/2;
    targ = mean(cat(3,noise,gabor),3);
    stimulus.targ = mglCreateTexture(targ);   
    
    %enable mouse moving
    robot = java.awt.Robot;
    robot.mouseMove(700, randsample((myscreen.resolution.screenHeight/2)-180:(myscreen.resolution.screenHeight/2)+180,1));
    
end

if (task.thistrial.thisseg==4) % record delay length
	task.thistrial.delay = task.thistrial.seglen(4);
end


if (task.thistrial.thisseg==6) % clear texture to free up mem duing RESP
    mglDeleteTexture(stimulus.targ);
    mglStencilSelect(0)

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = screenUpdateCallback(task, myscreen)

global stimulus

mglClearScreen;
if (task.thistrial.thisseg == 1) % iti
    mglGluDisk(0,0,.08,[.8 .8 .8],24,2);
end

if (task.thistrial.thisseg == 2) % cue
    mglGluDisk(0,0,.08,[.6 .6 .6],24,2);
end

if (task.thistrial.thisseg == 3) % target
    mglGluDisk(0,0,.08,[.6 .6 .6],24,2);
    mglBltTexture(stimulus.targ,[0 0]);
end

if (task.thistrial.thisseg == 4) % blank
    mglGluDisk(0,0,.08,[.6 .6 .6],24,2);
end
 
if (task.thistrial.thisseg == 5) % resp
    mouse = mglGetMouse;
    stimulus.resp = wrapTo360(mouse.y);
        %stimulus.resp = 225;

    probe = 256*(1+  mglMakeGrating(3,3,1.5,stimulus.resp,0)  /30)/2;
    stimulus.probe = mglCreateTexture(probe);
    %probe = zeros(121)*255/2;
    %probe(61,:) = 200;
    %stimulus.probe = mglCreateTexture(probe);
    mglBltTexture(stimulus.probe,[0 0],0,0,stimulus.resp); %THIS STUPID LINE IS DOUBLING THE ORIENTATION!!
    mglGluDisk(0,0,.08,[.6 .6 .6],24,2);

end

if (task.thistrial.thisseg == 6) % confresp
    mglGluDisk(0,0,.08,[.6 .6 .6],24,2);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%
%    responseCallback    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = responseCallback(task,myscreen)
global stimulus
%record resp to randVars 
if (task.thistrial.thisseg == 5) % resp
    if task.thistrial.whichButton == 1 
        dist = rad2deg(circ_dist2(deg2rad(task.thistrial.orientation),[deg2rad(stimulus.resp),deg2rad(stimulus.resp+180)]));
        [v,p] = min(abs(dist));
        task.thistrial.dist = dist(p);
        task.thistrial.resp = stimulus.resp;
        task = jumpSegment(task); 
    end
end

if (task.thistrial.thisseg == 6) % conf resp
    if task.thistrial.whichButton == 2 || task.thistrial.whichButton == 3 || task.thistrial.whichButton == 4 || task.thistrial.whichButton == 5
       task.thistrial.conf=task.thistrial.whichButton-1;
       task = jumpSegment(task); 
    end
end


