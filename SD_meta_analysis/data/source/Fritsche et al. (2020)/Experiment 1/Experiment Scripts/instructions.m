function instructions(Gparams)

showtext('Welcome to the experiment!\n\n\nIn this experiment you will be presented with\npatterns that have a particular orientation.\nIt will be your task to reproduce this orientation.\n\nLet''s look at an example first.\n\n\nPress a key for the next instructions...', Gparams);


% Generate example stimulus
[im, log.stim] = grating(...
    45,...
    Gparams.stim.f_spat1,...
    0.25,...
    Gparams.stim.sigma,...
    Gparams.stim.color,...
    [0 0 0],...
    Gparams.stim.michelson,...
    Gparams.background,...
    Gparams);


% Stimulus position
yPos = (Gparams.ScreenResY-1)/2;
xPos = (Gparams.ScreenResX - 1)/2 + deg2pix(Gparams.stim.eccentricity,Gparams);
[s1, s2, s3] = size(im);
baseRect = [0 0 s1 s2];
dstRect =  CenterRectOnPointd(baseRect, xPos, yPos);



% Display stimulus
pTexture = Screen('MakeTexture', Gparams.pWindow, im);
Screen('DrawTexture', Gparams.pWindow, pTexture,[],dstRect);
Screen('DrawDots', Gparams.pWindow, [(Gparams.ScreenResX - 1)/2, (Gparams.ScreenResY-1)/2], Gparams.FixDotSize, 255, [], 1);
DrawFormattedText(Gparams.pWindow, ...
    'This is an example of a pattern.\n\nThe pattern will be displayed only shortly, after which\nit is replaced by a pattern without any orientation.\n\n Within each block of the experiment the patterns will be presented either right or left of the central fixation point.\nBefore each block you will be informed on which side the patterns will appear. \n\n\nPress a key for the next instructions...' ...
, 'center', 1, 255);
Screen('Flip', Gparams.pWindow);
KbReleaseWait;
KbWait;


% Generate transparency mask
[trans_mask] = transparency_mask(...
    Gparams.dial.sd,...
    Gparams.background,...
    Gparams);

% Mask position
yPos = (Gparams.ScreenResY-1)/2;
xPos = (Gparams.ScreenResX - 1)/2 + deg2pix(Gparams.stim.eccentricity,Gparams);
[s1, s2, s3] = size(trans_mask);
baseRect = [0 0 s1 s2];
dstRect =  CenterRectOnPointd(baseRect, xPos, yPos);


pTexture = Screen('MakeTexture', Gparams.pWindow, trans_mask);


S = exp(1i*(-25 * (pi/180)))*Gparams.dial.r_line;
Screen('DrawLines', Gparams.pWindow, [real([S, -S]); -imag([S, -S])], Gparams.dial.line_width, 255, [(Gparams.ScreenResX - 1)/2 + deg2pix(Gparams.stim.eccentricity,Gparams), (Gparams.ScreenResY-1)/2], 1);
Screen('DrawDots', Gparams.pWindow, [(Gparams.ScreenResX - 1)/2, (Gparams.ScreenResY-1)/2], Gparams.FixDotSize, 255, [], 1);
Screen('DrawTexture', Gparams.pWindow, pTexture,[],dstRect);
DrawFormattedText(Gparams.pWindow, ...
    'Then, a small line will be presented, such as below.\n\nYou will be able to rotate the line using the left and the right arrow key.\nYour job is to rotate the line such that it has the same orientation\nas the pattern that you just saw.\nWhen you are ready, you can press the space bar to continue.\n\nWe''ll practice this in a minute.\n\n\nPress a button for the next instructions...' ...
, 'center', 1, 255);
Screen('Flip', Gparams.pWindow);
KbReleaseWait;
KbWait;

cfg = [];

% position of stimuli in the practice block
if mod(Gparams.curSubj,2) == 0
    cfg.position = -deg2pix(Gparams.stim.eccentricity,Gparams); 
    side = 'left';
else 
    cfg.position = deg2pix(Gparams.stim.eccentricity,Gparams);
    side = 'right';
end

showtext(['In summary:\n\n\n- Observe the orientation of the pattern\n\n- Adjust the line such that it has the same orientation\nas the pattern, using the arrow keys\n\n- When you are finished rotating the line, you press the space bar to continue\n\n\n Importantly, please maintain fixation on the central point at all times it is visible \n\n\nWe can now begin a couple of practice trials.\nIf you have any questions, please ask the researcher now.\nOtherwise, you can proceed.\n\n\nIn the practice trials the patterns will be presented on the ' side ' side\n\nPress a key to begin the practice trials...'], Gparams);

% Practice trials
cfg.curBlock = 1;
cfg.reps = 2;
cfg.cbOrder = 2;

block(cfg, Gparams);


% Check accuracy and precision of responses in the practice block

tmp = load(['Data - Subj' ...
    num2str(Gparams.curSubj) '/data_loc_Subj' num2str(Gparams.curSubj) '_block1.mat']);

practice_data = cell2mat(tmp.data);
practice_data(:,4) = mod(-practice_data(:, 4) + 180,180) - 90;

diff_angles = mod(practice_data(1:end, 5) - practice_data(1:end, 4) + 90, 180) - 90;

mean_error = mean(abs(diff_angles));

disp(['Mean error: ' num2str(mean_error)])

if (mean_error > 30)
    
    Screen('CloseAll');
    ShowCursor;
    ListenChar(0);
    Priority(0);
    
    disp('Responses were too inaccurate!')
    
    error('Ended script because responses were too inaccurate')
    
end


showtext('Well done!\n\nYou are now ready to begin the actual experiment.\nThe experiment consists of 10 blocks.\nAfter each block you can take a break, if you wish.\nThe entire experiment will take about 1.5 hours.\n\nFinally, please keep your eyes fixated on the central fixation dot.\n\nIf you have any questions, please ask the researcher now.\nOtherwise, you can begin.\n\n\nPress a button to begin the real experiment...', Gparams);

end