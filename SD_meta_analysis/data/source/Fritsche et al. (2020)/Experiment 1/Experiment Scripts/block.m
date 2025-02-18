function [dataPath, logPath] = block(funcfg, Gparams)

log = [];
log.cfg = funcfg;

%% Initialize
dataPath = [Gparams.DataPath, 'data_loc_Subj', num2str(Gparams.curSubj), '_block', num2str(funcfg.curBlock)];
while(exist([dataPath, '.mat'], 'file'))
    dataPath = [dataPath, '_1'];
end
dataPath = [dataPath, '.mat'];

logPath = [Gparams.LogPath, 'log_loc_Subj', num2str(Gparams.curSubj), '_block', num2str(funcfg.curBlock)];
while(exist([logPath, '.mat'], 'file'))
    logPath = [logPath, '_1'];
end
logPath = [logPath, '.mat'];



%% Counterbalancing of spatial frequencies

freqSequence = carryoverCounterbalance(2,funcfg.cbOrder,funcfg.reps,0);

log.orientation = nan(size(freqSequence,2), 1);
log.position = nan(size(freqSequence,2), 1);
log.trialfreq = nan(size(freqSequence,2), 1);

data = cell(size(freqSequence,2), 1);
log.trial = cell(size(freqSequence,2), 1);


%% Global loop over trials
for n = 1:size(freqSequence,2)
	
    
    cfg = [];
	cfg.orientation = rand*180 - 90;
    cfg.position = funcfg.position;
    cfg.trialfreq = freqSequence(n);

	log.orientation(n) = cfg.orientation;
    log.position(n) = cfg.position;
    log.trialfreq(n) = cfg.trialfreq;

    
	log.trial{n} = trial(cfg, Gparams);
    
    start = GetSecs;

	data{n} = [ ...
		Gparams.curSubj, ...
		funcfg.curBlock, ...
		n, ...
		log.orientation(n), ...
		log.trial{n}.dial.alpha, ...
		log.trial{n}.dial.rt, ...
        log.trialfreq(n)
	];
		
    save(dataPath, 'data');

    % Interrupt trial?
    [~, ~, keyCode] = KbCheck;
    if (ismember(KbName('F8'), find(keyCode)))
        DrawFormattedText(Gparams.pWindow, 'Experiment temporarily interrupted by researcher.\n\nPlease wait...', 'center', 'center', 255);
        Screen('Flip', Gparams.pWindow);
        
        while(KbCheck); 
            WaitSecs(0.01);
        end;
           
        keyCode = 0;
        while(~ismember(KbName('F8'), find(keyCode)))
            [~, ~, keyCode] = KbCheck;
            WaitSecs(0.01);
        end
        
        Screen('Flip', Gparams.pWindow);

        WaitSecs(2);
    end
    
    stop = GetSecs;
    disp(stop-start)
end

save(logPath, 'log');


end
