function log = responsedial(cfg, Gparams)

t0 = GetSecs;

log = [];
log.cfg = cfg;
log.d_alpha_function = '@(X) (1./(1+exp(-X/15)) - .5) * (pi/45)';
d_alpha_function = eval(log.d_alpha_function);
log.d_alpha_x_max = 25;
log.button_CW = 'RightArrow';
log.button_CCW = 'LeftArrow';


if (isfield(cfg, 'maxtime'))
	log.maxtime = cfg.maxtime;
else
	log.maxtime = inf;
end
if (isfield(cfg, 'fadeout_time'))
	log.fadeout_time = cfg.fadeout_time;
else
	log.fadeout_time = 0;
end


% Generate transparency mask
[trans_mask,~] = transparency_mask(...
    Gparams.dial.sd,...
    Gparams.stim.imSize,...
    Gparams.bg,...
    Gparams);

% Stimulus position
yPos = (Gparams.ScreenResY-1)/2;
xPos = (Gparams.ScreenResX - 1)/2 + cfg.h_pos;
[s1, s2, ~] = size(trans_mask);
baseRect = [0 0 s1 s2];
dstRect =  CenterRectOnPointd(baseRect, xPos, yPos);

maskTexture = Screen('MakeTexture', Gparams.pWindow, trans_mask);

imSize = round(deg2pix(Gparams.stim.imSize,Gparams));

%create bar texture
bar = ones(imSize) * Gparams.bg;
bar(:,round(imSize/2 + 1 - Gparams.dial.line_width/2) : round(imSize/2 + 1 + Gparams.dial.line_width/2)) = 255;
bartex = Screen('MakeTexture', Gparams.pWindow, bar);

srcRect=[0 0 imSize+1 imSize+1];

if (isfield(cfg, 'respkey'))
	log.respkey = cfg.respkey;
else
	log.respkey = KbName('Space');
end
if (isfield(cfg, 'beginalpha'))
	log.beginalpha = (mod(-cfg.beginalpha + 180, 180) - 90) * (pi/180);
else
	log.beginalpha = rand*pi - pi/2;
end


alpha = log.beginalpha;
Keys = [];
d_alpha_x = 0;

button_CW = KbName(log.button_CW);
button_CCW = KbName(log.button_CCW);

% determine parameter for fixation dot fading
fix_dot_greyvalue = linspace(255,180,round(cfg.fadeout_time/Gparams.flipinterval)+1);

while(1)
    
    current_time = GetSecs - t0;
    
    if (current_time >= log.maxtime)
        log.rt = nan;
        log.alpha = mod(alpha * (180/pi) + 90 ,180) - 90;
        Screen('Close', maskTexture);
        Screen('Close', bartex);
        return;
    end
    
   
    oldKeys = Keys;
	[~, ~, keyCode] = KbCheck;
	Keys = find(keyCode, 1);

	if(isempty(Keys))
        d_alpha_x = 0;
    else
        if(oldKeys ~= Keys)
            d_alpha_x = 0;
        end

		switch Keys
            case button_CW
                alpha = alpha + d_alpha_function(d_alpha_x);
                alpha = mod(alpha + pi/2, pi) - pi/2;

                d_alpha_x = d_alpha_x + 1;
                if (d_alpha_x > (log.d_alpha_x_max))
                    d_alpha_x = log.d_alpha_x_max;
                end
            case button_CCW
                alpha = alpha - d_alpha_function(d_alpha_x);
                alpha = mod(alpha + pi/2, pi) - pi/2;

                d_alpha_x = d_alpha_x + 1;
                if (d_alpha_x > (log.d_alpha_x_max))
                    d_alpha_x = log.d_alpha_x_max;
                end
			case log.respkey
				log.rt = GetSecs - t0;
				log.alpha = mod(alpha * (180/pi) + 90 ,180) - 90;
				break;		
		end
    end

    
    %draw bar
    Screen('DrawTexture', Gparams.pWindow, bartex,srcRect,dstRect, alpha*180/pi); %last argument is angle in degrees
 
    Screen('DrawTexture', Gparams.pWindow, maskTexture,srcRect,dstRect,alpha*180/pi);
    
    if (current_time >= log.maxtime - cfg.fadeout_time)
        
        i = fix((current_time - (log.maxtime - cfg.fadeout_time)) / Gparams.flipinterval);
        Screen('DrawDots', Gparams.pWindow, [(Gparams.ScreenResX - 1)/2, (Gparams.ScreenResY-1)/2], Gparams.FixDotSize, fix_dot_greyvalue(i+1), [], 1);
          
    else 
        Screen('DrawDots', Gparams.pWindow, [(Gparams.ScreenResX - 1)/2, (Gparams.ScreenResY-1)/2], Gparams.FixDotSize, 255, [], 1);
    end
    
	Screen('Flip', Gparams.pWindow);
    WaitSecs(.001);
end


Screen('DrawDots', Gparams.pWindow, [(Gparams.ScreenResX - 1)/2, (Gparams.ScreenResY-1)/2], Gparams.FixDotSize, 255, [], 1);
Screen('Flip', Gparams.pWindow);

Screen('DrawDots', Gparams.pWindow, [(Gparams.ScreenResX - 1)/2, (Gparams.ScreenResY-1)/2], Gparams.FixDotSize, 255, [], 1);
Screen('Flip', Gparams.pWindow, t0 + log.maxtime - Gparams.flip_lag);

Screen('Close', maskTexture);
Screen('Close', bartex);

end
