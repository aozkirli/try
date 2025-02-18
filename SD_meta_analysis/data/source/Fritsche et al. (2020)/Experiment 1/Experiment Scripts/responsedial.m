function log = responsedial(cfg, Gparams)

t0 = GetSecs;

log = [];
log.cfg = cfg;
log.d_alpha_function = '@(X) (1./(1+exp(-X/5)) - .5) * (pi/45)';
d_alpha_function = eval(log.d_alpha_function);
log.d_alpha_x_max = 100000;
log.button_CW = 'RightArrow';
log.button_CCW = 'LeftArrow';

% Generate transparency mask
[trans_mask] = transparency_mask(...
    Gparams.dial.sd,...
    Gparams.background,...
    Gparams);

% Stimulus position
yPos = (Gparams.ScreenResY-1)/2;
xPos = (Gparams.ScreenResX - 1)/2 + cfg.position;
[s1, s2, s3] = size(trans_mask);
baseRect = [0 0 s1 s2];
dstRect =  CenterRectOnPointd(baseRect, xPos, yPos);

pTexture = Screen('MakeTexture', Gparams.pWindow, trans_mask);



if (isfield(cfg, 'maxtime'))
	log.maxtime = cfg.maxtime;
else
	log.maxtime = inf;
end
if (isfield(cfg, 'respkey'))
	log.respkey = cfg.respkey;
else
	log.respkey = KbName('Space');
end
if (isfield(cfg, 'beginalpha'))
	log.beginalpha = (mod(cfg.beginalpha + 90, 180) - 90) * (pi/180);
else
	log.beginalpha = rand*pi - pi/2;
end

pos_x = (Gparams.ScreenResX - 1) / 2 + cfg.position;
pos_y = (Gparams.ScreenResY - 1) / 2;

alpha = log.beginalpha;
Keys = [];
d_alpha_x = 0;

button_CW = KbName(log.button_CW);
button_CCW = KbName(log.button_CCW);

while(1)
	if ((GetSecs - t0) >= log.maxtime)
		log.rt = GetSecs - t0;
		log.alpha = nan;
		break;
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
                alpha = alpha - d_alpha_function(d_alpha_x);
                alpha = mod(alpha + pi/2, pi) - pi/2;

                d_alpha_x = d_alpha_x + 1;
                if (d_alpha_x > (log.d_alpha_x_max))
                    d_alpha_x = log.d_alpha_x_max;
                end
            case button_CCW
                alpha = alpha + d_alpha_function(d_alpha_x);
                alpha = mod(alpha + pi/2, pi) - pi/2;

                d_alpha_x = d_alpha_x + 1;
                if (d_alpha_x > (log.d_alpha_x_max))
                    d_alpha_x = log.d_alpha_x_max;
                end
			case log.respkey
				log.rt = GetSecs - t0;
				log.alpha = alpha * (180/pi);
				break;		
		end
	end

	S = exp(1i*alpha)*Gparams.dial.r_line;
	Screen('DrawLines', Gparams.pWindow, [real([S, -S]); -imag([S, -S])], Gparams.dial.line_width, 255, [pos_x, pos_y], 1);
    Screen('DrawDots', Gparams.pWindow, [(Gparams.ScreenResX - 1)/2, (Gparams.ScreenResY-1)/2], Gparams.FixDotSize, 255, [], 1);
    Screen('DrawTexture', Gparams.pWindow, pTexture,[],dstRect);
	Screen('Flip', Gparams.pWindow);
    WaitSecs(.001);
end

Screen('Close', pTexture);

end
