function log = afc_response

t0 = GetSecs;

log = [];
log.button_right = 'RightArrow';
log.button_left  = 'LeftArrow';

Keys = [];

button_right = KbName(log.button_right);
button_left = KbName(log.button_left);

while(1)
    
    oldKeys = Keys;
	[~, ~, keyCode] = KbCheck;
	Keys = find(keyCode, 1);
    
    if(~isempty(Keys))
        
        switch Keys
            case button_right
                log.response = 1;
                log.rt = GetSecs - t0;
                break;
            case button_left
                log.response = -1;
                log.rt = GetSecs - t0;
                break;
        end
    end
end
  
end

