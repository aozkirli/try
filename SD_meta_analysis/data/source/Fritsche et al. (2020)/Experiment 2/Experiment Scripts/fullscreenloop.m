function fullscreenloop(Gparams)

flag = 0;
keycode = 1;

while(find(keycode) ~= KbName('C'))
    switch flag
        case 0
            text = 'Check';
            flag = 1;
        case 1
            text = 'Test';
            flag = 0;
    end
    DrawFormattedText(Gparams.pWindow, text, 'center', 'center', 255);
    Screen('Flip', Gparams.pWindow);
    
    KbReleaseWait;
    [~, keycode] = KbWait;
end

end
