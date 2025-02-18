function showtext(text, Gparams)

DrawFormattedText(Gparams.pWindow, text, 'center', 'center', 255);
Screen('Flip', Gparams.pWindow);

KbReleaseWait;
KbWait;

end
