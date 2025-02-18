function [pixels] = deg2pix(alpha, Gparams)
% This function assumes that the angle to be transformed originates at (i.e.
% is not centered around) the center.

pixels = tan(alpha * (pi / 180)) * Gparams.DistToScreen * (Gparams.ScreenResX / Gparams.ScreenWidth);

end


