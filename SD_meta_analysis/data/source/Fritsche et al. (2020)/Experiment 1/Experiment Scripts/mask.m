function [image] = mask(sigma,sd_smooth,color1,color2,background,Gparams)

imSize = 8*sigma;

imSize_px = round(deg2pix(imSize,Gparams));
sigma_px = round(deg2pix(sigma,Gparams));
sd_smooth_px  = round(deg2pix(sd_smooth,Gparams));



trim = .005; % trim off gaussian values smaller than this



noisepatch = rand(imSize_px);


X = 1:imSize_px;
X0 = (X / imSize_px) - .5;  % rescale X -> -.5 to .5
[x_kernel y_kernel] = meshgrid(X0, X0);     % 2D matrices


s = sd_smooth_px / imSize_px;

kernel = exp(-(x_kernel.^2 + y_kernel.^2)/(2*s^2))/(2*pi*s^2);

ftkernel = fft2(kernel);
ftnoise = fft2(noisepatch);
ftfilterednoise = ftkernel .* ftnoise;

filterednoise = real(ifft2(ftfilterednoise));
filterednoise = (filterednoise - min(filterednoise(:))) / range(filterednoise(:)) * 2 - 1;


% Gaussian window
s = sigma_px / imSize_px;
gauss = exp( -(((x_kernel.^2)+(y_kernel.^2)) ./ (2* s.^2)) ); % formula for 2D gaussian
gauss(gauss < trim) = 0; 

filterednoise = (filterednoise+1)./2; % rescale between 0 and 1

r = color1(1) * filterednoise + color2(1) * (1-filterednoise);
g = color1(2) * filterednoise + color2(2) * (1-filterednoise);
b = color1(3) * filterednoise + color2(3) * (1-filterednoise);

r = round(r.* gauss + background(1)*(1-gauss));
g = round(g.* gauss + background(2)*(1-gauss));
b = round(b.* gauss + background(3)*(1-gauss));

image = cat(3,r,g,b);

