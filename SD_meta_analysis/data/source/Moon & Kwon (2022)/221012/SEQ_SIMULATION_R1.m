clearvars;
% close all;
rng default
plot_or_not = true;

%% Load Data
dataIndv = cell(32,1);
dataCombined = [];
subjs = 1:32;
for s = subjs
    rawdata = []; nbackdata = [];
    for b = 1:15
        load(['E_subj' num2str(s) '_b' num2str(b) '.mat'])
        % [stdProp, pSame, amp0, adaptFactor, stdAdapt]
        par = [15 .9 4 .35 20];
        array(:,4) = IdealObserverModel( array(:,2), par );
        array(:,3) = bound_180( array(:,2) + array(:,4) );

        rawdata(length(array)*(b-1)+1:length(array)*b,:) = array;
    end
    %% Data Preprocessing
    currdata = rawdata;
    % Calculate relative stimulus, response and error on nth previous trial
    nBack = 1;
    for n = 1:nBack
        for i = 1+n:length(currdata)
            currdata(i,5+2*n-1) = bound_180(sum( currdata(i-n+1:i,1) ));
        end
        currdata(1+n:end,5+2*n) = bound_180(currdata(1:end-n,3) - currdata(1+n:end,2));
    end
    currdata(2:length(currdata),end+1) = currdata(1:length(currdata)-1,4);
    % Exclude first n trials
    for b = 1:15
        nbackdata( (length(array)-nBack)*(b-1)+1:(length(array)-nBack)*b,:) = currdata( length(array)*(b-1)+1+nBack:length(array)*b, : );
    end
    % Correct some labels
    if s==4 || s==8
        nbackdata( nbackdata(:,1) == 180, 1 ) = -180;
    end
    % Outlier correction
    % Exclude trials in which error was further than 2.5 s.d. away from the mean
    outliers = abs( nbackdata(:,4) - c_mean(nbackdata(:,4)) ) > 2.5*c_std(nbackdata(:,4));
    data = nbackdata( ~outliers & ~([false(nBack,1); outliers(1:end-nBack)]), : );
    nExcludedTrials(s) = length(nbackdata) - length(data);
    % Finalize data
    dataIndv{s} = data;
    dataCombined(end+1:end+length(data),:) = data;
end

%% Analysis & Plot
% For all analyses
binWidth = 10;
xBin = -90:binWidth:90;
% For joint bias plot
map = nan(length(xBin),length(xBin),length(subjs));
indvmap = nan(length(xBin),length(xBin),length(subjs));
isShortIndv = cell(1,length(subjs));
numEmpty = zeros(1,length(subjs));
numShort = zeros(1,length(subjs));
mapAll = nan(length(xBin));
numEmptyAll = 0;
numShortAll = 0;
[X,Y] = meshgrid(1:length(xBin)); X = X(:); Y = Y(:);
%% Subject Loop
for s = subjs
    % Load either empirical or simulation data
    data = dataIndv{s};
    mapdata = dataCombined;
    %% Joint Bias
    % Individual
    for i = 1:length(xBin)
        for j = 1:length(xBin)
            l = ( data(:,1)==xBin(i) ) & ( data(:,7)>=xBin(j)-binWidth/2 & data(:,7)<xBin(j)+binWidth/2 );
            if sum(l) == 0
                numEmpty(s) = numEmpty(s) + 1;
            elseif sum(l) < 5
                numShort(s) = numShort(s) + 1;
                map(i,j,s) = c_mean(data(l,4));
            else
                indvmap(i,j,s) = c_mean(data(l,4));
                map(i,j,s) = c_mean(data(l,4));
            end
        end
    end
    isShortIndv{s} = rot90(isnan(indvmap(:,:,s)));

    % Mean
    if s == subjs(end)
        indvmap = rot90(indvmap);
        map = rot90(map);
        meanmap = c_mean(map,[],3);
        isShort = sum(~isnan(map),3)<16;
    end

end


%% :: PLOT ::

%% :: PLOT :: General Settings
load bwr.mat
set(groot,'DefaultFigureColormap',bwr)
set(groot,'DefaultFigureColor', 'w')
set(groot,'DefaultAxesLineWidth', 1.5)
set(groot,'DefaultAxesXColor', 'k')
set(groot,'DefaultAxesYColor', 'k')
set(groot,'DefaultAxesFontUnits', 'points')
set(groot,'DefaultAxesFontSize', 20) % for ppt: 18, for publication: 25
set(groot,'DefaultAxesFontName', 'Helvetica')
set(groot,'DefaultLineLineWidth', 1)
set(groot,'DefaultTextFontUnits', 'Points')
set(groot,'DefaultTextFontSize', 20)
set(groot,'DefaultTextFontName', 'Helvetica')
set(groot,'DefaultAxesBox', 'off')
set(groot,'DefaultAxesTickLength', [0.02 0.025]);
set(groot,'DefaultAxesTickDir','out');
set(groot,'DefaultAxesTickDirMode','manual');
set(groot,'DefaultAxesLabelFontSizeMultiplier',1);



if plot_or_not
    %% :: PLOT :: Joing Bias Plot -- Mean
    figure; 
    clim = [-8 8];
    imagesc(meanmap,clim); hold on;
    cbar = colorbar;
    patch((repmat(X(isShort(:)),1,4)+[-.5 .5 .5 -.5])',(repmat(Y(isShort(:)),1,4)+[-.5 -.5 .5 .5])',[.9 .9 .9],'edgecolor',[.9 .9 .9])
    text([4 16],[4 16],'NA','horizontalalignment','center','fontsize',8/.45)
    % Setting
    ax = gca; pbaspect(ax, [1 1 1])
    ax.FontSize = 8/.45;
    ax.XLabel.String = 'Relative direction of previous stimulus (\circ)';
    ax.YLabel.String = 'Relative direction of previous response (\circ)';
    ax.LineWidth = 1; ax.TickDir = 'in'; ax.TickLength = [.01 .025];
    ax.XTick = 1:3:19; ax.XTickLabel = {'−90','−60','−30','0','30','60','90'};
    ax.YTick = 1:3:19; ax.YTickLabel = {'90','60','30','0','−30','−60',' −90'};
    cbar.LineWidth = 1; cbar.Color = 'k'; cbar.FontSize = 8/.45;
    cbar.Ticks = clim(1):range(clim)/4:clim(2); cbar.TickLength = 0;
    cbar.TickLabels = {'−8  ','−4','0','4','8'};
    cbar.Label.String = 'Response error (\circ)';
    cbar.Label.VerticalAlignment = 'middle'; cbar.Label.Rotation = 270;

end


%% :: FUNCTION ::

%% :: FUNCTION :: IdealObserverModel
function Error = IdealObserverModel( thStim, par )
    %% Model Parameters
    stdProp     = deg2rad(par(1));
    pSame       = par(2); % 0 = uniform, 1 = vm
    amp0        = par(3);
    adaptFactor = par(4); % 0 = no adaptaion, 1 = large adaptation
    stdAdapt    = deg2rad(par(5));

    %% Indices
    th          = deg2rad(-180:179);
    [TH,~]      = meshgrid(th,th);
    thStim      = deg2rad(thStim);
    EthGr       = zeros(length(thStim),1);
    % encoder population
    respFun     = @(th,thPreferred,a,stdSensory) a.*exp( 1./stdSensory.^2 .* (cos(th-thPreferred')-1) );
    N           = 20;
    thPreferred = deg2rad(linspace(-180,180,N+1)); thPreferred = thPreferred(1:end-1);
    amp         = amp0;
    stdSensory0 = deg2rad(30);
    stdSensory  = stdSensory0;
    lambda      = respFun( th, thPreferred, amp0, stdSensory0 );
    % propagation noise
    vmProp      = exp(1./stdProp.^2.*(cos(TH-TH')-1))./(2*pi.*besseli(0,1./stdProp.^2,1));
    unifProp    = unifpdf(TH,min(th),max(th));
    pthGthtm1   = pSame * vmProp + (1-pSame) * unifProp;
    pthGthtm1   = pthGthtm1./(sum(pthGthtm1)); % normalize
    % initial prediction (uniform)
    pth         = unifpdf(th,min(th),max(th));
    pth         = pth./sum(pth); % normalize

    %% Trial loop
    for trl = 1:length(thStim)
        % population responses
        r           = poissrnd( respFun( thStim(trl), thPreferred, amp, stdSensory ) );
        % likelihood p(r(t)|th(t))
        prGth       = prod( exp(-lambda).*lambda.^r./factorial(r) );
        prGth       = prGth./(sum(prGth)); % normalize
        % integration
        pthGr       = pth .* prGth;
        pthGr       = pthGr./sum(pthGr);
        % thh
        EthGr(trl)  = c_mean(th,pthGr); % L2 norm

        % adaptation - this will affect the next trial
        amp         = amp0*(1-adaptFactor*respFun(thStim(trl),thPreferred,1,stdAdapt));
        % prediction - this will affect the next trial
        pth         = pthGr * pthGthtm1; % convolve
        pth         = pth./sum(pth); % normalize
    end

    %% Finalize
    Error = rad2deg(bound_pi( EthGr-thStim ));

end

%% :: FUNCTION :: c_mean
function mu = c_mean(alpha, w, dim)

if nargin < 3
  dim = find(size(alpha) > 1, 1, 'first');
  if isempty(dim)
    dim = 1;
  end
end

if nargin < 2 || isempty(w)
  % if no specific weighting has been specified
  % assume no binning has taken place
	w = ones(size(alpha));
else
  if size(w,2) ~= size(alpha,2) || size(w,1) ~= size(alpha,1) 
    error('Input dimensions do not match');
  end 
end

% deg2rad
alpha = deg2rad(alpha);

% compute weighted sum of cos and sin of angles
r = sum(w.*exp(1i*alpha),dim,'omitnan');

% obtain mean by
mu = angle(r);

% rad2deg
mu = rad2deg(mu);

end

%% :: FUNCTION :: c_std
function s0 = c_std(alpha, w, d, dim)

if nargin < 4
  dim = find(size(alpha) > 1, 1, 'first');
  if isempty(dim)
    dim = 1;
  end  
end

if nargin < 3 || isempty(d)
  % per default do not apply correct for binned data
  d = 0;
end

if nargin < 2 || isempty(w)
  % if no specific weighting has been specified
  % assume no binning has taken place
	w = ones(size(alpha));
else
  if size(w,2) ~= size(alpha,2) || size(w,1) ~= size(alpha,1) 
    error('Input dimensions do not match');
  end 
end

% deg2rad
alpha = deg2rad(alpha);

% compute mean resultant vector length
r = circ_r(alpha,w,d,dim);

s0 = sqrt(-2*log(r));   % 26.21

% rad2deg
s0 = rad2deg(s0);

end

%% :: FUNCTION :: circ_r
function r = circ_r(alpha, w, d, dim)

if nargin < 4
  dim = find(size(alpha) > 1, 1, 'first');
  if isempty(dim)
    dim = 1;
  end
end

if nargin < 2 || isempty(w) 
  % if no specific weighting has been specified
  % assume no binning has taken place
	w = ones(size(alpha));
else
  if size(w,2) ~= size(alpha,2) || size(w,1) ~= size(alpha,1) 
    error('Input dimensions do not match');
  end 
end

if nargin < 3 || isempty(d)
  % per default do not apply correct for binned data
  d = 0;
end

% compute weighted sum of cos and sin of angles
r = sum(w.*exp(1i*alpha),dim,'omitnan');

% obtain length 
r = abs(r)./sum(w,dim);

% for data with known spacing, apply correction factor to correct for bias
% in the estimation of r (see Zar, p. 601, equ. 26.16)
if d ~= 0
  c = d/2/sin(d/2);
  r = c*r;
end
end

%% :: FUNCTION :: bound_pi
function vec = bound_pi(vec)

vec( vec > pi ) = vec( vec > pi ) - 2*pi;
vec( vec < -pi ) = vec( vec < -pi ) + 2*pi;

end

%% :: FUNCTION :: bound_180
function vec = bound_180( vec )

vec( vec > 180 ) = vec( vec > 180 ) - 360;
vec( vec < -180 ) = vec( vec < -180 ) + 360;

end

%% :: FUNCTION :: offsetAxes
function offsetAxes( ax, x, y )
if ~exist('ax', 'var'), ax = gca; end
if ~exist('x', 'var')
    x = [ax.XLim(1)+(ax.XTick(2)-ax.XTick(1))/4 ax.XLim(2)-(ax.XTick(2)-ax.XTick(1))/4];
elseif numel(x)<2
    x = [x ax.XLim(2)-(x-ax.XLim(1))];
end
if ~exist('y', 'var')
    y = [ax.YLim(1)+(ax.YTick(2)-ax.YTick(1))/4 ax.YLim(2)-(ax.YTick(2)-ax.YTicm(1))/4];
elseif numel(y)<2
    y = [y ax.YLim(2)-(y-ax.YLim(1))];
end
% remove the x and y ticks if they are outside the visible X and YLim
ax.XTick = ax.XTick( ax.XTick >= x(1) & ax.XTick <= x(2) );
ax.YTick = ax.YTick( ax.YTick >= y(1) & ax.YTick <= y(2) );
% position the X and Y labels on the center of the visible axes
% ax.XLabel.Position(1) = mean(x);
% ax.YLabel.Position(2) = mean(y);
% this will keep the changes constant even when resizing axes
addlistener (ax, 'MarkedClean', @(obj,event)resetVertex(ax,x,y));
end
function resetVertex( ax, x, y )
% extract the x axis vertext data
% X, Y and Z row of the start and end of the individual axle.
ax.XRuler.Axle.VertexData(1,1) = x(1);
ax.XRuler.Axle.VertexData(1,2) = x(2);
% repeat for Y (set 2nd row)
ax.YRuler.Axle.VertexData(2,1) = y(1);
ax.YRuler.Axle.VertexData(2,2) = y(2);
end

%% :: FUNCTION :: mysigstar
function h = mysigstar(ax, xpos, ypos, pval, whichWay, howMuch)
% replaces sigstar, which doesnt work anymore in matlab 2014b

if ~exist('ax', 'var'); ax = gca; end
if ~exist('whichWay', 'var'), whichWay = 'none'; end
if ~exist('howMuch', 'var'), howMuch = .05*range(get(gca,'ylim')); end

if numel(ypos) > 1
    assert(ypos(1) == ypos(2), 'line wont be straight!');
    ypos = ypos(1);
end

% draw line
hold on;
if numel(xpos) > 1
    % plot the horizontal line
    p = plot(ax, [xpos(1), xpos(2)], ...
        [ypos ypos], '-', 'LineWidth', 1.5, 'color', 'k');
    
    % also add small downward ticks
    switch whichWay
        case 'down'
            plot([xpos(1) xpos(1)], [ypos ypos-howMuch], '-k', 'LineWidth', 1.5);
            plot([xpos(2) xpos(2)], [ypos ypos-howMuch], '-k', 'LineWidth', 1.5);
        case 'up'
            plot([xpos(1) xpos(1)], [ypos ypos+howMuch], '-k', 'LineWidth', 1.5);
            plot([xpos(2) xpos(2)], [ypos ypos+howMuch], '-k', 'LineWidth', 1.5);
        case 'none'
    end
    txtBg = 'w';
else
    txtBg = 'none';
end

fz = 33; fontweight = 'normal';
ypos = ypos-0.025*range(get(gca,'ylim'));
if pval < 1e-3
    txt = '***';
elseif pval < 1e-2
    txt = '**';
elseif pval < 0.05
    txt = '*';
elseif ~isnan(pval)
    % this should be smaller
    txt = 'n.s.';
    %txt = '';
%     fz = fz * 2/3;
    fz = 18;
    ypos = ypos+0.025*range(get(gca,'ylim'));
    fontweight = 'normal';
else
    return
end

% draw the stars in the bar
h = text(mean(xpos), mean(ypos), txt, ...
    'horizontalalignment', 'center', 'backgroundcolor', ...
    txtBg, 'margin', 6, 'fontsize', fz, 'fontweight', fontweight, 'Parent', ax);
end