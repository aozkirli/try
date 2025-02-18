clearvars;
% close all;

%% Load Data
rawdataIndv = cell(32,1);
dataIndv = cell(32,1);
parSinu = zeros(32,3);
polyCoeff = zeros(32,11);
gprMdl = cell(1,32);

subjs = 1:32;
for s = subjs
    %% Load Data
    rawdata = []; nbackdata = [];
    for b = 1:15
        load(['E_subj' num2str(s) '_b' num2str(b) '.mat'])
        rawdata(length(array)*(b-1)+1:length(array)*b,:,1) = array;
    end
    rawdataIndv{s} = rawdata;
    %% Cardinal Bias Correction
    % Fit GP model
    currdata = rawdata;
    outliers = abs(bound_180( currdata(:,4) - c_mean(currdata(:,4)) )) > 2.5*c_std(currdata(:,4));
    currdata = currdata( ~outliers, : );
    tic
    parSinu(s,:) = fminsearch(@(par) cal_SSE_mySinusoid(currdata(:,2),currdata(:,4),par), zeros(1,3) );
    polyCoeff(s,:) = polyfit(currdata(:,2),currdata(:,4),10);
    gprMdl{s} = fitrgp( currdata(:,2), currdata(:,4) );
    toc
    % Error correction
    currdata = rawdata;
    currdata(:,4) = bound_180( currdata(:,4) - predict( gprMdl{s}, currdata(:,2) ) );
    % Response correction
    currdata(:,3) = bound_360( currdata(:,2) + currdata(:,4) );
    %% Data Preprocessing & Labeling
    % Calculate relative stimulus, response and error on nth previous trial
    nBack = 1;
    for n = 1:nBack
        currdata(1+n:end,5+2*n-1) = bound_180( currdata(1:end-n,2) - currdata(1+n:end,2) );
        currdata(1+n:end,5+2*n) = bound_180( currdata(1:end-n,3) - currdata(1+n:end,2) );
    end
    % Exclude first n trials
    for b = 1:15
        nbackdata( (length(array)-nBack)*(b-1)+1:(length(array)-nBack)*b,:) = currdata( length(array)*(b-1)+1+nBack:length(array)*b, : );
    end
    % Correct some labels
    if s==4 || s==8
        nbackdata( nbackdata(:,1) == 180, 1 ) = -180;
    end
    %% Outlier correction
    % Exclude trials in which error was further than 2.5 s.d. away from the mean
    outliers = abs(bound_180( nbackdata(:,4) - c_mean(nbackdata(:,4)) )) > 2.5*c_std(nbackdata(:,4));
    data = nbackdata( ~outliers & ~([false(nBack,1); outliers(1:end-nBack)]), : );
    nExcludedTrials(s) = length(nbackdata) - length(data);
    % Finalize data
    dataIndv{s} = data;
end


%% Save GP models
save SEQ_GPRmodels.mat gprMdl

%%
set(groot,'DefaultFigureColor', 'w')
set(groot,'DefaultAxesLineWidth', 1)
set(groot,'DefaultAxesXColor', 'k')
set(groot,'DefaultAxesYColor', 'k')
set(groot,'DefaultAxesFontUnits', 'points')
set(groot,'DefaultAxesFontSize', 16) 
set(groot,'DefaultAxesFontName', 'Helvetica')
set(groot,'DefaultLineLineWidth', 2/3)
set(groot,'DefaultTextFontUnits', 'Points')
set(groot,'DefaultTextFontSize', 15)
set(groot,'DefaultTextFontName', 'Helvetica')
set(groot,'DefaultAxesBox', 'off')
set(groot,'DefaultAxesTickLength', [0.025 0.025]);
set(groot,'DefaultAxesTickDir','out');
set(groot,'DefaultAxesTickDirMode','manual');
set(groot,'DefaultAxesLabelFontSizeMultiplier',1);

%% :: Figure - :: Cardinal Bias Correction -- Individual
xi = 0:360;
rawavg = zeros(length(subjs),length(xi));
rawavgstd = zeros(length(subjs),length(xi));
rawavgci = zeros(length(subjs),length(xi));
avg = zeros(length(subjs),length(xi));
avgstd = zeros(length(subjs),length(xi));
avgci = zeros(length(subjs),length(xi));
sinusoidal_fit = zeros(length(subjs),length(xi));
polynomial_fit = zeros(length(subjs),length(xi));
gp_fit = zeros(length(subjs),length(xi));
cardinal_loss = zeros(length(subjs),3);
for s = subjs
    rawdata = rawdataIndv{s};
    data = dataIndv{s};
    windowWidth = 22.5;
    for i = 1:length(xi)
        l = abs( rawdata(:,2) - xi(i) ) <= windowWidth/2 | ...
            abs( rawdata(:,2) - 360 - xi(i) ) <= windowWidth/2 | ...
            abs( rawdata(:,2) + 360 - xi(i) ) <= windowWidth/2;
        if sum(l) > 0
            rawavg(s,i) = c_mean(rawdata(l,4));
            rawavgstd(s,i) = c_std(rawdata(l,4));
            rawavgci(s,i) = tinv(.975,sum(l)-1)*c_std(rawdata(l,4))/sqrt(sum(l));
        end
        l = abs( data(:,2) - xi(i) ) <= windowWidth/2 | ...
            abs( data(:,2) - 360 - xi(i) ) <= windowWidth/2 | ...
            abs( data(:,2) + 360 - xi(i) ) <= windowWidth/2;
        if sum(l) > 0
            avg(s,i) = c_mean(data(l,4));
            avgstd(s,i) = c_std(data(l,4));
            avgci(s,i) = tinv(.975,sum(l)-1)*c_std(data(l,4))/sqrt(sum(l));
        end
    end
    % Fit
    sinusoidal_fit(s,:) = mySinusoid(xi,parSinu(s,:));
    polynomial_fit(s,:) = polyval(polyCoeff(s,:),xi);
    gp_fit(s,:) = predict(gprMdl{s},xi');
    % Loss
    cardinal_loss(s,1) = sqrt(mean(( rawdata(:,4) - mySinusoid(rawdata(:,2),parSinu(s,:)) ).^2));
    cardinal_loss(s,2) = sqrt(mean(( rawdata(:,4) - polyval(polyCoeff(s,:),rawdata(:,2)) ).^2));
    cardinal_loss(s,3) = sqrt(mean(( rawdata(:,4) - predict(gprMdl{s},rawdata(:,2)) ).^2));
end

%% :: Figure S1 :: Cardinal Bias Correction      
% Representative subjects
eg = [16 12];
figure; orient('landscape')
set(gcf,'PaperUnits','centimeters','PaperSize',[34,21]);
for s = 1:length(eg)
    % Raw data
    subplot(2,3,s); hold on
    plot([-20 360],[0 0],'k--','linewidth',2/3);
    shadedErrorBar(xi,rawavg(eg(s),:),rawavgstd(eg(s),:),'patchsaturation',.1,'lineprops',{'k','linewidth',1})
    p3 = plot(xi,sinusoidal_fit(eg(s),:),'linewidth',2.5,'color','#EDB120');
    p4 = plot(xi,polynomial_fit(eg(s),:),'linewidth',2.5,'color','#77AC30');
    p5 = plot(xi,gp_fit(eg(s),:),'linewidth',2.5,'color','#4DBEEE');
    p6 = plot(xi,rawavg(eg(s),:),'k','linewidth',1);
    % Setting
    ax = gca; pbaspect([1 .75 1])
    ax.FontSize = 16;
    if s == 1, ax.YLabel.String = 'Response error (\circ)'; end
    ax.XLim = [-20 380]; ax.YLim = [-35 35];
    ax.XTick = 0:90:360; ax.YTick = -30:10:30;
    ax.XTickLabel = {'0','90','180','270','360'};
    ax.YTickLabel = {[],'  −20',[],'0',[],'20',[]};
    offsetAxes(ax,0,-30)
    % Legend
    if s == 2
        lgd = legend([p6 p3 p4 p5],{'Data','Sinusoid','Polynomial','GP'},'fontsize',14,'numcolumns',2); legend('boxoff');
        lgd.Position = [0.4175    0.8319    0.2055    0.0773]; % top
    end
    % Title
    text(180,40,['Subject #' num2str(eg(s))],'horizontalalignment','center','fontsize',16)

    % Bias-corrected data
    subplot(2,3,s+3); hold on
    plot([-20 360],[0 0],'k--','linewidth',2/3);
    shadedErrorBar(xi,avg(eg(s),:),avgstd(eg(s),:),'patchsaturation',.1,'lineprops',{'k','linewidth',1})    
    ax = gca; pbaspect([1 .75 1])
    ax.FontSize = 16;
    ax.XLabel.String = 'Stimulus direction (\circ)';
    if s == 1, ax.YLabel.String = 'Residual response error (\circ)'; end
    ax.XLim = [-20 380]; ax.YLim = [-35 35];
    ax.XTick = 0:90:360; ax.YTick = -30:10:30;
    ax.XTickLabel = {'0','90','180','270','360'};
    ax.YTickLabel = {[],'  −20',[],'0',[],'20',[]};
    offsetAxes(ax,0,-30)
end

% Group
% Raw data
subplot(233); hold on
plot([0 360],[0 0],'k--','linewidth',2/3)
shadedErrorBar(xi,c_mean(rawavg),tinv(.975,size(rawavg,1)-1)*c_std(rawavg)/sqrt(size(rawavg,1)),'lineprops',{'k','linewidth',2})
ax = gca; pbaspect([1 .75 1])
ax.FontSize = 16;
ax.XLim = [-20 380]; ax.YLim = [-17.5 17.5];
ax.XTick = 0:90:360; ax.YTick = -15:5:15;
ax.XTickLabel = {'0','90','180','270','360'};
ax.YTickLabel = {[],'  −10',[],'0',[],'10',[]};
offsetAxes(ax,0,-15)
% Title
text(180,40/35*ax.YLim(2),'Group mean','horizontalalignment','center','fontsize',16)

% Bias-corrected data
subplot(236); hold on
plot([0 360],[0 0],'k--','linewidth',2/3)
shadedErrorBar(xi,c_mean(avg),tinv(.975,size(avg,1)-1)*c_std(avg)/sqrt(size(avg,1)),'lineprops',{'k','linewidth',2})    
ax = gca; pbaspect([1 .75 1])
ax.FontSize = 16;
ax.XLabel.String = 'Stimulus direction (\circ)';
ax.XLim = [-20 380]; ax.YLim = [-17.5 17.5];
ax.XTick = 0:90:360; ax.YTick = -15:5:15;
ax.XTickLabel = {'0','90','180','270','360'};
ax.YTickLabel = {[],'  −10',[],'0',[],'10',[]};
offsetAxes(ax,0,-15)


%% :: FUNCTION ::

%% :: FUNCTION :: 
function SSE = cal_SSE_mySinusoid( x, y, par )

yhat = mySinusoid( x, par );

SSE = sum(( y - yhat ).^2);

if abs(par(2)) > 1
    SSE = 1e10;
end
end

function yhat = mySinusoid( x, par )

a = par(1);
v = par(2);
y0 = par(3);

yhat = a.*sign(mod(x,180)-90).*(sind(4*abs(mod(x,180)-90)-asind(v))+v)+y0;

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

%% :: FUNCTION :: bound_180
function vec = bound_180( vec )

vec( vec > 180 ) = vec( vec > 180 ) - 360;
vec( vec < -180 ) = vec( vec < -180 ) + 360;

end

%% :: FUNCTION :: bound_360
function vec = bound_360( vec )

vec( vec > 360 ) = vec( vec > 360 ) - 360;
vec( vec < 0 ) = vec( vec < 0 ) + 360;

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

