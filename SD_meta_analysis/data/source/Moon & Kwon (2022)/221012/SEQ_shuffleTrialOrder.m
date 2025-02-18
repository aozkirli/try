clearvars;
% close all;
rng('default')
plot_or_not = true;

%% Shuffle Once
subjs = 1:32;
arrayIndv = cell(32,1);
rawdataIndv = cell(32,1);
dataIndv = cell(32,2);
% Loop through subjects
for s = subjs
    %% Load Data
    rawdata = []; nbackdata = [];
    for b = 1:15
        load(['E_subj' num2str(s) '_b' num2str(b) '.mat'])
        % Shuffle trial order for inspection
        array = array( randperm(size(array,1)), : );
        rawdata(length(array)*(b-1)+1:length(array)*b,:,1) = array;
    end
    arrayIndv{s} = array;
    rawdataIndv{s} = rawdata;
end

%% Analyse the Shuffled Data with/without the Cardinal Correction
load SEQ_GPRmodels.mat
% For all analyses
binWidth = 10;
xBin = -90:binWidth:90;
% For marginal bias plot
xMoving = -90:90;
MarginalRespBias = nan(length(subjs),length(xMoving),2);
for k = 1:2
    for s = subjs
        %% Cardinal Bias Correction
        array = arrayIndv{s};
        currdata = rawdataIndv{s};
        switch k
            case 1
                % Mean centering without cardinal correction
                currdata(:,4) = bound_180( currdata(:,4) - c_mean(currdata(:,4)) );
            case 2
                % Cardinal correction using GP regression
                currdata(:,4) = bound_180( currdata(:,4) - predict( gprMdl{s}, currdata(:,2) ) );
        end
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
        outliers = abs( nbackdata(:,4) - c_mean(nbackdata(:,4)) ) > 2.5*c_std(nbackdata(:,4));
        data = nbackdata( ~outliers & ~([false(nBack,1); outliers(1:end-nBack)]), : );
        nExcludedTrials(s) = length(nbackdata) - length(data);
        % Finalize data
        dataIndv{s,k} = data;
    end


    %% Analysis & Plot
    for s = subjs
        % Load either empirical or simulation data
        data = dataIndv{s,k};
        %% Marginal Bias
        % Marginal response bias
        for j = 1:length(xMoving)
            l = abs( data(:,7) - xMoving(j) ) <= binWidth/2;
            if sum(l) > 10
                MarginalRespBias(s,j,k) = c_mean(data(l,4));
            end
        end

    end
end

%% BF
% cd([pwd,'/MCMC/BMCBiol_R1/R_shuffleTrialOrder/rng(1) 0426 (current)'])
% rng('default')
% gammaMode = .5; gammaSD = 2;
% gammaScale = 2*gammaSD.^2 ./ (gammaMode + sqrt( gammaMode.^2 + 4*gammaSD.^2 )); % gammaRate = 1./gammaScale;
% gammaShape = 1 + gammaMode./gammaScale;
% prior_mu = unifrnd(-pi,pi,4000,1);
% prior_sigma = gamrnd(gammaScale,gammaShape,4000,1);
% figure; hold on
% histogram(prior_mu./prior_sigma,-5:.1:5,'normalization','pdf')
% load('par_R_shuffle_uncorrected.mat')
% histogram(samples_all(end-6+1,:)./samples_all(end-3+1,:),-5:.1:5,'normalization','pdf')
% BF10_uncorrected = ksdensity(prior_mu./prior_sigma,0) / ksdensity(samples_all(end-6+1,:)'./samples_all(end-3+1,:)',0);
% load('par_R_shuffle_corrected.mat')
% histogram(samples_all(end-6+1,:)./samples_all(end-3+1,:),-5:.1:5,'normalization','pdf')
% BF10_corrected = ksdensity(prior_mu./prior_sigma,0) / ksdensity(samples_all(end-6+1,:)'./samples_all(end-3+1,:)',0);
% cd ../../..

%% :: PLOT ::

%% :: PLOT :: General Settings
set(groot,'DefaultFigureColor', 'w')
set(groot,'DefaultAxesLineWidth', 1)
set(groot,'DefaultAxesXColor', 'k')
set(groot,'DefaultAxesYColor', 'k')
set(groot,'DefaultAxesFontUnits', 'points')
set(groot,'DefaultAxesFontSize', 8/.45) % current: 15
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

if plot_or_not
    % Axis limits
    XLimCB = [-45 45]; YLimCB = [-10 10]; XTickCB = [-40 40]; YTickCB = [-8 8];
    XLimMB = [-100 100]; YLimMB = YLimCB; XTickMB = [-90 90]; YTickMB = YTickCB; nXTick = 7;

    %% :: Figure S3AB :: Marginal Response Bias with/without the Cardinal Correction
    % Marginal response bias
    figure; 
    for k = 1:2
        subplot(1,2,k); hold on
        plot([XLimMB(1) 0; XTickMB(2) 0],[0 YLimMB(1); 0 YTickMB(2)],'k--','linewidth',2/3);
        meanMRB = c_mean(MarginalRespBias(:,:,k));
        ciMRB = tinv(0.975,sum(~isnan(MarginalRespBias(:,:,k)))-1).*c_std(MarginalRespBias(:,:,k),~isnan(MarginalRespBias(:,:,k)))./sqrt(sum(~isnan(MarginalRespBias(:,:,k))));
        shadedErrorBar(xMoving,meanMRB,ciMRB,'lineprops',{'-','color','k','linewidth',2});
        % Setting
        ax = gca; pbaspect(ax, [1 1.1 1])
        ax.FontSize = 16;
        ax.XLabel.String = ['Relative direction of' newline 'previous response (\circ)'];
        switch k
            case 1
                ax.YLabel.String = 'Response error (\circ)';
            case 2
                ax.YLabel.String = 'Residual response error (\circ)';
        end
        ax.XLim = XLimMB; ax.YLim = YLimMB;
        ax.XTick = linspace(XTickMB(1),XTickMB(2),nXTick); ax.YTick = linspace(YTickMB(1),YTickMB(2),5);
        ax.XTickLabel = {'−90','','','0','','','90'}; ax.YTickLabel = {'−8','−4','    0','4','8'};
        offsetAxes(ax,XTickMB,YTickMB)
    end
    
    %% :: Figure S3C :: Response Model fit with/without the Cardinal Correction
    load par_R_shuffle.mat
    figure; hold on
    subplot(121); hold on
    for k = 1:2
        b = bar(k,par{k}.modeGroup(1));
        b.FaceColor = [.8 .8 .8];
        b.EdgeColor = 'none';
        b.BarWidth = .5;
        b.ShowBaseLine = 'off';
        plot(normrnd(k,.1,32,1),par{k}.modeIndv(:,1),'o','markeredgecolor',.5*ones(1,3),'linewidth',1,'markersize',8)
        errorbar(k,par{k}.modeGroup(1),par{k}.hdi95(1,1)-par{k}.modeGroup(1),par{k}.hdi95(1,2)-par{k}.modeGroup(1),'k','linewidth',3,'capsize',0)
        [~,p] = ttest(par{k}.modeIndv(:,1));
        mysigstar(gca, k, 5.5, p);
    end
    plot([.75 2.25],[0 0],'k-','linewidth',1)
    [~,p] = ttest(par{1}.modeIndv(:,1),par{2}.modeIndv(:,1));
    mysigstar(gca, [1 2], 6.5, p, 'none', .15);
    % Setting
    ax = gca; pbaspect([1,1.35,1]) 
    ax.FontSize = 16;
    ax.YLabel.String = 'Bias (\circ)';
    ax.XTickLabel = {'Uncorrected','Corrected'};
    ax.TickLength = [.03 .025];
    ax.XLim = [.5 2.5]; ax.YLim = [-3.5 7.5];
    ax.XTick = 1:2; ax.YTick = -3:3:6;
    ax.YTickLabel = {'  −3','0','3','6'};
    offsetAxes(ax,1,[ax.YTick(1) ax.YTick(end)])
    
end



%% :: FUNCTION ::

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
        [ypos ypos], '-', 'LineWidth', 1, 'color', 'k');
    
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

fz = 30; fontweight = 'normal';
ypos = ypos-0.025*range(get(gca,'ylim'));
if pval < 1e-3
    txt = '***';
elseif pval < 1e-2
    txt = '**';
elseif pval < 0.05
    txt = '*';
elseif ~isnan(pval)
    % this should be smaller
    txt = 'NS';
    %txt = '';
%     fz = fz * 2/3;
    fz = 16;
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

