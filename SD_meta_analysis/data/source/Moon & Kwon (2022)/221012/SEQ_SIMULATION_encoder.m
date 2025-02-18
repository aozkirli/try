
%% :: PLOT :: General Settings
set(groot,'DefaultFigureColor', 'w')
set(groot,'DefaultAxesLineWidth', 1)
set(groot,'DefaultAxesXColor', 'k')
set(groot,'DefaultAxesYColor', 'k')
set(groot,'DefaultAxesFontUnits', 'points')
set(groot,'DefaultAxesFontSize', 16) % for ppt: 18, for publication: 25
set(groot,'DefaultAxesFontName', 'Helvetica')
set(groot,'DefaultLineLineWidth', 2/3)
set(groot,'DefaultTextFontUnits', 'Points')
set(groot,'DefaultTextFontSize', 16)
set(groot,'DefaultTextFontName', 'Helvetica')
set(groot,'DefaultAxesBox', 'off')
set(groot,'DefaultAxesTickLength', [0.03 0.025]);
set(groot,'DefaultAxesTickDir','out');
set(groot,'DefaultAxesTickDirMode','manual');
set(groot,'DefaultAxesLabelFontSizeMultiplier',1);


%%
for k = 1:1
clearvars -except k
rng default

%% Parameter manipulations
% [stdProp, pSame, amp0, adaptFactor, stdAdapt]
if k == 1
    par = [15 0 4 .35 30]; % pSame = .9;
%     par = [15 0 4 .9 30]; % pSame = .9;
else
    par = [15 0 4 .35 30]; % model: [15 0 4 .35 20]; current plot: [15 0 4 .35 30];
    par = [15 0 5 0 30]; % model: [15 0 4 .35 20]; current plot: [15 0 4 .35 30];
end
th_t = 210;

%% Model Parameters
% Changes in gain: [ 5~10 0.95 5 20 0.5 30 ]
stdProp     = deg2rad(par(1));
pSame       = par(2); % 0 = uniform, 1 = vm
amp0        = par(3);
adaptFactor = par(4); % 0 = no adaptaion, 1 = large adaptation
stdAdapt    = deg2rad(par(5));

%% Indices
th          = deg2rad(0:359);
[TH,M]      = meshgrid(th,th);
if k == 1
%     nTrials = 1000;
%     thStim  = deg2rad(th_t)*ones(1,nTrials);
    nTrials = 2000;
    thStim  = repmat(deg2rad([th_t 180]),1,nTrials/2);
else
    nTrials = 2000;
    thStim  = repmat(deg2rad([th_t 180]),1,nTrials/2);
end
prGth       = zeros(nTrials,length(th));
EthGr       = zeros(nTrials,1);
% encoder population
respFun     = @(th,thPreferred,a,stdSensory) a.*exp( 1./stdSensory.^2 .* (cos(th-thPreferred')-1) );
N           = 20;
thPreferred = deg2rad(linspace(-180,180,N+1)); thPreferred = thPreferred(1:end-1);
% thPreferred = deg2rad(180)*ones(1,20);
amp         = amp0*ones(N,1);
stdSensory0 = deg2rad(30);
% stdSensory0 = deg2rad(20);
stdSensory  = stdSensory0;
lambda      = respFun( th, thPreferred, amp0, stdSensory0 );
sqrtFI      = @(amp) sqrt(sum( gradient( respFun(th,thPreferred,amp,stdSensory0) ).^2 ./ respFun(th,thPreferred,amp,stdSensory0) ));


%% Trial loop
for trl = 1:nTrials
% population responses
r           = poissrnd( respFun( thStim(trl), thPreferred, amp(:,trl), stdSensory ) );
% likelihood p(r(t)|th(t))
prGth(trl,:)= prod( exp(-lambda).*lambda.^r./factorial(r) );
prGth(trl,:)= prGth(trl,:)./sum(prGth(trl,:)); % normalize

% adaptation - this will affect the next trial
if mod(trl,2) == 0
if k == 1
amp(:,trl+1)= amp0*1./(1+r).^adaptFactor;
elseif k == 2
% stdSensory  = 1./(1+r).^adapt_factor * stdSensory0;
% amp         = 1./(1+r).^adapt_factor * amp0;
amp(:,trl+1)= amp0*(1-adaptFactor*respFun(thStim(trl),thPreferred,1,stdAdapt));
end
else
amp(:,trl+1)= amp0*ones(N,1);
end

end

%% :: PLOT :: 
figure
subplot(2,2,k); hold on;
lambda_new = respFun( th, thPreferred, mean(amp(:,1:2:end),2), stdSensory0 );
plot(lambda_new','k','linewidth',1.25)
annotation('arrow',.3025*ones(1,2),[.89 .84],'linewidth',1.5,'headstyle','vback2','headwidth',15)
% Setting
ax = gca; pbaspect([2.5 1 1])
ax.FontSize = 16;
ax.XLabel.String = 'Direction (\circ)';
ax.YLabel.String = 'Firing rates';
ax.XLim = [-15 360]; ax.YLim = [-.25 4];
ax.XTick = 0:180:360; ax.YTick = 0;
ax.XTickLabel = {'0',[],'360'};
ax.YTickLabel = '0';
offsetAxes(ax,[0 360],[0 4])

subplot(2,2,2); hold on;
bkgrpu = [0 0 0; 74 155 122; 117 112 179]/255;
if k == 1
%     lik = mean(prGth);
%     stdlik = std(prGth);
    lik = mean(prGth(1:2:end,:));
    stdlik = std(prGth(1:2:end,:));
else
    lik = mean(prGth(1:2:end,:));
    stdlik = std(prGth(1:2:end,:));
end
plot([th_t th_t],[0 lik(th_t)],'k:','linewidth',2)
% plot(rad2deg(th),lik,'color','k','linewidth',2)
% shadedErrorBar(rad2deg(th),lik,stdlik,'lineprops',{'color',bkgrpu(4-k,:),'linewidth',2})
shadedErrorBar(rad2deg(th),lik,stdlik,'lineprops',{'color','k','linewidth',2})
% plot([th_t th_t],[0 1.5*max(lik)],'-','color',[.25 .25 .25],'linewidth',2)
text(240,max(lik),' \itp\rm(\bfr\rm|\theta)')
% Setting
ax = gca; pbaspect([2 1 1])
ax.FontSize = 16;
ax.XLabel.String = 'Direction (\circ)';
ax.YLabel.String = 'Likelihood';
ax.XLim = [-15 360]; ax.YLim = [0 1.5*max(lik)];
ax.XTick = [0 180 th_t 360]; ax.YTick = 0;
% ax.XTickLabel = {'0',[],'\theta_{\itt}','360'};
ax.XTickLabel = {'0',[],[],'360'};
ax.YTickLabel = '0';
offsetAxes(ax,[0 360],[0 1.5*max(lik)])

% subplot(2,2,3); hold on;
% % ls = circ_mean(repmat(th',1,1000),prGth(1:2:end,:)'./sum(prGth(1:2:end,:)'));
% % histogram( rad2deg(ls) )
% [~,idx] = max(prGth(1:2:end,:)');
% mle = rad2deg(th(idx));
% histogram( mle )

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

%% :: FUNCTION :: vm
function p = vm(x,mu,sigma)
kappa = 1./sigma.^2;
p = exp( kappa.*(cos(x-mu)-1) ) ./ (2*pi.*besseli(0,kappa,1));
% p = p/max(p);
end

%% :: FUNCTION :: bound_pi
function vec = bound_pi(vec)
vec( vec > pi ) = vec( vec > pi ) - 2*pi;
vec( vec < -pi ) = vec( vec < -pi ) + 2*pi;
end


