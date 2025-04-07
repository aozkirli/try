function figure3A()
FontName = 'Arial';
%% Figure 3A: model predictions
figure('Units','normalized','position',[.05 .4 .34 .5]);
addpath('model')
nanunique = @(x) unique(x(~isnan(x)));
Sd          = 9;    % sigma (witdh) of decoding transition distribution, fixed
Pt          = 0.5;  % probability of integration
tau         = 10;    % temporal decay parameter, fixed
nback       = 20;
meas        = 0;
theta       = datasample(0:179,1000)';
opt_scat   = @(delta,sigma) sqrt(1./(2+(delta/sigma).^2).^2.*sigma.^2+(1-1./(2+(delta/sigma).^2)).^2.*sigma^2);
bincode = {'iso' 'mid' 'ortho'};
sds         = 5:1:20;
colors = blue2pink(length(sds));

for k = 1:numel(sds)
    % bayes
    theta       = datasample(0:179,1000000)';
    Se          = sds(k);
    [model]     = SD_ma_model_bayesian(Se,Sd,Pt,tau,theta,nback,meas);
    delta       = nanunique(model.delta);

    % optimal int
    sigma       = sds(k);
    wp          = 1./(2+(delta/sigma).^2);
    bias        = wp.*delta;
    sigma_dec   = sqrt(wp.^2.*sigma.^2+(1-wp).^2.*sigma^2);


    % plot bayes
    subplot(223)
    model.error(model.delta<0 & model.delta~=-90) = -model.error(model.delta<0 & model.delta~=-90);
    plot(unique(abs(delta)),grpstats(model.error,abs(model.delta)),'-','linewidth',2,'color',colors(k,:));hold on
    g           = grpstats(model.sigma,abs(model.delta),'mean'); % based on posterior

    subplot(224)
    plot(unique(abs(delta)),g./Se,'LineWidth',2,'Color',colors(k,:));hold on

    % plot optInt
    subplot(221)
    plot(delta(delta>=0),bias(delta>=0),'LineWidth',2,'Color',colors(k,:));
    hold on
    subplot(222)
    plot(delta(delta>=0),(sigma_dec(delta>=0))./sigma,'LineWidth',2,'Color',colors(k,:));
    hold on
end

subplot(221)
xticks(0:45:90);
xlim([-5 95]);xticklabels(bincode)
ylim([-1 10])
ylabel('Bias (°)')
xlabel('|Δ| (°)')
box on;
grid on;
set(gca,'FontName',FontName,'FontSize',16,'GridLineWidth', 1)

subplot(222)
plot([-90 90],[sigma sigma]./sigma,'LineWidth',2,'Color','k','LineStyle','--')
xticks(0:45:90);
xlim([-5 95]);xticklabels(bincode)
ylim([0.65 1.2])
ylabel('Normalized scatter (a.u.)')
xlabel('|Δ| (°)')
box on;
grid on;
set(gca,'FontName',FontName,'FontSize',16,'GridLineWidth', 1)

subplot(223)
xticks(0:45:90);
xlim([-5 95]);xticklabels(bincode)
ylim([-1 10])
ylabel('Bias (°)')
xlabel('|Δ| (°)')
box on;
grid on;
set(gca,'FontName',FontName,'FontSize',16,'GridLineWidth', 1)

subplot(224)
plot([-90 90],[sigma sigma]./sigma,'LineWidth',2,'Color','k','LineStyle','--')
xticks(0:45:90);
xlim([-5 95]);xticklabels(bincode)
ylim([0.65 1.2])
ylabel('Normalized scatter (a.u.)')
xlabel('|Δ| (°)')
box on;
grid on;
set(gca,'FontName',FontName,'FontSize',16,'GridLineWidth', 1)

colormap(colors);clrbr =colorbar;
% Customize the tick labels
clrbr.Ticks = [0:1]; % Specify tick positions
clrbr.TickLabels = string(round([sds(1); sds(end)],2))'; % Set custom tick labels

% Adjust colorbar size and position (optional, if needed)
clrbr.Position = [0.95, 0.565, 0.03, 0.35]; % [x, y, width, height]

% Add text on top of the colorbar using normalized coordinates
annotation('textbox', [clrbr.Position(1)-0.02, clrbr.Position(2) + clrbr.Position(4) + 0.02, 0.05, 0.05], ...
    'String', 'Baseline scatter', 'HorizontalAlignment', 'center', ...
    'FontName', FontName, 'FontSize', 16, 'FontWeight', 'bold', 'EdgeColor', 'none', 'Units', 'normalized');

set(gcf,'PaperOrientation','landscape')
set(gcf, 'PaperUnits', 'normalized');
exportgraphics(gcf,['figures' filesep 'main_model_predictions.pdf'],'BackgroundColor','none')
savefig(gcf, ['figures' filesep 'main_model_predictions.fig']);

end