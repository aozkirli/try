clear all;
close all;
clc;
addpath('utils')
load(['..' filesep '..' filesep 'data' filesep 'datasets' filesep 'SD_ma_master_table.mat'])
tbl.theta(tbl.theta==180) = 0;
tbl.delta(tbl.delta==-90) = 90;
nBoot = []; % if empty, gives true mean and non-smoothed confidence intervals
rep_size = 3;
bin_size = 21;
FontName = 'Arial';
%% model predictions:
figure('Units','normalized','position',[.05 .4 .6 .4]);
addpath('model')
nanunique = @(x) unique(x(~isnan(x)));
Sd          = 9;    % sigma (witdh) of decoding transition distribution, fixed
Pt          = 0.5;  % probability of integration
tau         = 1;    % temporal decay parameter, fixed
nback       = 1;
meas        = 0;
theta       = datasample(0:179,10000)';
opt_scat   = @(delta,sigma) sqrt(1./(2+(delta/sigma).^2).^2.*sigma.^2+(1-1./(2+(delta/sigma).^2)).^2.*sigma^2);
bincode = {'iso' 'mid' 'ortho'};
figure('position',[40 180 1000 570])
sds         = 5:1:20;

% Define colors
dark_blue = [0.172, 0, 0.627];  % Dark blue (pastel)
purple = [0.494, 0, 0.835];    % Purple (Pastel)
dark_red = [0.627, 0.129, 0.302]; % Dark red (pastel)

% Define transition points (dark pink -> violet -> purple)
colors = [
    dark_blue;
    purple;
    dark_red;
    ];
% Interpolate to create the colormap
x = linspace(1, size(colors, 1), length(sds)); % Target range
xi = 1:size(colors, 1); % Original range
colors = interp1(xi, colors, x);

for k = 1:numel(sds)
    % bayes
    theta       = datasample(0:179,1000000)';
    Se          = sds(k); 
    Sd          = 9;
    Pt          = .5;
    tau         = 1;
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

    % save bins for later
    orthoB(k)    = mean(model.sigma(ismember(abs(model.delta),90)));
    isoB(k)      = mean(model.sigma(ismember(abs(model.delta), 0)));
    midB(k)      = mean(model.sigma(ismember(abs(model.delta),45)));

    orthoC(k)    = mean(sigma_dec(ismember(abs(delta),90)));
    isoC(k)      = mean(sigma_dec(ismember(abs(delta), 0)));
    midC(k)      = mean(sigma_dec(ismember(abs(delta),45)));
   
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
% Assuming you already have your plot
% clrbr.Location = 'eastoutside'; % Move colorbar to the right of the axis

% Customize the tick labels
clrbr.Ticks = [0:1]; % Specify tick positions
clrbr.TickLabels = string(round([sds(1); sds(end)],2))'; % Set custom tick labels

% Adjust colorbar size and position (optional, if needed)
clrbr.Position = [0.95, 0.565, 0.03, 0.35]; % [x, y, width, height]

% Add text on top of the colorbar using normalized coordinates
annotation('textbox', [clrbr.Position(1)-0.02, clrbr.Position(2) + clrbr.Position(4) + 0.02, 0.05, 0.05], ...
    'String', 'Baseline scatter', 'HorizontalAlignment', 'center', ...
    'FontName', FontName, 'FontSize', 16, 'FontWeight', 'bold', 'EdgeColor', 'none', 'Units', 'normalized');

%% summarize datasets:
datasets     = unique(tbl.studynum);  % Get unique dataset codes
n_datasets   = numel(datasets);
variable     = 'error_ori_deb';  % Variable for which to compute trial counts
report_table = table();  % Initialize report table

% Loop over each dataset and compute trial counts
for i = 1:n_datasets
    tmp          = tbl(tbl.studynum==i, :);  % Subset data for the current dataset

    % Prepare report row with dataset info
    [ID, Study, Datasets, Stimulus, N, AllTrials, CleanTrials] = deal(...
        {num2str(i)}, tmp.study(1), length(unique(tmp.codenum)), tmp.stimulus(1), ...
        numel(unique(tmp.obsid)), length(tmp.(variable)), sum(isfinite(tmp.(variable)) & abs(tmp.delta)<=90));

    % Append report row to report table
    report_i     = table(ID, Study, Stimulus, Datasets,  N, AllTrials, CleanTrials);
    report_table = vertcat(report_table, report_i);
end
[ID, Study, Stimulus, Datasets, N, AllTrials, CleanTrials] = deal({' '},{' '},{'Total'},sum(report_table.Datasets),sum(report_table.N),sum(report_table.AllTrials),sum(report_table.CleanTrials));
report_i     = table(ID, Study, Stimulus, Datasets,  N, AllTrials, CleanTrials);
report_table = vertcat(report_table, report_i);
writetable(report_table,['tables' filesep 'summary_studies.csv'])


%% bias cleaning check
% Check how cleaning went:
figure('Units','normalized','position',[.05 .4 .6 .4]);
tmp = tbl(abs(tbl.delta)<=90 & tbl.theta<=180,:);

% Orientation bias removal
rep_size = 3;
bin_size = 21;

subplot(121)
p0 = plot_mean_ci(tmp.theta,tmp.error_norm,rep_size,bin_size,'k','mean');
hold on;
p1 = plot_mean_ci(tmp.theta,tmp.error_iqr_norm,rep_size,bin_size,'r','mean');
p2 = plot_mean_ci(tmp.theta,tmp.error_ori_deb_norm,rep_size,bin_size,'b','mean');
p3 = plot_mean_ci(tmp.theta,tmp.error_ori_deb_sd_deb_norm,rep_size,bin_size,'g','mean');
ylim([-.5 .5])
xticks(0:45:180)
xlabel('θ(°)')
ylabel('stimulus-specific bias')
set(gca, 'FontSize', 20, 'LineWidth', 2, 'GridLineWidth', 1, 'FontName', FontName);

subplot(122)
p0 = plot_mean_ci(tmp.delta,tmp.error_norm,rep_size,bin_size,'k','mean');
hold on;
p1 = plot_mean_ci(tmp.delta,tmp.error_iqr_norm,rep_size,bin_size,'r','mean');
p2 = plot_mean_ci(tmp.delta,tmp.error_ori_deb_norm,rep_size,bin_size,'b','mean');
p3 = plot_mean_ci(tmp.delta,tmp.error_ori_deb_sd_deb_norm,rep_size,bin_size,'g','mean');
ylim([-.5 .5])
xticks(-90:45:90)
xlabel('Δ (°)')
ylabel('serial dependence bias')
set(gca, 'FontSize', 20, 'LineWidth', 2, 'GridLineWidth', 1, 'FontName', FontName);
legend([p0 p1 p2 p3],{'raw data' 'outliers removed' 'stimulus bias removed' 'serial dependence bias removed'}, 'Location', 'best')

set(gcf,'PaperOrientation','landscape')
set(gcf, 'PaperUnits', 'normalized');
exportgraphics(gcf,['figures' filesep 'SI_removalChecks.pdf'],'BackgroundColor','none')

%% Show ortho vs. iso ori biases
error_variable = 'error_iqr_norm';
n_bins = 10;
intersect_thetas = [0:181];

tmp    = tbl(abs(tbl.delta)<=90 & tbl.theta<=180 & ~isnan(tbl.delta),:);
rep_size = 3;
bin_size = 31;
for i = 0:180
    dd(i+1) = nanmean(tmp.error_iqr_norm(tmp.theta == i))/nanstd(tmp.error_iqr_norm(tmp.theta == i));
end
[abs_dd_sorted,idx]=sort(abs(dd));

% Create a new figure
figure('Units','normalized','position',[0 0 1 1]);
% Define positions for 2x2 on the left-hand side
left_positions = [0 .065 0 0]+[
    0.1, 0.5, 0.2, 0.35;  % Top-left
    0.1, 0.075, 0.2, 0.35;  % Bottom-left
    0.4, 0.5, 0.2, 0.35;  % Top-right
    0.4, 0.075, 0.2, 0.35   % Bottom-right
    ];

% Define colors
dark_blue = [0.172, 0, 0.627];  % Dark blue (pastel)
purple = [0.494, 0, 0.835];    % Purple (Pastel)
dark_red = [0.627, 0.129, 0.302]; % Dark red (pastel)

% Define transition points (dark pink -> violet -> purple)
colors = [
    dark_blue;
    purple;
    dark_red;
    ];
% Interpolate to create the colormap
x = linspace(1, size(colors, 1), n_bins); % Target range
xi = 1:size(colors, 1); % Original range
pastel_cmap = interp1(xi, colors, x);

superiority = nan(n_bins,1);
sd_peak     = nan(n_bins,1);
cohenD_mean = nan(n_bins,1);

% Define positions for the 2x2 grid (left)
for i = 1:4
    subplot('Position', left_positions(i,:));
    for b = 1:n_bins
        highlight_idx = intersect(idx((b-1)*180/n_bins+1:b*180/n_bins),intersect_thetas);
        cohenD_mean(b) = nanmean(abs_dd_sorted((b-1)*180/n_bins+1:b*180/n_bins));
        tt = tmp(ismember(abs(tmp.theta),highlight_idx-1),:);
        switch i
            case 1
                if b ==1
                    p1 = plot_mean_ci(tmp.theta,tmp.(error_variable),rep_size,bin_size,'k','mean',1);
                    xlabel('θ (°)')
                    ylabel('Bias Amplitude')
                    ylim([-.4 .4])
                    xticks(0:45:180)
                    title({'STIMULUS-SPECIFIC','(with respect to θ)'})
                    colormap(pastel_cmap);clrbr =colorbar;
                end
                h1=plot(highlight_idx-1,zeros(size(highlight_idx)),'Marker','square','Color',pastel_cmap(b,:),'MarkerSize',10,'MarkerFaceColor',pastel_cmap(b,:),'LineStyle','none');
            case 2
                if b ==1
                    p1 = plot_mean_ci(tmp.theta,tmp.(error_variable),rep_size,bin_size,'k','std',0);
                    xlabel('θ (°)')
                    ylabel('Error Scatter')
                    ylim([.88 1.12])
                    yticks(.9:.1:1.1)
                    xticks(0:45:180);
                end

            case 3
                p1m = plot_mean_ci(tt.delta,tt.(error_variable),rep_size,bin_size,pastel_cmap(b,:),'mean',1);
                xlabel('Δ (°)')
                ylabel('Bias Amplitude')
                ylim([-.4 .4])
                xticks(-90:45:90);
                sd_peak(b) = max(p1m.YData);
                if b ==1
                    title({'SERIAL DEPENDENCE','(with respect to Δ)'})
                end
            case 4
                [p1s,f1s] = plot_mean_ci(tt.delta,tt.(error_variable),rep_size,bin_size,pastel_cmap(b,:),'std',1);
                xlabel('Δ (°)')
                ylabel('Error Scatter')
                ylim([.88 1.12])
                yticks(.9:.1:1.1)
                xticks(-90:45:90);
                superiority(b) = (p1s.YData(1)-p1s.YData(91));
        end
    end
    if i==4
        [p1,f1]=plot_mean_ci(tmp.delta,tmp.(error_variable),rep_size,bin_size,'k','std',1);
        p1.LineWidth = 4;
    end
end
% Assuming you already have your plot
% clrbr.Location = 'eastoutside'; % Move colorbar to the right of the axis

% Customize the tick labels
clrbr.Ticks = [0; 1]; % Specify tick positions
clrbr.TickLabels = string(round([cohenD_mean(1) cohenD_mean(end)],2))'; % Set custom tick labels

% Adjust colorbar size and position (optional, if needed)
clrbr.Position = [0.33, 0.565, 0.01, 0.35]; % [x, y, width, height]

% Add text on top of the colorbar using normalized coordinates
annotation('textbox', [clrbr.Position(1)-0.02, clrbr.Position(2) + clrbr.Position(4) + 0.02, 0.05, 0.05], ...
    'String', 'Cohen''s d for θs included', 'HorizontalAlignment', 'center', ...
    'FontName', FontName, 'FontSize', 16, 'FontWeight', 'bold', 'EdgeColor', 'none', 'Units', 'normalized');

% Define positions for 3x1 on the right-hand side
right_positions = [
    0.675, 0.75, 0.12, 0.22;  % Top-right
    0.675, 0.425, 0.12, 0.22;  % Middle-right
    0.675, 0.1, 0.12, 0.22   % Bottom-right
    ];

% % Define positions for the 3x1 grid (right)
for i = 1:3
    % Adjust position manually for the 3x1 grid
    subplot('Position', right_positions(i,:));
    for b=1:n_bins
        switch i
            case 1
                scatter(cohenD_mean(b),sd_peak(b),200,'filled','o','MarkerFaceColor',pastel_cmap(b,:),'MarkerFaceAlpha',.8);hold on;
                xlabel({'Stimulus-specific Bias' 'Effect Size'})
                ylabel({'Serial Dependence' 'Bias Peak'})
                ll1 = plot(cohenD_mean, polyval(polyfit(cohenD_mean, sd_peak, 1), cohenD_mean),'Color',[.5 .5 .5], 'LineWidth', 2);
                [rr,pp] = corr(cohenD_mean,sd_peak); legend(ll1,['r=' num2str(round(rr,2)) repmat('*', 1, (pp < 0.001) * 3 + (pp >= 0.001 & pp < 0.01) * 2 + (pp >= 0.01 & pp < 0.05) * 1)],'Location','northeast')
                set(gca, 'FontSize', 20, 'LineWidth', 2, 'GridLineWidth', 1, 'FontName', FontName);
                grid on;
                box on;
            case 2
                scatter(cohenD_mean(b),superiority(b),200,'filled','o','MarkerFaceColor',pastel_cmap(b,:),'MarkerFaceAlpha',.8);hold on;
                xlabel({'Stimulus-specific Bias' 'Effect Size'})
                ylabel({'Superiority' '(ortho-iso)'})
                ll1 = plot(cohenD_mean, polyval(polyfit(cohenD_mean, superiority, 1), cohenD_mean),'Color',[.5 .5 .5], 'LineWidth', 2);
                [rr,pp] = corr(cohenD_mean,superiority); legend(ll1,['r=' num2str(round(rr,2)) repmat('*', 1, (pp < 0.001) * 3 + (pp >= 0.001 & pp < 0.01) * 2 + (pp >= 0.01 & pp < 0.05) * 1)],'Location','southeast')
                set(gca, 'FontSize', 20, 'LineWidth', 2, 'GridLineWidth', 1, 'FontName', FontName);
                grid on;
                box on;
            case 3
                scatter(sd_peak(b),superiority(b),200,'filled','o','MarkerFaceColor',pastel_cmap(b,:),'MarkerFaceAlpha',.8);hold on;
                xlabel({'Serial Dependence' 'Bias Peak'})
                ylabel({'Superiority' '(ortho-iso)'})
                ll1 = plot(sd_peak, polyval(polyfit(sd_peak, superiority, 1), sd_peak),'Color',[.5 .5 .5], 'LineWidth', 2);
                [rr,pp] = corr(sd_peak,superiority); legend(ll1,['r=' num2str(round(rr,2)) repmat('*', 1, (pp < 0.001) * 3 + (pp >= 0.001 & pp < 0.01) * 2 + (pp >= 0.01 & pp < 0.05) * 1)],'Location','northeast')
                set(gca, 'FontSize', 20, 'LineWidth', 2, 'GridLineWidth', 1, 'FontName', FontName);
                grid on;
                box on;
        end
    end
    set(gca,'FontSize',18)
end

set(gcf,'PaperOrientation','landscape');
set(gcf, 'PaperUnits', 'normalized');
exportgraphics(gcf, ['figures' filesep 'SI_oriBias_scatter.pdf'], 'BackgroundColor', 'none', 'ContentType', 'vector');

%% dataset superiority vs. serial dependence vs. stimulus bias peak:
tmp = tbl(abs(tbl.delta)<=90,:);
datasets = unique(tmp.code);
colors          = linspecer(length(datasets));
sd_bias = cell(length(datasets),1);
ss_bias = cell(length(datasets),1);

sd_bias_p = [];
ss_bias_p = [];
sup_p     = [];
mfreethr = [20 30];
ssfreethr = [10:25 100:115; 65:80 155:170];
figure('Units','normalized','position',[0.9 0.3 .3 .4]);
c = 1;
for i = 1:length(datasets)
    tbl_i               = tmp(tmp.codenum==i,:); 
    observers = unique(tbl_i.obs);
    sdbias = nan(length(observers),1);
    ssbias = nan(length(observers),1);
    sup    = nan(length(observers),1);
    for o = 1:length(observers)
        tbl_io               = tbl_i(tbl_i.obs==o,:);  
        sdbias(o) =  mean(tbl_io.error_iqr(abs(tbl_io.delta)<=mfreethr(2) & abs(tbl_io.delta)>=mfreethr(1)).*...
            sign(tbl_io.delta(abs(tbl_io.delta)<=mfreethr(2) & abs(tbl_io.delta)>=mfreethr(1))),'omitnan');
        ssbias(o) = mean(tbl_io.error_iqr(ismember(tbl_io.theta,ssfreethr(1,:))),'omitnan')-...
            mean(tbl_io.error_iqr(ismember(tbl_io.theta,ssfreethr(2,:))),'omitnan');
        sup(o) = std(tbl_io.error_iqr(tbl_io.bin==3),'omitnan')-std(tbl_io.error_iqr(tbl_io.bin==1),'omitnan');
    end
    idx = ~isnan(sup) & ~isnan(ssbias) & ~isnan(sdbias);
    sd_bias_p = [sd_bias_p; (sdbias(idx)-nanmean(sdbias(idx)))./nanstd(sdbias(idx))];
    ss_bias_p = [ss_bias_p; (ssbias(idx)-nanmean(ssbias(idx)))./nanstd(ssbias(idx))];
    sup_p     = [sup_p; (sup(idx)-nanmean(sup(idx)))./nanstd(sup(idx))];
    sd_bias{i}=mean(sdbias,'omitnan');
    ss_bias{i}=mean(ssbias,'omitnan');
end

% Filter table columns to include only numeric and logical types
var_types       = varfun(@class, tmp, 'OutputFormat', 'cell');
num_log_vars    = ismember(var_types, {'double', 'logical'});
filtered_tbl    = tmp(:, num_log_vars);

% Aggregate and group data by specific columns
tbl_scatter     = grpstats(filtered_tbl, {'studynum' 'codenum' 'obsid' 'bin'}, 'std');
tbl_mean          = grpstats(filtered_tbl, {'studynum' 'codenum' 'obsid'}, 'mean');
tbl_mean.sup      = tbl_scatter.std_error_iqr(tbl_scatter.bin==3) - tbl_scatter.std_error_iqr(tbl_scatter.bin==1);
tbl_datasets_mean = grpstats(tbl_mean, {'codenum'}, 'mean');
tbl_datasets_mean.sd_bias = cellfun(@mean, sd_bias);
tbl_datasets_mean.ss_bias = cellfun(@mean, ss_bias);

% subplot(131)
% m                 = fitlm(tbl_datasets_mean, 'sd_bias ~ mean_mean_stim_bias_peak');
% x_values     = linspace(min(tbl_datasets_mean.ss_bias), max(tbl_datasets_mean.ss_bias), 100)';
% tmp          = table;
% tmp.mean_mean_stim_bias_peak = x_values;
% [ypred, yci] = predict(m, tmp);
% 
% % Plot data using scatter with different colors for each dot
% numPoints = height(tbl_datasets_mean);
% for i = 1:numPoints
%     scatter(tbl_datasets_mean.ss_bias(i), tbl_datasets_mean.sd_bias(i), ...
%             50, 'filled', 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7);
%     hold on; % Keep adding to the same plot
% end
% % Plot prediction line and confidence interval
% l = plot(x_values, ypred, '-k', 'LineWidth', 3);       
% fill([x_values; flipud(x_values)], [yci(:,1); flipud(yci(:,2))], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); % Shaded CI
% [rr,pp] = corr(tbl_datasets_mean.ss_bias,tbl_datasets_mean.sd_bias) 
% legend(l,['r=' num2str(round(rr,2)) ],'Location','northeast')
% title({'At Dataset Level'});
% xlabel('Mean Stimulus-specific Bias Peak (°)')
% ylabel({'Mean Serial Dependence Bias Peak (°)'})
% grid on; box on;
% set(gca, 'FontSize', 20, 'LineWidth', 2, 'GridLineWidth', 1 'FontName', FontName);
% 
% subplot(132)
m                 = fitlm(tbl_datasets_mean, 'mean_sup ~ ss_bias');
x_values     = linspace(min(tbl_datasets_mean.ss_bias), max(tbl_datasets_mean.ss_bias), 100)';
tmp          = table;
tmp.ss_bias  = x_values;
[ypred, yci] = predict(m, tmp);

% Plot data using scatter with different colors for each dot
numPoints = height(tbl_datasets_mean);
for i = 1:numPoints
    scatter(tbl_datasets_mean.ss_bias(i), tbl_datasets_mean.mean_sup(i), ...
            50, 'filled', 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7);
    hold on; % Keep adding to the same plot
end
% Plot prediction line and confidence interval
l = plot(x_values, ypred, '-k', 'LineWidth', 3);       
fill([x_values; flipud(x_values)], [yci(:,1); flipud(yci(:,2))], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); % Shaded CI
[rr,pp] = corr(tbl_datasets_mean.ss_bias,tbl_datasets_mean.mean_sup); 
legend(l,['r=' num2str(round(rr,2)) repmat('*', 1, (pp < 0.001) * 3 + (pp >= 0.001 & pp < 0.01) * 2 + (pp >= 0.01 & pp < 0.05) * 1)],'Location','southeast')
title({'At Dataset Level'});
xlabel('Mean Stimulus-specific Bias Peak (°)')
ylabel({'Mean Superiority (°)' '(ortho-iso)'})
grid on; box on;
set(gca, 'FontSize', 20, 'LineWidth', 2, 'GridLineWidth', 1, 'FontName', FontName);

% subplot(133)
% m                 = fitlm(tbl_datasets_mean, 'mean_sup ~ sd_bias');
% x_values     = linspace(min(tbl_datasets_mean.sd_bias), max(tbl_datasets_mean.sd_bias), 100)';
% tmp          = table;
% tmp.sd_bias = x_values;
% [ypred, yci] = predict(m, tmp);
% 
% % Plot data using scatter with different colors for each dot
% numPoints = height(tbl_datasets_mean);
% for i = 1:numPoints
%     scatter(tbl_datasets_mean.sd_bias(i), tbl_datasets_mean.mean_sup(i), ...
%             50, 'filled', 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7);
%     hold on; % Keep adding to the same plot
% end
% % Plot prediction line and confidence interval
% l = plot(x_values, ypred, '-k', 'LineWidth', 3);       
% fill([x_values; flipud(x_values)], [yci(:,1); flipud(yci(:,2))], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); % Shaded CI
% [rr,pp] = corr(tbl_datasets_mean.sd_bias,tbl_datasets_mean.mean_sup); 
% legend(l,['r=' num2str(round(rr,2)) ],'Location','southeast')
% title({'At Dataset Level'});
% xlabel('Mean Serial Dependence Bias Peak (°)')
% ylabel({'Mean Superiority (°)' '(ortho-iso)'})
% grid on; box on;
% set(gca, 'FontSize', 20, 'LineWidth', 2, 'GridLineWidth', 1 'FontName', FontName);

set(gcf,'PaperOrientation','landscape')
set(gcf, 'PaperUnits', 'normalized');
exportgraphics(gcf, ['figures' filesep 'SI_oriBias_scatter_dataset.pdf'], 'BackgroundColor', 'none', 'ContentType', 'vector');

%% Show bias cleaning effects on scatter (as done in Fritsche et al. 2020)
figure('Units','normalized','position',[0 0 1 1]);
tmp = tbl(abs(tbl.delta)<=90,:);
tmp.theta= round(tmp.theta);
tmp.delta = round(tmp.delta);
datasets = unique(tmp.codenum);
bin_size = 21;
rep_size = 3;

for i = 1:length(datasets)
    subplot(7,7,i)

    tt = tmp(tmp.codenum==datasets(i),:);
    obsids = unique(tt.obsid);
    scatt_subj_iqr = nan(91,length(obsids));
    scatt_subj_ori_deb = nan(91,length(obsids));

    for oo = 1:length(unique(tt.obsid))
        tto = tt(tt.obsid == obsids(oo),:);
        tto.error_iqr(tto.delta<0) = -tto.error_iqr(tto.delta<0);
        scatt = nan(1,91);
        scatt(ismember([0:90],unique(abs(tto.delta)))) = grpstats(tto.error_iqr,abs(tto.delta),'std');scatt(scatt==0)=NaN;
        scatt_subj_iqr(:,oo)  = movmean(scatt,30,'omitnan')-mean(movmean(scatt,30,'omitnan'));

        tto.error_ori_deb(tto.delta<0) = -tto.error_ori_deb(tto.delta<0);
        scatt = nan(1,91);
        scatt(ismember([0:90],unique(abs(tto.delta)))) = grpstats(tto.error_ori_deb,abs(tto.delta),'std');scatt(scatt==0)=NaN;
        scatt_subj_ori_deb(:,oo)  = movmean(scatt,30,'omitnan')-mean(movmean(scatt,30,'omitnan'));
    end
    [p2,f2] = plot_mean_ci(repmat(0:90,1,length(obsids)),scatt_subj_iqr(:),1,1,'r','mean',0);
    [p3,f3] = plot_mean_ci(repmat(0:90,1,length(obsids)),scatt_subj_ori_deb(:),1,1,'b','mean',0);
    set(f2,'Vertices',[f2.Vertices(:,1), f2.Vertices(:,2)- p2.YData(1)])
    set(p2, 'YData', p2.YData-p2.YData(1));
    set(f3,'Vertices',[f3.Vertices(:,1), f3.Vertices(:,2)- p3.YData(1)])
    set(p3, 'YData', p3.YData-p3.YData(1));

    title(tt.code(1))
    xlabel('|Δ| (°)')
    ylabel('Error Scatter (°)')
    xticks(0:45:90)
    plot([0 90], [0 0],'k--','LineWidth',.5)
    grid on;
    set(gca,'FontSize',11)
end
set(gcf,'PaperOrientation','landscape')
set(gcf, 'PaperUnits', 'normalized');
exportgraphics(gcf, ['figures' filesep 'SI_allDatasetsScatter.pdf'], 'BackgroundColor', 'none', 'ContentType', 'vector');

%% Initialize variables for polynomial fitting
nanunique = @(x) unique(x(~isnan(x)));
bin_size = 21;
tbl_scatter = [];
all_curve = [];
pol_deg = [];
aggregated_mvav  = [];
aggregated_stimType  = [];

aggregated_delta = [];
aggregated_scatt = [];
aggregated_counts = [];
obsnum = [];
isOrientationDataset = find(grpstats(tbl(:,{'codenum','stimtype'}),{'codenum'},'mean').mean_stimtype);
oo = 1;
error_variable = 'error_ori_deb'; % error_iqr   error_ori_deb   error_iqr_norm   error_ori_deb_norm
% Loop over studies
for i = 1:length(unique(tbl.studynum))
    % if i ==5 || i ==12
    %     continue
    % end
    tbl_i = tbl(tbl.studynum==i,:);
    tbl_i = tbl_i(abs(tbl_i.delta) <= 90, :);  % Restrict to |delta| <= 90
    tbl_i.theta = round(tbl_i.theta);         % Round angles
    tbl_i.delta = round(tbl_i.delta);

    % fold if sd_deb is not removed already...
    if ~contains(error_variable,'sd_deb')
        err = tbl_i.(error_variable);
        err(tbl_i.delta<0 & tbl_i.delta~=-90) = -err(tbl_i.delta<0 & tbl_i.delta~=-90);
        tbl_i.(error_variable) = err;
    end
    
    % Loop over experiments and conditions
    nexps = max(tbl_i.expnum);
    for k = 1:nexps
        tbl_i_k = tbl_i(tbl_i.expnum==k,:);
        cond = unique(tbl_i_k.cond);
        ncond = numel(cond);

        for j = 1:ncond
            tbl_i_k_j = tbl_i_k(tbl_i_k.cond==cond(j),:);
            obs = unique(tbl_i_k_j.obs);
            nobs = numel(obs);

            for o = 1:nobs
                tbl_i_k_j_o = tbl_i_k_j(tbl_i_k_j.obs==obs(o),:);

                % Polynomial fitting
                scatt = grpstats(tbl_i_k_j_o.(error_variable), abs(tbl_i_k_j_o.delta), 'std');
                [counts, grpp] = groupcounts(abs(tbl_i_k_j_o.delta));
                delta_fit = nanunique(abs(tbl_i_k_j_o.delta));
                idx   = scatt ~= 0 & ~isnan(scatt);
                scatt_fit = scatt(idx);
                delta_fit = delta_fit(idx);
                counts = counts(idx);
                weights = counts./max(counts) ;
                
                % Determine polynomial degree range
                if length(counts) <= 10
                    degrees = [2:min(length(counts) - 1, 3)];
                else
                    degrees = [2:min(length(counts) - 1, 3)];
                end

                % BIC calculation
                BIC_values = [];
                for degree = degrees
                    fitOptions = fitoptions('Method', 'LinearLeastSquares', 'Weights', weights);
                    model = fit(delta_fit, scatt_fit, ['poly' num2str(degree)], fitOptions);

                    % Calculate residuals and BIC
                    y_fit = feval(model, delta_fit);
                    residuals = scatt_fit - y_fit;
                    SSR = sum((sqrt(weights) .* residuals).^2./sum(weights));
                    numParams = degree + 1;  % Polynomial degree + 1 for intercept
                    n = length(delta_fit);
                    BIC = n * log(SSR / n) + numParams * log(n);
                    BIC_values = [BIC_values; BIC];
                end

                % Select best polynomial degree using BIC
                [sortedBIC, sortedIndex] = sort(BIC_values);
                bestDegree = degrees(min(sortedIndex(sortedBIC < sortedBIC(1) + 2)));

                best_model = fit(delta_fit, scatt_fit, ['poly' num2str(bestDegree)], fitOptions);
                % Predict over range 0:90
                best_fit = feval(best_model, [0:90]');

                % % Skip if scatter becomes negative or NaN
                if any(best_fit([1, 46, 91]) <= 0 | any(isnan(best_fit([1, 46, 91]))))
                    best_fit([1, 46, 91])
                    continue
                end

                % Collect results
                pol_deg = [pol_deg; bestDegree];
                % p1s= plot_mean_ci(tbl_i_k_j_o.delta,tbl_i_k_j_o.(error_variable),rep_size,bin_size,'r','std',1);hold off;
                mvav = nan(91,1); mvav(delta_fit+1) = scatt_fit; mvav = [mvav(end-1:-1:2);mvav]; mvav=repmat(mvav,3,1);
                weights = nan(91,1); weights(delta_fit+1) = counts; weights = [weights(end-1:-1:2);weights]; weights=repmat(weights,3,1);
                weighted_moving_avg = movsum(mvav.*weights,bin_size,'omitmissing')./movsum(weights,bin_size,'omitmissing');
                weighted_moving_avg = weighted_moving_avg(180+[90:180]);

       
                
                aggregated_mvav = [aggregated_mvav weighted_moving_avg]; % p1s.YData(91:end)
                aggregated_stimType = [aggregated_stimType; strcmp(tbl_i_k_j_o.stimulus(1),'Orientation')];

                all_curve = [all_curve, best_fit];
                aggregated_delta = [aggregated_delta; delta_fit];
                aggregated_scatt = [aggregated_scatt; scatt_fit];
                aggregated_counts = [aggregated_counts; counts];
                obsnum = [obsnum; oo * ones(size(scatt_fit))];
                oo = oo + 1;

                % best_fit = normalize(best_fit);
                % Store polynomial fit results
                tmp = repmat(tbl_i_k_j_o(1, ismember(tbl_i_k_j_o.Properties.VariableNames, {'obsid', 'codenum'})), 3, 1);
                tmp.SI = categorical({'iso' 'mid' 'ortho'})';
                tt = grpstats(tbl_i_k_j_o(:, {'codenum', 'stimtype', 'obsid', 'bin', error_variable}), ...
                    {'codenum', 'stimtype', 'obsid', 'bin'}, 'std');
                tmp.ES = best_fit([1, 46, 91]);
                tmp.bin_scatter = tt.(['std_' error_variable]);
                tmp.ntrials     = repmat(length(counts),3,1);
                tmp.isoutlierDataset  = repmat(ismember(tbl_i_k_j_o.codenum(1),isOrientationDataset),3,1);
                tbl_scatter = [tbl_scatter; tmp];
            end
        end
    end
end
save('tbl_scatter.mat','tbl_scatter');

%% Plot results
bin_size = 5;
x = 0:90;
figure('Units','normalized','position',[.05 .4 .9 .4]);
for i = 1:3
    subplot(1,3,i)
    switch i
        case 1
            % All datasets
            out = zeros(size(all_curve,2),1);
            out2 = zeros(size(aggregated_trials,1),1);
            whichDatasets = 'All';
        case 2
            % Identify orientation datasets
            out = ~tbl_scatter.isoutlierDataset(1:3:end);
            out2 = aggregated_stimType == 0;
            whichDatasets = 'Orientation';
        case 3
            % Identify motion datasets
            out = tbl_scatter.isoutlierDataset(1:3:end);
            out2 = aggregated_stimType == 1;
            whichDatasets = 'Motion';
    end

    % Filter aggregated data based on outliers
    weighted_agg_scatt = aggregated_scatt(ismember(obsnum, find(~out)));
    weighted_agg_delta = aggregated_delta(ismember(obsnum, find(~out)));
    weighted_agg_counts = aggregated_counts(ismember(obsnum, find(~out)));
    aggregated_fit_curve = all_curve(:, ~out); aggregated_fit_curve = aggregated_fit_curve(:)';
    aggregated_mvav_curve = aggregated_mvav(:,~out2); aggregated_mvav_curve = aggregated_mvav_curve(:)';

    % Calculate average scatter for each delta
    all_scatter = grpstats(weighted_agg_scatt, weighted_agg_delta, 'mean');
    % Scatter plot of aggregated scatter data
    hold on;
    % scatter(unique(weighted_agg_delta), all_scatter,15*grpstats(weighted_agg_counts, weighted_agg_delta, 'mean'), 'filled', 'MarkerFaceColor', '#2B9AA9', 'MarkerFaceAlpha', .7);
    weights = grpstats(weighted_agg_counts, weighted_agg_delta, 'sum')./sum(grpstats(weighted_agg_counts, weighted_agg_delta, 'mean'));
    weights = repmat([weights(end-1:-1:2);weights],3,1);
    all_scatter = repmat([all_scatter(end-1:-1:2);all_scatter],3,1);

    p=plot_mean_ci([repmat(0:90,1,sum(~out2))],aggregated_mvav_curve,1,1,'b','mean',[]);
    p=plot_mean_ci([repmat(0:90,1,sum(~out))],aggregated_fit_curve,1,1,'k','mean',[]);
    p.LineWidth = 3;

    % Configure plot
    xlabel('|Δ| (°)');
    ylabel('Error Scatter (°)');
    title([whichDatasets ' Datasets']); % ' (' num2str(sum(~out)) '/' num2str(oo-1) ') code:obsid'
    xticks([0 45 90]);
    set(gca, 'FontSize', 20, 'LineWidth', 2, 'GridLineWidth', 1, 'FontName', FontName);

    % xticklabels({'iso' 'mid' 'ortho'});
end
set(gcf,'PaperOrientation','landscape')
set(gcf, 'PaperUnits', 'normalized');
exportgraphics(gcf,['figures' filesep 'main_results_polyfit.pdf'],'BackgroundColor','none')


%% Figure for bin comparison with boxplots (all datasets) 
addpath('model')
Sd          = 9;    % sigma (witdh) of decoding transition distribution, fixed
Pt          = 0.5;  % probability of integration
tau         = 1;    % temporal decay parameter, fixed
nback       = 1;
meas        = 0;
theta       = datasample(0:179,10000)';
opt_scat   = @(delta,sigma) sqrt(1./(2+(delta/sigma).^2).^2.*sigma.^2+(1-1./(2+(delta/sigma).^2)).^2.*sigma^2);
bincode = {'iso' 'mid' 'ortho'};
for i = 1:max(tbl_scatter.codenum)
    tmp = tbl_scatter(tbl_scatter.codenum==i,:);
    simulationScat = mean(tmp.bin_scatter);
    model = SD_ma_model_bayesian(simulationScat,Sd,Pt,tau,theta,nback,meas);
    optIntSig = opt_scat(model.delta,simulationScat);
    for k = 1:3
        bayes(i,k) = nanmean(model.sigma(ismember(abs(model.delta),45*(k-1))));
        optInt(i,k) = nanmean(optIntSig(ismember(abs(model.delta),45*(k-1))));
        study_scatter(i,k) = nanmean(tmp.bin_scatter(tmp.SI == bincode{k}));
    end
end
%%
figure('Units','normalized','position',[.05 .4 .9 .4]);
subplot(131)
plot(100*(bayes./bayes(:,1)-1)','Color',[0, 0, 0, 0.1])
hold on;
errorbar(1:3,100*nanmean(bayes./bayes(:,1)-1),nanstd(100*(bayes./bayes(:,1)-1)),'LineWidth',2,'Color','#023858')
plot([.5 3.5],[0 0],'k--','LineWidth',1)
title('Bayesian Integration')
ylabel('Percent change w.r.t. iso bin')
ylim([-20 60]);xlim([0.5 3.5]);
xticks(1:3);xticklabels(bincode)
set(gca, 'FontSize', 20, 'LineWidth', 2, 'GridLineWidth', 1, 'FontName', FontName);

subplot(132)
plot(100*(optInt./optInt(:,1)-1)','Color',[0, 0, 0, 0.1])
hold on;
errorbar(1:3,100*nanmean(optInt./optInt(:,1)-1),nanstd(100*(optInt./optInt(:,1)-1)),'LineWidth',2,'Color','#74A9CF')
plot([.5 3.5],[0 0],'k--','LineWidth',1)
title('Optimal Cue Integration')
ylabel('Percent change w.r.t. iso bin')
ylim([-20 60]);xlim([0.5 3.5]);
xticks(1:3);xticklabels(bincode)
set(gca, 'FontSize', 20, 'LineWidth', 2, 'GridLineWidth', 1, 'FontName', FontName);


subplot(133)
plot(100*(study_scatter./study_scatter(:,1)-1)','Color',[0, 0, 0, 0.1])
hold on;
errorbar(1:3,100*nanmean(study_scatter./study_scatter(:,1)-1),nanstd(100*(study_scatter./study_scatter(:,1)-1)),'LineWidth',2,'Color','k')
plot([.5 3.5],[0 0],'k--','LineWidth',1)
title('Real data')
ylabel('Percent change w.r.t. iso bin')
ylim([-20 60]);xlim([0.5 3.5]);
xticks(1:3);xticklabels(bincode)
set(gca, 'FontSize', 20, 'LineWidth', 2, 'GridLineWidth', 1, 'FontName', FontName);
set(gcf,'PaperOrientation','landscape')
set(gcf, 'PaperUnits', 'normalized');
exportgraphics(gcf, ['figures' filesep 'SI_comparisons.pdf'], 'BackgroundColor', 'none', 'ContentType', 'vector');
%%
datasets        = unique(tbl.code);
row_spacing     = 18;
x_spacing       = 1:row_spacing:row_spacing*max(tbl.codenum);
colors          = linspecer(max(tbl.codenum));%
% main figure options
optfg           = [];
optfg.dotshift  = 8;
optfg.dotssize  = 8;
optfg.jitfactor = .5;
optfg.showdots  = true;
optfg.boxshift  = 0.05;
optfg.rescaling = 18;
optfg.connline  = false;
optfg.facecolor = 'sameas_noalpha';
optfg.viorefine = false;
optfg.ksdensamp = 10;

figure('Units','normalized','position',[0 0 .7 1]);
tiledlayout(1,3,'TileSpacing','none','Padding','tight')
compnames       = {'\bf Iso vs Mid','\bf Iso vs Ortho','\bf Ortho vs Mid'};
compindex       = [1 2; 1 3; 3 2];
medcomp_all     = cell(3,1);
n_datasets      = max(tbl_scatter.codenum);
% scatter_bin     = {iso,mid,ort};

for k = 1:3
    nexttile
    % collect all medians
    medcomp         = nan(max(tbl_scatter.codenum),1);
    effs            = medcomp;
    for i = 1:n_datasets
        tmp = tbl_scatter(tbl_scatter.codenum==i,:);
        scatter_bin = {tmp.ES(tmp.SI == 'iso')  tmp.ES(tmp.SI == 'mid')  tmp.ES(tmp.SI == 'ortho')};
        comparison  = scatter_bin{compindex(k,1)} - scatter_bin{compindex(k,2)};
        compBayes(i)  = bayes(i,compindex(k,1))-bayes(i,compindex(k,2));
        compOptInt(i)  = optInt(i,compindex(k,1))-optInt(i,compindex(k,2));
        hedgesg     = computeHedges_g(scatter_bin{compindex(k,1)},scatter_bin{compindex(k,2)});
        optfg.x     = x_spacing(i);
        optfg.color = colors(i,:);
        plot_violin(comparison,optfg);
        hold on
        medcomp(i)  = median(comparison,'omitnan');
        effs(i)     = hedgesg;
        % add the effect size (Hedge's g)
        if abs(hedgesg)>=.2
            scatter(optfg.x,median(comparison),abs(hedgesg)*50,'markerfacecolor',optfg.color,'markeredgecolor',[.2 .2 .2],'LineWidth',1);
            % small effect size criterion as reference
            scatter(optfg.x,median(comparison),.2*50,'markerfacecolor',[1 1 1],'markeredgecolor',[1 1 1],'LineWidth',1);
        else
            % smaller than the effect size criterion of reference
            scatter(optfg.x,median(comparison),.2*50,'markerfacecolor',[1 1 1],'markeredgecolor',[1 1 1],'LineWidth',1);
            % scatter(optfg.x,median(comparison),abs(ef.hedgesg)*50,'markerfacecolor',optfg.color,'markeredgecolor',[.2 .2 .2]);
        end
        set(gca, 'FontSize', 14, 'LineWidth', 2, 'GridLineWidth', 1, 'FontName', FontName,'TickLength', [0 0]); % Removes both major ticks

    end

    optfg.x         = x_spacing(i)+40;
    optfg.color     = [0.01 0.44 0.65];
    optfg.ksdensamp = 15;
    % then add the median statistic
    plot_violin(compBayes',optfg);
    hold on
    % add the median effect size (Hedge's g)
    scatter(optfg.x,median(compBayes),abs(computeHedges_g(compBayes,zeros(size(compBayes)))*50),'markerfacecolor',optfg.color,'markeredgecolor',[.2 .2 .2],'LineWidth',1);
    % small effect size criterion as reference
    scatter(optfg.x,median(compBayes),.2*50,'markerfacecolor',[1 1 1],'markeredgecolor',[1 1 1]);

    optfg.x         = x_spacing(i)+60;
    optfg.color     = '#74A9CF';
    optfg.ksdensamp = 15;
    % then add the median statistic
    plot_violin(compOptInt',optfg);
    hold on
    % add the median effect size (Hedge's g)
    scatter(optfg.x,median(compOptInt),abs(computeHedges_g(compOptInt,zeros(size(compOptInt)))*50),'markerfacecolor',optfg.color,'markeredgecolor',[.2 .2 .2],'LineWidth',1);
    % small effect size criterion as reference
    scatter(optfg.x,median(compOptInt),.2*50,'markerfacecolor',[1 1 1],'markeredgecolor',[1 1 1]);
    
    optfg.x         = x_spacing(i)+80;
    optfg.color     = [0.5 0.5 0.5];
    optfg.ksdensamp = 15;
    % then add the median statistic
    plot_violin(medcomp,optfg)
    hold on
    % add the median effect size (Hedge's g)
    scatter(optfg.x,median(medcomp),abs(median(effs)*50),'markerfacecolor',optfg.color,'markeredgecolor',[.2 .2 .2],'LineWidth',1);
    % small effect size criterion as reference
    scatter(optfg.x,median(medcomp),.2*50,'markerfacecolor',[1 1 1],'markeredgecolor',[1 1 1]);
    medcomp_all{k} = medcomp;
    

    % aestetichs
    box on;
    % grid on;
    view([90 -90])
    set(gca,'ytick',-8:4:8)
    ylim([-15 15])
    if k==1
        datname                   = strrep(datasets,'_',' ');
        datname{n_datasets+1}     = '\bf Bayesian Model Prediction';
        datname{n_datasets+2}     = '\bf Optimal Integration Prediction';
        datname{n_datasets+3}     = '\bf Median (Data)';
        set(gca,'xtick',[x_spacing x_spacing(end)+(40:20:80)],'xticklabel',datname,'XDir','reverse','XGrid','on')
    else
        set(gca,'xtick',[x_spacing x_spacing(end)+(40:20:80)],'xticklabel','','XDir','reverse','XGrid','on')
    end
    hline(0,'k--');
    xlim([-20 x_spacing(end)+95]);
    % box on
    tl                        = title(compnames{k});     tl.FontSize = 16; 
end
set(gcf,'PaperOrientation','landscape')
set(gcf, 'PaperUnits', 'normalized');
exportgraphics(gcf, ['figures' filesep 'main_comparisons.pdf'], 'BackgroundColor', 'none', 'ContentType', 'vector');


%% 
figure('Units', 'normalized', 'Position', [0 0 0.8 0.75]);
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'tight');

titles = {'\bf Bias', '\bf Scatter'};
n_datasets = max(tbl_scatter.codenum);
medscat = nan(n_datasets, 1);
datasets = unique(tbl.code);
mfreethr = [0 45];

for k = 1:2
    medvalues = nan(n_datasets, 1);
    nexttile;
    
    for i = 1:n_datasets
        tmp = tbl(tbl.codenum == i, :);
        optfg.x = x_spacing(i);
        optfg.color = colors(i, :);
        
        if k == 1 % Bias
            obs = unique(tmp.obsid);
            values = nan(size(obs));
            
            for o = 1:length(obs)
                tbl_io = tmp(tmp.obsid == obs(o), :);
                valid_delta = abs(tbl_io.delta) <= mfreethr(2) & abs(tbl_io.delta) >= mfreethr(1);
                values(o) = mean(tbl_io.error_ori_deb(valid_delta) .* sign(tbl_io.delta(valid_delta)), 'omitnan');
            end

            hg = computeHedges_g(values, zeros(size(values)));
            plot_violin(values, optfg);
            hold on;

            % Add the effect size (Hedge's g)
            scatter(optfg.x, median(values), abs(hg * 50), 'MarkerFaceColor', optfg.color, 'MarkerEdgeColor', [0.2 0.2 0.2]);
            scatter(optfg.x, median(values), 0.2 * 50, 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [1 1 1]);
            medvalues(i) = median(values);
        else % Scatter
            values = grpstats(tmp.error_ori_deb_sd_deb, tmp.obsid, 'std');
            plot_violin(values, optfg);
            hold on;

            % Add small effect size criterion as reference
            scatter(optfg.x, median(values), 0.2 * 50, 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [1 1 1]);
            medvalues(i) = median(values);
        end
    end
    disp(['Hedge''s g on the population level: ' num2str(computeHedges_g(medvalues,zeros(size(medvalues))))])
    % Plot median statistics
    optfg.x = x_spacing(i) + 60;
    optfg.color = [0.5 0.5 0.5];
    optfg.ksdensamp = 15;
    plot_violin(medvalues, optfg);
    hold on;

    if k == 1
        % Add median effect size (Hedge's g)
        hg_med = median(computeHedges_g(medvalues, zeros(size(medvalues))));
        scatter(optfg.x, median(medvalues), abs(hg_med * 50), 'MarkerFaceColor', optfg.color, 'MarkerEdgeColor', [0.2 0.2 0.2], 'LineWidth', 1);
    end

    scatter(optfg.x, median(medvalues), 0.2 * 50, 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [1 1 1]);

    % Aesthetics
    box on;
    view([90 -90]);
    set(gca, 'FontSize', 14, 'LineWidth', 2, 'FontName', FontName,'GridLineWidth', 1, 'TickLength', [0 0]);

    if k == 1
        datname = strrep(datasets, '_', ' ');
        datname{n_datasets + 1} = '\bf Median';
        set(gca, 'XTick', [x_spacing x_spacing(end) + 60], 'XTickLabel', datname, 'XDir', 'reverse', 'XGrid', 'on');
        ylim([-7 7]);
        yticks(-6:3:6);
    else
        set(gca, 'XTick', [x_spacing x_spacing(end) + 60], 'XTickLabel', '', 'XDir', 'reverse', 'XGrid', 'on');
        ylim([0 62]);
        yticks(0:10:60);
    end

    hline(0, 'k--');
    xlim([-20 x_spacing(end) + 85]);
    tl = title(titles{k});
    tl.FontSize = 22;
end

set(gcf, 'PaperOrientation', 'landscape');
set(gcf, 'PaperUnits', 'normalized');
exportgraphics(gcf, ['figures' filesep 'main_bias_scatter.pdf'], 'BackgroundColor', 'none', 'ContentType', 'vector');


%% Plot Results
%
% figure
% hist(pol_deg,50)
% set(gca,'FontSize',15)
% xlabel('Best Degree (Lowest BIC)')
% ylabel('Count')
% set(gcf,'PaperOrientation','landscape')
% set(gcf, 'PaperUnits', 'normalized');
% saveas(gcf,'degree_histogram.eps')
% saveas(gcf,'degree_histogram.png')


%% Visualization of results
figure('Units','normalized','position',[.05 .4 .9 .4]);
titles = {'iso', 'mid', 'ortho'};
for i = 1:3
    subplot(1, 3, i);
    poly = tbl_scatter.ES(i:3:end);
    binned = tbl_scatter.bin_scatter(i:3:end);
    idx = ~isnan(poly) & ~isnan(binned);
    [rr, pp] = corr(binned(idx), poly(idx));
    scatter(binned(idx), poly(idx), 'filled','ko');
    hold on;
    l1 = lsline;
    l1.Color = 'r';
    l1.LineWidth = 3;
    xlabel('Binned Scatter Estimate');
    ylabel('Polynomial Scatter Estimate');
    title([titles{i}]);
    legend(l1, ['r = ', num2str(round(rr, 2)),],'Location','northwest')
    set(gca, 'FontSize', 24);
    set(gca,'LineWidth',2)
    set(gca,'FontName',FontName)
    grid on;
    box on;
end
set(gcf,'PaperOrientation','landscape')
set(gcf, 'PaperUnits', 'normalized');
exportgraphics(gcf,['figures' filesep 'SI_polyBinCorrelations.pdf'],'BackgroundColor','none')
%% LMM s
delete(['tables' filesep 'LME.txt'])
diary(['tables' filesep 'LME.txt'])
tbl_scatter.SI = categorical(tbl_scatter.SI);
lmm             = fitlme(tbl_scatter,'ES ~ SI + (1+SI|obsid) + (1+SI|codenum) + (1+SI|codenum:obsid)');
disp(lmm)
tbl_scatter.new_SI = repmat(categorical([1 2 1]'),height(tbl_scatter)/3,1);
reduced_lmm = fitlme(tbl_scatter, 'ES ~ new_SI + (1+new_SI|obsid) + (1+new_SI|codenum) + (1+new_SI|codenum:obsid)');

bic_full = lmm.ModelCriterion.BIC;
bic_reduced = reduced_lmm.ModelCriterion.BIC;
dBIC = bic_reduced-bic_full;
BF_reduced_over_full = 1/exp((bic_reduced-bic_full)/2);
disp(['dBIC = ' num2str(dBIC) ', and BF_reduced_over_full = ' num2str(BF_reduced_over_full)])
diary off;

lmm             = fitlme(tbl_scatter,'bin_scatter ~ SI + (1+SI|obsid) + (1+SI|codenum) + (1+SI|codenum:obsid)');
disp(lmm)
tbl_scatter.new_SI = repmat(categorical([1 2 1]'),height(tbl_scatter)/3,1);
reduced_lmm = fitlme(tbl_scatter, 'bin_scatter ~ new_SI + (1+new_SI|obsid) + (1+new_SI|codenum) + (1+new_SI|codenum:obsid)');

bic_full = lmm.ModelCriterion.BIC;
bic_reduced = reduced_lmm.ModelCriterion.BIC;
dBIC = bic_reduced-bic_full;
BF_reduced_over_full = 1/exp((bic_reduced-bic_full)/2);

disp(['dBIC = ' num2str(dBIC) ', and BF_reduced_over_full = ' num2str(BF_reduced_over_full)])


