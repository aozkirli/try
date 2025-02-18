clc  % Clear command window

% Get the full path of the current script
scriptPath = mfilename('fullpath');

% Extract directory path of the current script
[currentDir, ~, ~] = fileparts(scriptPath);

% Navigate to the current directory
cd(currentDir)

% Add in-house toolbox (must be replaced with the shared folder path containing necessary functions)
addpath(genpath('C:\Users\chare\Google Drive\Work\09_Code\BEIM_toolbox\matlabv\'))
addpath(genpath('/Users/ayberkozkirli/Documents/GitHub/BEIM_toolbox'))

% Change to 'datasets' directory
cd(['..' filesep 'data' filesep 'datasets'])
load('SD_ma_master_table.mat')

total_trials            = sum(nnz(~isnan(tbl.error_ori_deb_sd_deb)));
total_datasets          = numel(unique(tbl.codenum));

fprintf('\nTotal number of trials: %.d\n',total_trials)
fprintf('Total number of subjects: %.d\n',max(tbl.obsid))

%% Estimate the Scatter per bin
% Filter table columns to include only numeric and logical types
var_types       = varfun(@class, tbl, 'OutputFormat', 'cell');
num_log_vars    = ismember(var_types, {'double', 'logical'});
filtered_tbl    = tbl(:, num_log_vars);

% Aggregate and group data by specific columns
tbl_scatter     = grpstats(filtered_tbl, {'bin', 'studynum', 'obsid', 'codenum'}, 'std');
tbl_scatter_count = grpstats(filtered_tbl, {'bin', 'studynum', 'obsid', 'codenum'}, 'numel');

% Filter out rows with insufficient observations
low_count       = tbl_scatter_count.numel_error_ori_deb_sd_deb < 15;
tbl_scatter.numel_error_ori_deb_sd_deb(low_count, :) = nan;

%% Estimate Bias and Scatter (overall) across datasets
bias                    = cell(total_datasets,1);
typeeffect              = 'hedgesg';
effect                  = nan(total_datasets,1);
scatters                = cell(total_datasets,1);
datasets                = unique(tbl.code);
[iso,mid,ort]           = deal(bias,bias,bias);
for i = 1:total_datasets
    tbl_i               = tbl_subset(tbl,'codenum',i); 
    % collect bias
    bias{i}             = sdp_model_free_bias(tbl_i,'variables',{'delta','error_ori_deb'});
    % store effect size
    ef                  = mes(bias{i},zeros(size(bias{i})),'hedgesg','isDep',1);
    effect(i)           = ef.(typeeffect);
    % collect scatter
    tbl_i.dummy         = ones(size(tbl_i.obs));
    stat                = groupstat(tbl_i,'x','dummy','y','error_ori_deb_sd_deb','g','obs','stat','std');
    scatters{i}         = stat.std;
    % collect scatter per bin
    iso{i}              = tbl_scatter.std_error_ori_deb_sd_deb(tbl_scatter.codenum==i & tbl_scatter.bin==1);
    mid{i}              = tbl_scatter.std_error_ori_deb_sd_deb(tbl_scatter.codenum==i & tbl_scatter.bin==2);
    ort{i}              = tbl_scatter.std_error_ori_deb_sd_deb(tbl_scatter.codenum==i & tbl_scatter.bin==3);
end

%% Summary plot of Bias and Scatter
row_spacing     = 18;
x_spacing       = 1:row_spacing:row_spacing*total_datasets;
colors          = linspecer(total_datasets);%
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

%% Bias 
fig1            = figure('position',[240 140 850 700]);
tiledlayout(1,2,'TileSpacing','none','Padding','tight')
nexttile

% collect all medians
medbias         = nan(total_datasets,1);
for i = 1:total_datasets
    optfg.x     = x_spacing(i);
    optfg.color = colors(i,:);
    plot_violin(bias{i},optfg);
    hold on
    medbias(i)  = median(bias{i});
    % add the effect size (Hedge's g)
    scatter(optfg.x,median(bias{i}),abs(effect(i)*50),'markerfacecolor',optfg.color,'markeredgecolor',[.2 .2 .2]);
    % small effect size criterion as reference
    scatter(optfg.x,median(bias{i}),.2*50,'markerfacecolor',[1 1 1],'markeredgecolor',[1 1 1]);
end

optfg.x         = x_spacing(i)+40;
optfg.color     = [.5 .5 .5];
optfg.ksdensamp = 15;
% then add the median statistic
plot_violin(medbias,optfg)
hold on
% add the median effect size (Hedge's g)
scatter(optfg.x,median(medbias),median(effect)*50,'markerfacecolor',optfg.color,'markeredgecolor',[.2 .2 .2]);
% small effect size criterion as reference
scatter(optfg.x,median(medbias),.2*50,'markerfacecolor',[1 1 1],'markeredgecolor',[1 1 1]);

% aestetichs
box on;
view([90 -90])
set(gca,'ytick',-10:5:10)
datname                 = strrep(datasets,'_',' ');
datname{total_datasets+1}     = 'Average';
set(gca,'xtick',[x_spacing x_spacing(end)+40],'xticklabel',datname,'XDir','reverse','XGrid','on')
hline(0,'k--');
xlim([-20 x_spacing(end)+65]);
ylim([-18 18])
format_figure([],[],[],[],[],[],12,[],[]);
tl                      = title('Bias'); tl.FontWeight = 'normal';

%% Scatter 
medscat         = nan(total_datasets,1);
nexttile
for i = 1:total_datasets
    optfg.x     = x_spacing(i);
    optfg.color = colors(i,:);
    plot_violin(scatters{i},optfg);
    hold on
    medscat(i)  = median(scatters{i});
    % add the median
    scatter(optfg.x,median(scatters{i}),30,'markerfacecolor',[1 1 1],'markeredgecolor',[.2 .2 .2]);
end

optfg.x         = x_spacing(i)+40;
optfg.color     = [.5 .5 .5];
optfg.ksdensamp = 15;
% then add the median statistic
plot_violin(medscat,optfg)
hold on
scatter(optfg.x,median(medscat),30,'markerfacecolor',[1 1 1],'markeredgecolor',[.2 .2 .2]);

% aestetichs
box on;
set(gca,'ytick',0:10:60)
ylim([0 60])
view([90 -90])
set(gca,'xtick',[x_spacing x_spacing(end)+40],'xticklabel','','XDir','reverse','XGrid','on')
xlim([-20 x_spacing(end)+65]);
format_figure([],[],[],[],[],[],12,[],[]);
tl                      = title('Scatter'); tl.FontWeight = 'normal';
fig1.PaperPositionMode = 'manual';

%% Figure for bin comparison (all datasets)

fig2            = figure('position',[240 140 1250 700]);
tiledlayout(1,3,'TileSpacing','none','Padding','tight')
compnames       = {'Iso vs Mid','Iso vs Ortho','Ortho vs Mid'};
compindex       = [1 2; 1 3; 3 2];
medcomp_all     = cell(3,1);
scatter_bin     = {iso,mid,ort};

for k = 1:3
    nexttile

    % collect all medians
    medcomp         = nan(total_datasets,1);
    effs            = medcomp;
    for i = 1:total_datasets
        comparison  = scatter_bin{compindex(k,1)}{i} - scatter_bin{compindex(k,2)}{i};
        ef          = mes(comparison,zeros(size(comparison)),'hedgesg','isDep',1);
        optfg.x     = x_spacing(i);
        optfg.color = colors(i,:);
        plot_violin(comparison,optfg);
        hold on
        medcomp(i)  = median(comparison,'omitnan');
        effs(i)     = ef.hedgesg;
        % add the effect size (Hedge's g)
        if abs(ef.hedgesg)>=.2
        scatter(optfg.x,median(comparison),abs(ef.hedgesg)*50,'markerfacecolor',optfg.color,'markeredgecolor',[.2 .2 .2]);
        % small effect size criterion as reference
        scatter(optfg.x,median(comparison),.2*50,'markerfacecolor',[1 1 1],'markeredgecolor',[1 1 1]);
        else
        % smaller than the effect size criterion of reference
        scatter(optfg.x,median(comparison),.2*50,'markerfacecolor',[1 1 1],'markeredgecolor',[1 1 1]);
        % scatter(optfg.x,median(comparison),abs(ef.hedgesg)*50,'markerfacecolor',optfg.color,'markeredgecolor',[.2 .2 .2]);
        end
        set(gca, 'FontName', 'Arial');
    end

    optfg.x         = x_spacing(i)+40;
    optfg.color     = [0 0 0];
    optfg.ksdensamp = 15;
    % then add the median statistic
    plot_violin(medcomp,optfg)
    hold on
    % add the median effect size (Hedge's g)
    scatter(optfg.x,median(medcomp),abs(median(effs)*50),'markerfacecolor',optfg.color,'markeredgecolor',[.5 .5 .5]);
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
    datname{total_datasets+1}     = 'Median (Data)';
    set(gca,'xtick',[x_spacing x_spacing(end)+40],'xticklabel',datname,'XDir','reverse','XGrid','on')
    else
    set(gca,'xtick',[x_spacing x_spacing(end)+40],'xticklabel','','XDir','reverse','XGrid','on')
    end
    hline(0,'k--');
    xlim([-20 x_spacing(end)+65]);
    % box on
    format_figure([],[],[],[],[],[],12,[],[]);
    tl                        = title(compnames{k}); tl.FontWeight = 'normal';
end

%% Linear Mixed Model
formula              = 'std_error_ori_deb_sd_deb ~ bin + (1|obsid) + (1|codenum)';
% make sure that insufficient observations have been removed at the
% beginning of the code, make all necessary categorical variables
tbl_scatter.bin      = categorical(tbl_scatter.bin);
tbl_scatter.obsid    = categorical(tbl_scatter.obsid);
tbl_scatter.codenum  = categorical(tbl_scatter.codenum);
tbl_scatter.stimtype = categorical(tbl_scatter.std_stimtype);
lmm                  = fitlme(tbl_scatter,formula);
disp(lmm)
lmm.Rsquared

% Extract fixed effect estimates and 95% confidence intervals
[beta, betaNames, statsTable] = fixedEffects(lmm, 'alpha', 0.05);
% Extract the lower and upper confidence intervals from statsTable
ciLower      = statsTable.Lower;
ciUpper      = statsTable.Upper;

% Create a new figure
figure;
hold on;

beta(2:3)    = beta(2:3)+beta(1);
ciLower(2:3) = beta(1)+ciLower(2:3);
ciUpper(2:3) = beta(1)+ciUpper(2:3);

% Plot each fixed effect with error bars
numEffects = numel(beta);
for i = 1:numEffects
    % Plot the fixed effect estimate with its confidence interval
    errorbar(i, beta(i), beta(i) - ciLower(i), ciUpper(i) - beta(i), ...
             'o', 'MarkerSize', 8, 'LineWidth', 1.5, 'CapSize', 10, 'Color', 'b');
end

% Customize plot
set(gca, 'XTick', 1:numEffects, 'XTickLabel', {'Intercept(iso)','\beta(mid)','\beta(ortho)'});
ylabel('Fixed Effect Estimate');
title('Linear Mixed Model Results');
xlim([0 4])
format_figure(beta(1),nan,'','Fixed Effects');
hold on;

%% LMM resampling
numBoot         = 100;
% Storage for bootstrap estimates
bootEstimates   = zeros(numBoot, length(lmm.Coefficients.Estimate));
rng('default'); % For reproducibility

for i = 1:numBoot
    % Resample the data with replacement
    resampledIdx  = randsample(height(tbl_scatter), height(tbl_scatter), true);
    resampledData = tbl_scatter(resampledIdx, :);

    % Fit the model to the resampled data
    lmmResampled  = fitlme(resampledData, formula);

    % Store the estimated fixed effects coefficients
    bootEstimates(i, :) = lmmResampled.Coefficients.Estimate';

    % Print progress every 10%
    if mod(i, numBoot / 10) == 0
        fprintf('Progress (resampling LMM): %d%%\n', round(i / numBoot * 100));
    end
end

%% Plot
bootBeta        = bootEstimates;
bootBeta(:,2:3) = bootBeta(:,2:3)+bootBeta(:,1);

fig3            = figure('position',[240 140 540 700]);
optfg           = [];
optfg.color     = [.5 .5 .5];
optfg.showdots  = 1;
optfg.ksdensamp = .1;
optfg.rescaling = .25;
optfg.dotssize  = 10;
optfg.dotshift  = .2;
optfg.jitfactor = .05;
optfg.connline  = false;
optfg.facecolor = 'sameas';
optfg.facealpha = .7;

plot_violin(bootBeta,optfg);
for i = 1:numEffects
    % Plot the fixed effect estimate (only slopes) with its confidence interval
    errorbar(i, beta(i), beta(i) - ciLower(i), ciUpper(i) - beta(i), ...
             'o', 'MarkerSize', 8, 'LineWidth', 1.5, 'CapSize', 10, 'Color', 'k');
end
% Customize plot
set(gca, 'XTick', 1:numEffects, 'XTickLabel', {'Intercept(iso)','\beta(mid)','\beta(ortho)'});
ylabel('Fixed Effect Estimate');
title('Linear Mixed Model Results');
xlim([0 4])
format_figure(beta(1),nan,'','Fixed Effects');
hold on;

%% Across datasets, relationship between superiority and stim-specific biases
% Calculate group means
tbl_mean          = grpstats(filtered_tbl, {'studynum', 'obsid', 'codenum'}, 'mean');
tbl_mean.sup      = tbl_scatter.std_errorsd(tbl_scatter.bin=="3") - tbl_scatter.std_errorsd(tbl_scatter.bin=="1");
tbl_datasets_mean = grpstats(tbl_mean, {'codenum'}, 'mean');
% Add the average bias (SD) too
tbl_datasets_mean.sd_bias = cellfun(@mean, bias);


% Fit linear model with the peak of stimulus-specific bias as predictor
% m                 = fitlm(tbl_datasets_mean.mean_mean_stim_bias_peak, tbl_datasets_mean.mean_sup);
m                 = fitlm(tbl_datasets_mean, 'mean_sup ~ mean_mean_stim_bias_peak');

% Plot data 
figure;
subplot(121)
numPoints = height(tbl_datasets_mean);
% Plot data using scatter with different colors for each dot
for i = 1:numPoints
    scatter(tbl_datasets_mean.mean_mean_stim_bias_peak(i), tbl_datasets_mean.mean_sup(i), ...
            50, 'filled', 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7);
    hold on; % Keep adding to the same plot
end

% Extract fitted values and confidence intervals
x_values     = linspace(min(tbl_datasets_mean.mean_mean_stim_bias_peak), max(tbl_datasets_mean.mean_mean_stim_bias_peak), 100)';
tmp          = table;
tmp.mean_mean_stim_bias_peak = x_values;
[ypred, yci] = predict(m, tmp);

% Plot prediction line and confidence interval
plot(x_values, ypred, '-k', 'LineWidth', 1.5);        % Red prediction line
fill([x_values; flipud(x_values)], [yci(:,1); flipud(yci(:,2))], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); % Shaded CI

% Customize plot appearance
title({'Relationship Between Superiority', 'and Stimulus-Specific Biases Across Datasets'});
% legend('Data', 'Prediction', '95% CI');
% grid on;
box on;
format_figure(0,nan,'Stimulus-specific peak bias','Superiority')
hold off;


% Fit linear model with the SD bias as predictor
% m                 = fitlm(tbl_datasets_mean.mean_mean_stim_bias_peak, tbl_datasets_mean.mean_sup);
m                 = fitlm(tbl_datasets_mean, 'mean_sup ~ sd_bias');

% Plot data 
subplot(122)
numPoints = height(tbl_datasets_mean);
% Plot data using scatter with different colors for each dot
for i = 1:numPoints
    scatter(tbl_datasets_mean.sd_bias(i), tbl_datasets_mean.mean_sup(i), ...
            50, 'filled', 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7);
    hold on; % Keep adding to the same plot
end

% Extract fitted values and confidence intervals
x_values     = linspace(min(tbl_datasets_mean.sd_bias), max(tbl_datasets_mean.sd_bias), 100)';
tmp          = table;
tmp.sd_bias  = x_values;
[ypred, yci] = predict(m, tmp);

% Plot prediction line and confidence interval
plot(x_values, ypred, '-k', 'LineWidth', 1.5);        % Red prediction line
fill([x_values; flipud(x_values)], [yci(:,1); flipud(yci(:,2))], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); % Shaded CI

% Customize plot appearance
title({'Relationship Between Superiority', 'and SD Biases Across Datasets'});
% legend('Data', 'Prediction', '95% CI');
% grid on;
box on;
format_figure(0,nan,'Serial dependence bias','Superiority')
hold off;
