clear all;
close all;
clc
%% bias cleaning check
load(['..' filesep '..' filesep '..' filesep 'data' filesep 'datasets' filesep 'SD_ma_master_table.mat'])
tbl.codenum_obsid = strcat(string(tbl.codenum),'_',string(tbl.obsid));
% tbl = tbl(tbl.studynum~=8 & tbl.studynum~=8,:);
% Check how cleaning went:
figure('Units','normalized','position',[.05 .4 .6 .4]);
tmp = tbl(abs(tbl.delta)<=90 & tbl.theta<=180,:);
% Orientation bias removal
rep_size = 3;
bin_size = 21;

subplot(121)
p1 = plot_mean_ci_participant_level(tmp.theta,tmp.error_iqr_norm,tmp.codenum_obsid,rep_size,bin_size,'k');
hold on;
p2 = plot_mean_ci_participant_level(tmp.theta,tmp.error_ori_deb_norm,tmp.codenum_obsid,rep_size,bin_size,'b');
p3 = plot_mean_ci_participant_level(tmp.theta,tmp.error_ori_deb_sd_deb_norm,tmp.codenum_obsid,rep_size,bin_size,'r');
ylim([-.55 .55])
% legend([p1 p2 p3],{'outliers removed' 'stimulus bias removed' 'serial dependence bias removed'}, 'Location', 'best')
xlabel('theta')
ylabel('stimulus-specific bias')

subplot(122)
p1 = plot_mean_ci_participant_level(tmp.delta,tmp.error_iqr_norm,tmp.codenum_obsid,rep_size,bin_size,'k');
hold on;
p2 = plot_mean_ci_participant_level(tmp.delta,tmp.error_ori_deb_norm,tmp.codenum_obsid,rep_size,bin_size,'b');
hold on;
p3 = plot_mean_ci_participant_level(tmp.delta,tmp.error_ori_deb_sd_deb_norm,tmp.codenum_obsid,rep_size,bin_size,'r');
ylim([-.55 .55])
xlabel('delta')
ylabel('serial dependence bias')
legend([p1 p2 p3],{'outliers removed' 'stimulus bias removed' 'serial dependence bias removed'}, 'Location', 'best')

saveas(gcf,'removal_check.png')
saveas(gcf,'removal_check.eps')

%% Show ortho vs. iso ori biases
tmp    = tbl(abs(tbl.delta)<=90 & tbl.theta<=180,:);
iso    = tmp(ismember(abs(tmp.delta),0:10),:);
ortho  = tmp(ismember(abs(tmp.delta),80:90),:);
rep_size = 3;
bin_size = 35;

dd = nan(1000,181);
groups = unique(tmp.codenum_obsid);
for gg = 1:length(groups)
    tmpp = tmp(strcmp(tmp.codenum_obsid,groups(gg)),:);
    for i = 0:180
        if nanstd(tmpp.error_iqr_norm(tmpp.theta == i)) ~=0
            dd(gg,i+1) = nanmean(tmpp.error_iqr_norm(tmpp.theta == i))/nanstd(tmpp.error_iqr_norm(tmpp.theta == i));
        end
    end
end
[abs_dd_sorted,idx]=sort(abs(nanmean(dd)));
n_bins = 4;
figure('Units','normalized','position',[0 .4 1 .3]);
subplot(141)
p1 = plot_mean_ci_participant_level(iso.theta,iso.error_iqr_norm,iso.codenum_obsid,rep_size,bin_size,'k');
highlight = nan(181,1); highlight(idx(181-180/n_bins:181)) = 0;
h1=plot(0:180,highlight,'Color','#FF7F0E','LineWidth',5);
p2 = plot_mean_ci_participant_level(ortho.theta,ortho.error_iqr_norm,ortho.codenum_obsid,rep_size,bin_size,'r');
highlight = nan(181,1); highlight(idx(1:180/n_bins)) = 0;
h2=plot(0:180,highlight,'Color','#17BECF','LineWidth',5);
xlabel('theta')
ylabel('stimulus-specific bias')
legend([p1 p2 h1 h2],{'iso' 'ortho' 'strong shift' 'no shift'}, 'Location', 'best')

% check sort the data based on cohen's d for orientation bias:
card     = tmp(ismember(abs(tmp.theta),idx(181-180/n_bins:181)-1),:);
oblique  = tmp(ismember(abs(tmp.theta),idx(1:180/n_bins)-1),:);

subplot(142)
p1 = plot_mean_ci_participant_level(card.delta,card.error_iqr_norm,card.codenum_obsid,rep_size,bin_size,'#FF7F0E','mean');
hold on;
p2 = plot_mean_ci_participant_level(oblique.delta,oblique.error_iqr_norm,oblique.codenum_obsid,rep_size,bin_size,'#17BECF','mean');
xlabel('delta')
ylabel('serial dependence bias')
legend([p1 p2],{'strong shift' 'no shift'}, 'Location', 'best')

subplot(143)
p1 = plot_mean_ci_participant_level(card.delta,card.error_iqr_norm,card.codenum_obsid,rep_size,bin_size,'#FF7F0E','std');
hold on;
p2 = plot_mean_ci_participant_level(oblique.delta,oblique.error_iqr_norm,oblique.codenum_obsid,rep_size,bin_size,'#17BECF','std');
xlabel('delta')
ylabel('error scatter')
legend([p1 p2],{'strong shift' 'no shift'}, 'Location', 'best')

subplot(144)
tmp    = tbl(abs(tbl.delta)<=90,:);
var_types       = varfun(@class, tmp, 'OutputFormat', 'cell');
num_log_vars    = ismember(var_types, {'double', 'logical'});
filtered_tbl    = tmp(:, num_log_vars);
tmp_scatter     = grpstats(filtered_tbl, {'codenum' 'obsid' 'bin'}, 'std');
tbl_mean          = grpstats(filtered_tbl, {'codenum' 'obsid'}, 'mean');
tbl_mean.sup      = tmp_scatter.std_error_iqr_norm(tmp_scatter.bin==3) - tmp_scatter.std_error_iqr_norm(tmp_scatter.bin==1);
tbl_datasets_mean = grpstats(tbl_mean, {'codenum'}, 'mean');
colors          = linspecer(length(unique(tmp.codenum)));%
numPoints = height(tbl_datasets_mean);
% Plot data using scatter with different colors for each dot
for i = 1:numPoints
    scatter(tbl_datasets_mean.mean_mean_stim_bias_peak(i), tbl_datasets_mean.mean_sup(i), ...
            50, 'filled', 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7);
    hold on; % Keep adding to the same plot
end
m                 = fitlm(tbl_datasets_mean, 'mean_sup ~ mean_mean_stim_bias_peak');
% Extract fitted values and confidence intervals
x_values     = linspace(min(tbl_datasets_mean.mean_mean_stim_bias_peak), max(tbl_datasets_mean.mean_mean_stim_bias_peak), 100)';
tmp          = table;
tmp.mean_mean_stim_bias_peak = x_values;
[ypred, yci] = predict(m, tmp);

% Plot prediction line and confidence interval
plot(x_values, ypred, '-k', 'LineWidth', 1.5);        % Red prediction line
fill([x_values; flipud(x_values)], [yci(:,1); flipud(yci(:,2))], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); % Shaded CI
set(gca,'FontSize',15)

% Customize plot appearance
title('At dataset level');
xlabel('stimulus-specific peak bias')
ylabel('superiority')
hold off;
saveas(gcf,'iso-ortho-superiority-ori-bias.png')
saveas(gcf,'iso-ortho-superiority-ori-bias.eps')

%% Show bias cleaning effects on scatter
close all;
figure('Units','normalized','position',[0 0 1 1]);
rep_size = 3;
bin_size = 35;
tmp = tbl(abs(tbl.delta)<=90,:);
datasets = unique(tmp.codenum);
for i = 1:length(datasets)
    tt = tmp(tmp.codenum==datasets(i),:);
    subplot(7,7,i)
    p1 = plot_mean_ci_participant_level(tt.delta,tt.error_norm,tt.codenum_obsid,rep_size,bin_size,'k','std');
    p2 = plot_mean_ci_participant_level(tt.delta,tt.error_iqr_norm,tt.codenum_obsid,rep_size,bin_size,'b','std');
    p3 = plot_mean_ci_participant_level(tt.delta,tt.error_ori_deb_norm,tt.codenum_obsid,rep_size,bin_size,'r','std');
    title(tt.code(1))
    xlabel('delta')
    ylabel('scatter')
    if i == 1
        legend([p1 p2 p3],{'raw' '+no outliers' '+no ori bias'}, 'Location', 'best')
    end
    set(gca,'FontSize',11)
end
%% Initialize variables for polynomial fitting
nanunique = @(x) unique(x(~isnan(x)));
tbl_scatter = [];
all_curve = [];
pol_deg = [];
aggregated_delta = [];
aggregated_scatt = [];
aggregated_weights = [];
obsnum = [];
isOrientationDataset = find(grpstats(tbl(:,{'codenum','stimtype'}),{'codenum'},'mean').mean_stimtype);
oo = 1;
error_variable = 'error_iqr'; % error_iqr   error_ori_deb_sd_deb   error_iqr_norm   error_ori_deb_sd_deb_norm   
% Loop over studies
for i = 1:length(unique(tbl.studynum))
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

                % Normalize delta_fit
                delta_fit_mean = mean(delta_fit);
                delta_fit_std = std(delta_fit);
                delta_fit_normalized = (delta_fit - delta_fit_mean) / delta_fit_std;
               
                % BIC calculation
                BIC_values = [];
                for degree = degrees
                    fitOptions = fitoptions('Method', 'LinearLeastSquares', 'Weights', weights);
                    model = fit(delta_fit_normalized, scatt_fit, ['poly' num2str(degree)], fitOptions);

                    % Calculate residuals and BIC
                    y_fit = feval(model, delta_fit_normalized);
                    residuals = scatt_fit - y_fit;
                    SSR = sum((sqrt(weights) .* residuals).^2);
                    numParams = degree + 1;  % Polynomial degree + 1 for intercept
                    n = length(delta_fit);
                    BIC = n * log(SSR / n) + numParams * log(n);
                    BIC_values = [BIC_values; BIC];
                end
                
                % Select best polynomial degree using BIC
                [sortedBIC, sortedIndex] = sort(BIC_values);
                bestDegree = degrees(min(sortedIndex(sortedBIC < sortedBIC(1) + 2)));

                best_model = fit(delta_fit_normalized, scatt_fit, ['poly' num2str(bestDegree)], fitOptions);
                % Predict over range 0:90
                best_fit = feval(best_model, ([0:90]' - delta_fit_mean) / delta_fit_std);

                % % Skip if scatter becomes negative or NaN
                if any(best_fit([1, 46, 91]) <= 0 | any(isnan(best_fit([1, 46, 91]))))
                    continue
                end

                % Collect results
                pol_deg = [pol_deg; bestDegree];

                all_curve = [all_curve, best_fit];
                aggregated_delta = [aggregated_delta; delta_fit];
                aggregated_scatt = [aggregated_scatt; scatt_fit];
                aggregated_weights = [aggregated_weights; counts];
                obsnum = [obsnum; oo * ones(size(scatt_fit))];
                oo = oo + 1;

                % Store polynomial fit results
                tmp = repmat(tbl_i_k_j_o(1, ismember(tbl_i_k_j_o.Properties.VariableNames, {'obsid', 'codenum'})), 3, 1);
                tmp.bin = categorical(1:3)';
                tt = grpstats(tbl_i_k_j_o(:, {'codenum', 'stimtype', 'obsid', 'bin', error_variable}), ...
                    {'codenum', 'stimtype', 'obsid', 'bin'}, 'std');
                tmp.poly_scatter = best_fit([1, 46, 91]);
                tmp.bin_scatter = tt.(['std_' error_variable]);
                tmp.ntrials     = repmat(length(counts),3,1);
                tmp.isoutlierDataset  = repmat(ismember(tbl_i_k_j_o.codenum(1),isOrientationDataset),3,1);
                tbl_scatter = [tbl_scatter; tmp];
            end
        end
    end
end
save('tbl_scatter.mat','tbl_scatter');
%% Do paired t-test on iso-mid-ortho
% Start recording to a text file
delete('paired t-test.txt')
diary('paired t-test.txt');
% using polynomial:
tbl_scatter.bin =double(tbl_scatter.bin);
iso = tbl_scatter.poly_scatter(tbl_scatter.bin==1);
mid = tbl_scatter.poly_scatter(tbl_scatter.bin==2);
ortho = tbl_scatter.poly_scatter(tbl_scatter.bin==3);
weights = ones(size(tbl_scatter.ntrials(tbl_scatter.bin==1)));
clc

disp('%%%%%%%%%%%%%%%%% POLYNOMIAL ESTIMATE FOR ISO-MID-ORTHO %%%%%%%%%%%%%%%%%')
disp(' ')
% Iso vs. Mid
[cohen_d, pp, bf] = weighted_paired_ttest(iso, mid, weights);
fprintf('Iso vs. Mid   : Cohen''s d = %6.2f, p = %6.4f, mean difference = %6.2f, BF%s\n', cohen_d, pp, nanmean(iso - mid), bf);
% Iso vs. Ortho
[cohen_d, pp, bf] = weighted_paired_ttest(iso, ortho, weights);
fprintf('Iso vs. Ortho : Cohen''s d = %6.2f, p = %6.4f, mean difference = %6.2f, BF%s\n', cohen_d, pp, nanmean(iso - ortho), bf);
% Mid vs. Ortho
[cohen_d, pp, bf] = weighted_paired_ttest(mid, ortho, weights);
fprintf('Mid vs. Ortho : Cohen''s d = %6.2f, p = %6.4f, mean difference = %6.2f, BF%s\n', cohen_d, pp, nanmean(mid - ortho), bf);
disp(' ')


% Using bin size of 10:
tbl_scatter.bin = double(tbl_scatter.bin);
iso = tbl_scatter.bin_scatter(tbl_scatter.bin == 1);
mid = tbl_scatter.bin_scatter(tbl_scatter.bin == 2);
ortho = tbl_scatter.bin_scatter(tbl_scatter.bin == 3);
weights = ones(size(tbl_scatter.ntrials(tbl_scatter.bin==1)));

disp('%%%%%%%%%%%%%%%%%% DATA BIN ESTIMATE FOR ISO-MID-ORTHO %%%%%%%%%%%%%%%%%%')
disp(' ')
% Iso vs. Mid
[cohen_d, pp, bf] = weighted_paired_ttest(iso, mid, weights);
fprintf('Iso vs. Mid   : Cohen''s d = %6.2f, p = %6.4f, mean difference = %6.2f, BF%s\n', cohen_d, pp, nanmean(iso - mid), bf);
% Iso vs. Ortho
[cohen_d, pp, bf] = weighted_paired_ttest(iso, ortho, weights);
fprintf('Iso vs. Ortho : Cohen''s d = %6.2f, p = %6.4f, mean difference = %6.2f, BF%s\n', cohen_d, pp, nanmean(iso - ortho), bf);
% Mid vs. Ortho
[cohen_d, pp, bf] = weighted_paired_ttest(mid, ortho, weights);
fprintf('Mid vs. Ortho : Cohen''s d = %6.2f, p = %6.4f, mean difference = %6.2f, BF%s\n', cohen_d, pp, nanmean(mid - ortho), bf);
disp(' ')
diary off
%% Plot Results
% 
figure
hist(pol_deg,50)
set(gca,'FontSize',15)
xlabel('Best Degree (Lowest BIC)')
ylabel('Count')
saveas(gcf,'degree_histogram.eps')
saveas(gcf,'degree_histogram.png')

%% Define X-axis (stimulus dissimilarity index)
x = 0:90;
% Plot results
figure('Units','normalized','position',[.05 .4 .9 .4]);

for i = 1:3
    subplot(1,3,i)
    switch i
        case 1
            % All datasets
            out = zeros(size(all_curve,2),1);
            whichDatasets = 'All';
        case 2
            % Identify orientation datasets
            out = ~tbl_scatter.isoutlierDataset(1:3:end);
            whichDatasets = 'Orientation';
        case 3
            % Identify motion datasets
            out = tbl_scatter.isoutlierDataset(1:3:end);
            whichDatasets = 'Motion';
    end

    % Filter aggregated data based on outliers
    weighted_agg_scatt = aggregated_scatt(ismember(obsnum, find(~out)));
    weighted_agg_delta = aggregated_delta(ismember(obsnum, find(~out)));
    weighted_agg_counts = aggregated_weights(ismember(obsnum, find(~out)));
    aggregated_fit_curve = all_curve(:, ~out); aggregated_fit_curve = aggregated_fit_curve(:)';

    % Calculate average scatter for each delta
    all_scatter = grpstats(weighted_agg_scatt, weighted_agg_delta, 'mean');
    % Scatter plot of aggregated scatter data
    scatter(x, all_scatter(1:91),15*grpstats(weighted_agg_counts, weighted_agg_delta, 'mean'), 'filled', 'MarkerFaceColor', '#2B9AA9', 'MarkerFaceAlpha', .7);
    hold on;
    plot_mean_ci([repmat(0:90,1,sum(~out))],aggregated_fit_curve,1,1,'k');

    % Configure plot
    xlabel('Stimulus dissimilarity index');
    ylabel('Error Scatter');
    title([whichDatasets ' (' num2str(sum(~out)) '/' num2str(oo-1) ') code:obsid']);
    legend('Participant scatter average', '95% CI', 'Mean', 'Location', 'best');
    grid on;
    set(gca, 'FontSize', 15);
    xticks([0 45 90]);
end
saveas(gcf,'polyfit.eps')
saveas(gcf,'polyfit.png')


%% Visualization of results
figure('Units','normalized','position',[.05 .4 .9 .4]);
titles = {'iso', 'mid', 'ortho'};
for i = 1:3
    subplot(1, 3, i);
    poly = tbl_scatter.poly_scatter(i:3:end);
    binned = tbl_scatter.bin_scatter(i:3:end);
    idx = ~isnan(poly) & ~isnan(binned);
    [rr, pp] = corr(binned(idx), poly(idx));
    plot(binned(idx), poly(idx), 'ko');
    hold on;
    l1 = lsline;
    l1.Color = 'k';
    l1.LineWidth = 2;
    xlabel('Binned');
    ylabel('Polynomial Fit');
    title([titles{i}, ': r = ', num2str(round(rr, 2)), ', p = ', num2str(round(pp,3))]);
    set(gca, 'FontSize', 15);
end
saveas(gcf,'bin vs. polynomial.eps')
saveas(gcf,'bin vs. polynomial.png')

%% LMM s
delete('LME.txt')
diary('LME.txt')
tbl_scatter.bin = categorical(tbl_scatter.bin);
lmm             = fitlme(tbl_scatter,'poly_scatter ~ bin + (1+bin|obsid) + (1+bin|codenum) + (1+bin|codenum:obsid)');
disp(lmm)
tbl_scatter.new_bin = repmat(categorical([1 2 1]'),height(tbl_scatter)/3,1);
reduced_lmm = fitlme(tbl_scatter, 'poly_scatter ~ new_bin + (1+new_bin|obsid) + (1+new_bin|codenum) + (1+new_bin|codenum:obsid)');

bic_full = lmm.ModelCriterion.BIC;
bic_reduced = reduced_lmm.ModelCriterion.BIC;
dBIC = bic_reduced-bic_full;
BF_reduced_over_full = 1/exp((bic_reduced-bic_full)/2);
disp(['dBIC = ' num2str(dBIC) ', and BF_reduced_over_full = ' num2str(BF_reduced_over_full)])
diary off;

lmm             = fitlme(tbl_scatter,'bin_scatter ~ bin + (1+bin|obsid) + (1+bin|codenum) + (1+bin|codenum:obsid)');
disp(lmm)
tbl_scatter.new_bin = repmat(categorical([1 2 1]'),height(tbl_scatter)/3,1);
reduced_lmm = fitlme(tbl_scatter, 'bin_scatter ~ new_bin + (1+new_bin|obsid) + (1+new_bin|codenum) + (1+new_bin|codenum:obsid)');

bic_full = lmm.ModelCriterion.BIC;
bic_reduced = reduced_lmm.ModelCriterion.BIC;
dBIC = bic_reduced-bic_full;
BF_reduced_over_full = 1/exp((bic_reduced-bic_full)/2);

disp(['dBIC = ' num2str(dBIC) ', and BF_reduced_over_full = ' num2str(BF_reduced_over_full)])


