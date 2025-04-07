clear all;
close all;
clc;
load(['..' filesep '..' filesep 'data' filesep 'datasets' filesep 'SD_ma_master_table.mat'])
tbl.theta(tbl.theta==180) = 0;
tbl.delta(tbl.delta==-90) = 90;
nBoot = []; % if empty, gives true mean and non-smoothed confidence intervals
rep_size = 3;
bin_size = 21;
%% summarize datasets:
datasets     = unique(tbl.studynum);  % Get unique dataset codes
n_datasets   = numel(datasets);
variable     = 'error_ori_deb_sd_deb';  % Variable for which to compute trial counts
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
writetable(report_table,'summary_studies.csv')


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
ylabel('stimulus-specific bias (°)')

subplot(122)
p0 = plot_mean_ci(tmp.delta,tmp.error_norm,rep_size,bin_size,'k','mean');
hold on;
p1 = plot_mean_ci(tmp.delta,tmp.error_iqr_norm,rep_size,bin_size,'r','mean');
p2 = plot_mean_ci(tmp.delta,tmp.error_ori_deb_norm,rep_size,bin_size,'b','mean');
p3 = plot_mean_ci(tmp.delta,tmp.error_ori_deb_sd_deb_norm,rep_size,bin_size,'g','mean');
ylim([-.5 .5])
xticks(-90:45:90)
xlabel('Δ (°)')
ylabel('serial dependence bias (°)')
legend([p0 p1 p2 p3],{'raw data' 'outliers removed' 'stimulus bias removed' 'serial dependence bias removed'}, 'Location', 'best')

set(gcf,'PaperOrientation','landscape')
set(gcf, 'PaperUnits', 'normalized');
exportgraphics(gcf,['figures' filesep 'SI_removalChecks.pdf'],'BackgroundColor','none')

%% Show ortho vs. iso ori biases
clear superiority sd_peak cohenD_mean

n_bins = 10;
intersect_thetas = [0:181];
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

tmp    = tbl(abs(tbl.delta)<=90 & tbl.theta<=180,:);
rep_size = 3;
bin_size = 31;
for i = 0:180
    dd(i+1) = nanmean(tmp.error_iqr_norm(tmp.theta == i))/nanstd(tmp.error_iqr_norm(tmp.theta == i));
end
[abs_dd_sorted,idx]=sort(abs(dd));


figure('Units','normalized','position',[0 .4 .5 .9]);
subplot(6,3,1:3)
p1 = plot_mean_ci(tmp.theta,tmp.error_iqr_norm,rep_size,bin_size,'k','mean');
xlabel('θ (°)')
ylabel('Stimulus-specific Bias (°)')
ylim([-.6 .6])
xticks(0:45:180)
for b = 1:n_bins
    subplot(3,6,1:3)
    highlight_idx = intersect(idx((b-1)*180/n_bins+1:b*180/n_bins),intersect_thetas);
    h1=plot(highlight_idx-1,zeros(size(highlight_idx)),'Marker','square','Color',pastel_cmap(b,:),'MarkerSize',10,'MarkerFaceColor',pastel_cmap(b,:),'LineStyle','none');
    tt = tmp(ismember(abs(tmp.theta),highlight_idx-1),:);
    subplot(3,6,7:9)
    p1m = plot_mean_ci(tt.delta,tt.error_iqr_norm,rep_size,bin_size,pastel_cmap(b,:),'mean');
    xlabel('Δ (°)')
    ylabel('Serial Dependence Bias (°)')
    ylim([-.6 .6])
    xticks(-90:45:90);

    subplot(3,6,10:12)
    [p1s,f1s] = plot_mean_ci(tt.delta,tt.error_iqr_norm,rep_size,bin_size,pastel_cmap(b,:),'std');
    xlabel('Δ (°)')
    ylabel('Error Scatter (°)')
    ylim([.8 1.2])
    xticks(-90:45:90);
    % title(['Mean Error Scatter = ' num2str(round(mean(p1s.YData),2))])

    sd_peak(b) = max(p1m.YData);
    superiority(b) = (p1s.YData(1)-p1s.YData(91));
    cohenD_mean(b) = nanmean(abs_dd_sorted((b-1)*180/n_bins+1:b*180/n_bins));

    subplot(3,6,13:14)
    scatter(cohenD_mean(b),sd_peak(b),400,'filled','o','MarkerFaceColor',pastel_cmap(b,:),'MarkerFaceAlpha',.8);hold on;

    subplot(3,6,15:16)
    scatter(cohenD_mean(b),superiority(b),400,'filled','o','MarkerFaceColor',pastel_cmap(b,:),'MarkerFaceAlpha',.8);hold on;

    subplot(3,6,17:18)
    scatter(sd_peak(b),superiority(b),400,'filled','o','MarkerFaceColor',pastel_cmap(b,:),'MarkerFaceAlpha',.8);hold on;
end
subplot(3,6,4:6)
p1 = plot_mean_ci(tmp.theta,tmp.error_iqr_norm,rep_size,bin_size,'k','std');
xlabel('θ (°)')
ylabel('Error Scatter (°)')
ylim([.8 1.2])
xticks(0:45:180);
% title(['Mean Error Scatter = ' num2str(round(mean(p1.YData),2))])

subplot(3,6,13:14)
xlabel('Cohen''s d for Stimulus-specific Bias (°)')
ylabel('Serial Dependence Bias Peak (°)')
ll1 = plot(cohenD_mean, polyval(polyfit(cohenD_mean, sd_peak, 1), cohenD_mean),'Color',[.5 .5 .5], 'LineWidth', 2);
[rr,pp] = corr(cohenD_mean',sd_peak'); legend(ll1,['r=' num2str(round(rr,2))],'Location','northeast')
set(gca, 'FontSize', 20, 'LineWidth', 2, 'FontName', 'Times New Roman');
grid on;
box on;

subplot(3,6,15:16)
xlabel('Cohen''s d for Stimulus-specific Bias (°)')
ylabel('Superiority (ortho-iso) (°)')
ll1 = plot(cohenD_mean, polyval(polyfit(cohenD_mean, superiority, 1), cohenD_mean),'Color',[.5 .5 .5], 'LineWidth', 2);
[rr,pp] = corr(cohenD_mean',superiority'); legend(ll1,['r=' num2str(round(rr,2)) ],'Location','southeast')
set(gca, 'FontSize', 20, 'LineWidth', 2, 'FontName', 'Times New Roman');
grid on;
box on;

subplot(3,6,17:18)
xlabel('Serial Dependence Bias Peak (°)')
ylabel('Superiority (ortho-iso) (°)')
ll1 = plot(sd_peak, polyval(polyfit(sd_peak, superiority, 1), sd_peak),'Color',[.5 .5 .5], 'LineWidth', 2);
[rr,pp] = corr(sd_peak',superiority'); legend(ll1,['r=' num2str(round(rr,2)) ],'Location','southeast')
set(gca, 'FontSize', 20, 'LineWidth', 2, 'FontName', 'Times New Roman');
grid on;
box on;

set(gcf,'PaperOrientation','landscape')
set(gcf, 'PaperUnits', 'normalized');
exportgraphics(gcf, ['figures' filesep 'SI_oriBias_scatter.pdf'], 'BackgroundColor', 'none', 'ContentType', 'vector');


%% Show bias cleaning effects on scatter (as done in Fritsche et al. 2020)
figure('Units','normalized','position',[0 0 1 1]);
tmp = tbl(abs(tbl.delta)<=90,:);
tmp.theta= round(tmp.theta);
tmp.delta = round(tmp.delta);
datasets = unique(tmp.codenum);
bin_size = 21;
rep_size = 3;

% Define the ranges of theta and delta
theta_range = 0:360;
delta_range = -90:90;
% Create a table of all unique combinations of theta and delta
[theta_grid, delta_grid] = ndgrid(theta_range, delta_range); % Create all combinations
T = table(); 

for i = 1:length(datasets)
    subplot(7,7,i)

    tt = tmp(tmp.codenum==datasets(i),:);
    obsids = unique(tt.obsid);
    scatt_subj_iqr = nan(91,length(obsids));
    scatt_subj_ori_deb = nan(91,length(obsids));
    % for th = 0:179
    %     for dl = -90:90
    %         T = [T ; [i th dl sum(tt.theta==th & tt.delta==dl)]];
    %     end
    % end

    mytbl=table();
    mytbl.dataset = i*ones(size(theta_grid(:)));
    mytbl.theta = theta_grid(:);
    mytbl.delta = delta_grid(:);
    counts = accumarray([tt.theta + 1, tt.delta + 91], 1, [361, 181]);
    mytbl.count = counts(:); % 1-based indexing
    % Use accumarray to efficiently count occurrences
    T = [T;mytbl];
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
    xlabel('Δ (°)')
    ylabel('Error Scatter (°)')
    xticks(0:45:90)
    plot([0 90], [0 0],'k--','LineWidth',.5)
    grid off;
    set(gca,'FontSize',11)
end
writetable(T,'counts.csv')
set(gcf,'PaperOrientation','landscape')
set(gcf, 'PaperUnits', 'normalized');
exportgraphics(gcf, ['figures' filesep 'SI_allDatasetsScatter.pdf'], 'BackgroundColor', 'none', 'ContentType', 'vector');

%% Initialize variables for polynomial fitting
nanunique = @(x) unique(x(~isnan(x)));
tbl_scatter = [];
all_curve = [];
pol_deg = [];
aggregated_delta = [];
aggregated_scatt = [];
aggregated_counts = [];
obsnum = [];
isOrientationDataset = find(grpstats(tbl(:,{'codenum','stimtype'}),{'codenum'},'mean').mean_stimtype);
oo = 1;
error_variable = 'error_ori_deb'; % error_iqr   error_ori_deb   error_iqr_norm   error_ori_deb_norm   
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
%% Do paired t-test on iso-mid-ortho
% Start recording to a text file
delete('paired t-test.txt')
diary('paired t-test.txt');
% using polynomial:
tbl_scatter.SI =double(tbl_scatter.SI);
iso = tbl_scatter.ES(tbl_scatter.SI==1);
mid = tbl_scatter.ES(tbl_scatter.SI==2);
ortho = tbl_scatter.ES(tbl_scatter.SI==3);
weights = ones(size(tbl_scatter.ntrials(tbl_scatter.SI==1)));
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
tbl_scatter.SI = double(tbl_scatter.SI);
iso = tbl_scatter.bin_scatter(tbl_scatter.SI == 1);
mid = tbl_scatter.bin_scatter(tbl_scatter.SI == 2);
ortho = tbl_scatter.bin_scatter(tbl_scatter.SI == 3);
weights = ones(size(tbl_scatter.ntrials(tbl_scatter.SI==1)));

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
set(gcf,'PaperOrientation','landscape')
set(gcf, 'PaperUnits', 'normalized');
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
    weighted_agg_counts = aggregated_counts(ismember(obsnum, find(~out)));
    aggregated_fit_curve = all_curve(:, ~out); aggregated_fit_curve = aggregated_fit_curve(:)';

    % Calculate average scatter for each delta
    all_scatter = grpstats(weighted_agg_scatt, weighted_agg_delta, 'mean');
    % Scatter plot of aggregated scatter data
    scatter(unique(weighted_agg_delta), all_scatter,15*grpstats(weighted_agg_counts, weighted_agg_delta, 'mean'), 'filled', 'MarkerFaceColor', '#2B9AA9', 'MarkerFaceAlpha', .7);
    hold on;
    p=plot_mean_ci([repmat(0:90,1,sum(~out))],aggregated_fit_curve,1,1,'k','mean',[]);
    p.LineWidth = 3;

    % Configure plot
    xlabel('|Δ| (°)');
    ylabel('Error Scatter (°)');
    title([whichDatasets ' Datasets']); % ' (' num2str(sum(~out)) '/' num2str(oo-1) ') code:obsid'
    xticks([0 45 90]);
    % xticklabels({'iso' 'mid' 'ortho'});
end
set(gcf,'PaperOrientation','landscape')
set(gcf, 'PaperUnits', 'normalized');
exportgraphics(gcf,['figures' filesep 'main_results_polyfit.pdf'],'BackgroundColor','none')


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
    legend(l1, ['r = ', num2str(round(rr, 2)),],'Location','best')
    set(gca, 'FontSize', 24);
    set(gca,'LineWidth',2)
    set(gca,'FontName','Times New Roman')
    grid on;
    box on;
end
set(gcf,'PaperOrientation','landscape')
set(gcf, 'PaperUnits', 'normalized');
exportgraphics(gcf,['figures' filesep 'SI_polyBinCorrelations.pdf'],'BackgroundColor','none')
%% LMM s
delete('LME.txt')
diary('LME.txt')
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


