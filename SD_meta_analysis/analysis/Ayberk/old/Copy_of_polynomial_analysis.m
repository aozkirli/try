clear all;
close all;
clc;
load(['..' filesep '..' filesep 'data' filesep 'datasets' filesep 'SD_ma_master_table.mat'])
tbl.theta(tbl.theta==180) = 0;
tbl.delta(tbl.delta==-90) = 90;

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

%% parameters for plot:
nBoot = []; % if empty, gives true mean and non-smoothed confidence intervals
rep_size = 3;
bin_size = 21;
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
%% simulation for the 
figure('Units','normalized','position',[.05 .4 .6 .4]);
subplot(121)
tmp = tbl(strcmp(tbl.stimulus, 'Orientation'), :);
p1 = plot_mean_ci(tmp.theta,tmp.error_iqr_norm,rep_size,bin_size,'k','mean');
p2 = plot_mean_ci(tmp.theta,tmp.error_iqr_norm,rep_size,bin_size,'k','std');

ori_bias = p1.YData;
ori_bias_std = p2.YData;
clear counts;
dd = 1;
for d = 0:90
    tt = 1;
    for th = 0:179
        counts(dd,tt)=sum(tmp.theta==th & abs(tmp.delta) ==d);
        tt = tt+1;
    end
    dd = dd+1;
end
first_bin = [0:10:10];
second_bin = [80:10:90];

first_bin_counts = sum(counts(first_bin+1,:));
second_bin_counts = sum(counts(second_bin+1,:))/((length(second_bin)/length(first_bin))); % scaled to control for #delta
subplot(122)
b1=bar(0:179,first_bin_counts,'stacked','red','FaceAlpha',.5,'EdgeColor','none');
hold on;
b2=bar(0:179,second_bin_counts,'stacked','black','red','FaceAlpha',.5,'EdgeColor','none');
legend([b1 b2],{['Delta in [' num2str(first_bin(1)) ', ' num2str(first_bin(end)) ']'], ['Delta in [' num2str(second_bin(1)) ', ' num2str(second_bin(end)) ']']}, 'Location', 'northeast');

for nb = 1:1000
    sample1 = [];
    for i = 1:180
        sample1 = [sample1 ori_bias(i)+ori_bias_std(i)*randn(1,first_bin_counts(i))];
    end
    sample2 = [];
    for i = 1:180
        sample2 = [sample2 ori_bias(i)+ori_bias_std(i)*randn(1,ceil(second_bin_counts(i)))];
    end
    std1(nb)=std(sample1);
    std2(nb)=std(sample2);
end
[hh,~,~,stats]=ttest(std2,std1);
superiority = (std2-std1);
disp(['Ortho - iso : BF_{10} = ' num2str(t1smpbf(stats.tstat,180)) ' , Cohen''s d = ' num2str(computeCohen_d(std2,std1,'paired')) ...
    ', mean difference = ' num2str(mean(superiority))]);
%%
std1 = nan(1000,91);
for nb=1:1000
    for d=1:91
    sample1 = [];
    for i = 1:180
        sample1 = [sample1 ori_bias(i)+.0*randn(1,counts(d,i))];
    end
    std1(nb,d)=std(sample1);
    end
end
[hh,~,~,stats]=ttest(std2,std1);
superiority = (std2-std1);
disp(['Ortho - iso : BF_{10} = ' num2str(t1smpbf(stats.tstat,180)) ' , Cohen''s d = ' num2str(computeCohen_d(std2,std1,'paired')) ...
    ', mean difference = ' num2str(mean(superiority))]);
%% 
figure('Units', 'normalized', 'Position', [0.05 0.4 0.3 0.4]);
nbins = 90;
% leave one delta (10*k) out to avoid overestimation for either bin (0 and 90 has many more trials due to binned datasets)
first_bin = [0:10:20];
second_bin = [30:10:90];
% Filter and process data
tmp = tbl(strcmp(tbl.stimulus, 'Orientation'), :);
% delta -90 and 90 are the same, theta 180 and 0 are the same
tmp.theta(tmp.theta == 180) = 0;
tmp.delta(tmp.delta == -90) = 90;

hfirst= histogram(tmp.theta(ismember(abs(tmp.delta), first_bin)),nbins, 'FaceColor', 'red', 'FaceAlpha', .5, 'EdgeColor', 'none');
hold on;
hsecond=histogram(tmp.theta(ismember(abs(tmp.delta), second_bin)),nbins, 'FaceColor', 'black', 'FaceAlpha', .5, 'EdgeColor', 'none');
% Scale the second histogram counts
hsecond.BinCounts = hsecond.BinCounts / (length(second_bin)/length(first_bin));
ylim([min([hfirst.BinCounts,hsecond.BinCounts])/1.1 max([hfirst.BinCounts,hsecond.BinCounts])*1.1])

% Add labels and title
xlabel('Theta');
ylabel('Count');
title('Overlapping Histograms for Different Delta Ranges');
legend([hfirst hsecond],{['Delta in [' num2str(first_bin(1)) ', ' num2str(first_bin(end)) ']'], ['Delta in [' num2str(second_bin(1)) ', ' num2str(second_bin(end)) ']']}, 'Location', 'northeast');
hold off;
set(gca,'FontSize',15)

%% Show ortho vs. iso ori biases
tmp    = tbl(abs(tbl.delta)<=90 & tbl.theta<=180,:);
% tmp    = tbl(strcmp(tbl.stimulus,'Orientation'),:);
iso    = tmp(ismember(abs(tmp.delta),0:10),:);
ortho  = tmp(ismember(abs(tmp.delta),80:90),:);
rep_size = 3;
bin_size = 21;
for i = 0:180
    dd(i+1) = nanmean(tmp.error_iqr_norm(tmp.theta == i))/nanstd(tmp.error_iqr_norm(tmp.theta == i));
end
[abs_dd_sorted,idx]=sort(abs(dd));
n_bins = 4;

veryStrongOriBiasTheta = idx(181-180/n_bins:181)-1;
veryWeakOriBiasTheta = idx(1:180/n_bins)-1;
% veryStrongOriBiasTheta = setdiff(0:179,[5:20 90+[5:20] 90-[5:20] 180-[5:20]]);%intersect(veryStrongOriBiasTheta,[0:45 90:135]);
% veryWeakOriBiasTheta = [5:20 90+[5:20] 90-[5:20] 180-[5:20]];%setdiff(veryWeakOriBiasTheta,[40:50 130:140]);

figure('Units','normalized','position',[0 .4 1 .3]);
subplot(141)
p1 = plot_mean_ci(tmp.theta,tmp.error_iqr_norm,rep_size,bin_size,'k','mean');
highlight = nan(181,1); highlight(veryStrongOriBiasTheta+1) = 0;
h1=plot(0:180,highlight,'Color','#FF7F0E','LineWidth',5);
highlight = nan(181,1); highlight(veryWeakOriBiasTheta+1) = 0;
h2=plot(0:180,highlight,'Color','#17BECF','LineWidth',5);
xlabel('theta')
ylabel('stimulus-specific bias')
legend([h1 h2],{'stronger' 'weak'}, 'Location', 'northeast')
ylim([-.6 .6])

% check sort the data based on cohen's d for orientation bias:
strong     = tmp(ismember(abs(tmp.theta),veryStrongOriBiasTheta),:);
weak  = tmp(ismember(abs(tmp.theta),veryWeakOriBiasTheta),:);

subplot(142)
p1m = plot_mean_ci(strong.delta,strong.error_iqr_norm,rep_size,bin_size,'#FF7F0E','mean');
hold on;
p2m = plot_mean_ci(weak.delta,weak.error_iqr_norm,rep_size,bin_size,'#17BECF','mean');
xlabel('delta')
ylabel('serial dependence bias')
% legend([p1 p2],{'strong shift' 'no shift'}, 'Location', 'northeast')
ylim([-.6 .6])

subplot(143)
p1s = plot_mean_ci(strong.delta,strong.error_iqr_norm,rep_size,bin_size,'#FF7F0E','std');
hold on;
p2s = plot_mean_ci(weak.delta,weak.error_iqr_norm,rep_size,bin_size,'#17BECF','std');
xlabel('delta')
ylabel('error scatter')
ylim([.8 1.2])

subplot(144)
tmp    = tbl(abs(tbl.delta)<=90,:);
var_types       = varfun(@class, tmp, 'OutputFormat', 'cell');
num_log_vars    = ismember(var_types, {'double', 'logical'});
filtered_tbl    = tmp(:, num_log_vars);
tmp_scatter     = grpstats(filtered_tbl, {'codenum' 'obsid' 'bin'}, 'std');
tbl_mean          = grpstats(filtered_tbl, {'codenum' 'obsid'}, 'mean');
tbl_mean.sup      = tmp_scatter.std_error_iqr_norm(tmp_scatter.bin==3) - tmp_scatter.std_error_iqr_norm(tmp_scatter.bin==1);
tbl_datasets_mean = grpstats(tbl_mean, {'codenum'}, 'mean');
colors            = linspecer(length(unique(tmp.codenum)));%
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

set(gcf,'PaperOrientation','landscape')
set(gcf, 'PaperUnits', 'normalized');
exportgraphics(gcf, ['figures' filesep 'SI_oriBiasVSscatter.pdf'], 'BackgroundColor', 'none', 'ContentType', 'vector');


%% Show bias cleaning effects on scatter (as done in Fritsche et al. 2020)
figure('Units','normalized','position',[0 0 1 1]);
tmp = tbl(abs(tbl.delta)<=90,:);
datasets = unique(tmp.codenum);
bin_size = 21;
rep_size = 3;
fold = 0;
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
    xlabel('Δ (°)')
    ylabel('Error Scatter (°)')
    xticks(0:45:90)
    plot([0 90], [0 0],'k--','LineWidth',.5)
    grid off;
    set(gca,'FontSize',11)
end

set(gcf,'PaperOrientation','landscape')
set(gcf, 'PaperUnits', 'normalized');
exportgraphics(gcf, ['figures' filesep 'SI_allDatasetsScatter.pdf'], 'BackgroundColor', 'none', 'ContentType', 'vector');


%% heatmaps from David
% orientation
tbl_ori    = tbl(strcmp(tbl.stimulus,'Orientation'),:);
x = -90:90;
sldwin              = 10;
scatter_ori_iqr         = nan(numel(x),numel(x));
bias_ori_iqr            = scatter_ori_iqr;
scatter_ori_ori_deb = scatter_ori_iqr;
bias_ori_ori_deb = scatter_ori_iqr;
for k = 1:numel(x)
    k
    idx             = mod(x(k)+[-sldwin:sldwin],180);
    tmp             = tbl_ori(ismember(tbl_ori.theta,idx),:);
    mv              = sdp_mvav(tmp.delta,tmp.error_iqr_norm,sldwin*2+1,2);
    scatter_ori_iqr(k,:)= mv.std;
    bias_ori_iqr(k,:)   = mv.m;

    mv              = sdp_mvav(tmp.delta,tmp.error_ori_deb_norm,sldwin*2+1,2);
    scatter_ori_ori_deb(k,:)= mv.std;
    bias_ori_ori_deb(k,:)   = mv.m;
end
%%
% cmap = [linspace(0, 1, 50)', linspace(0, 1, 50)', ones(50, 1); ...
% ones(50, 1), linspace(1, 0, 50)', linspace(1, 0, 50)'];
cmap='parula'
figure('Units','normalized','Position',[0 0 .5 .6])
subplot(221)
imagesc(-90:90, 0:180, bias_ori_iqr-bias_ori_iqr(:,91));
colorbar; colormap(cmap);
xlabel('Delta'); ylabel('Theta'); title('Bias (error_iqr_norm)','Interpreter','none');
% clim([-1, 1]);
set(gca,'FontSize',15);

subplot(222)
imagesc(-90:90, 0:180, 100*scatter_ori_iqr./mean(scatter_ori_iqr(:,80:100),2)-100);
xlabel('Delta'); ylabel('Theta'); title('Scatter (error_iqr_norm)','Interpreter','none');
colorbar; colormap(cmap);
clim([-12 12]);
set(gca,'FontSize',15);

subplot(223)
imagesc(-90:90, 0:180, bias_ori_ori_deb-bias_ori_ori_deb(:,91));
colorbar; colormap(cmap);
xlabel('Delta'); ylabel('Theta'); title('Bias (error_ori_deb_norm)','Interpreter','none');
% clim([-1, 1]);
set(gca,'FontSize',15);

subplot(224)
imagesc(-90:90, 0:180, 100*scatter_ori_ori_deb./mean(scatter_ori_ori_deb(:,80:100),2)-100);
% imagesc(-90:90, 0:180, 100*scatter_ori_ori_deb./scatter_ori_ori_deb(:,91)-100);
xlabel('Delta'); ylabel('Theta'); title('Scatter (error_ori_deb_norm)','Interpreter','none');
colorbar; colormap(cmap);
% clim([.7 1.2]);
set(gca,'FontSize',15);
clim([-12 12]);

figure
imagesc(-90:90, 0:180,(100*scatter_ori_ori_deb./mean(scatter_ori_ori_deb(:,80:100),2)-100)- (100*scatter_ori_iqr./mean(scatter_ori_iqr(:,80:100),2)-100))
colorbar; colormap(cmap);
xlabel('Delta'); ylabel('Theta'); title('Decrease in superiority after orientation bias removal','Interpreter','none');
clim([-6, 6]);
set(gca,'FontSize',15);

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


