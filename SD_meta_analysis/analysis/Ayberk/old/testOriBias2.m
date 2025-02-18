toplot = tbl(strcmp(tbl.stimulus,'Orientation'),:);
toplot0 = toplot(ismember(abs(toplot.delta),[5]),:);
toplot90 = toplot(ismember(abs(toplot.delta),[15]),:);

figure
mv = sdp_mvav(toplot0.theta-90,toplot0.error_norm,21,0,[]).m; 
plot(0:180,mv,'Color','k','LineWidth',1);
hold on
mv = sdp_mvav(toplot90.theta-90,toplot90.error_norm,21,0,[]).m; 
plot(0:180,mv,'Color','r','LineWidth',1);

hold on
legend({'delta=15','delta=75'})
set(gca,'FontSize',15)
xlabel('Theta')
ylabel('Orientation bias')

%%
toplot = tbl(strcmp(tbl.stimulus,'Orientation'),:);
clear bias bias_all_theta 
num_delta_per_bin = 10;
bins = linspace(0,90,num_delta_per_bin+1);
color = [1 0 0].*linspace(0,1,length(bins))';
figure
% bins = [15 15] ;
for i  = 1:length(bins)-1
    thresholded_toplot = toplot(ismember(abs(toplot.delta),round([bins(i):ceil(bins(i+1))])),:);
    mv = sdp_mvav(thresholded_toplot.theta-90,thresholded_toplot.error_norm,21,0,[]).m;
    bias(i,:) = mv;
    bias_all_theta(i) = max(mv);
    plot(0:180,mv,'Color',color(i,:))
    hold on
end
figure
plot(bias_all_theta)

%% 
toplot = tbl(strcmp(tbl.stimulus,'Orientation'),:);
iso = toplot(abs(toplot.delta)<=10,:);
ortho = toplot(abs(toplot.delta)>=80,:);
figure
iso_ori_bias = sdp_mvav(iso.theta-90,iso.error_norm,21,0,[]).m;
ortho_ori_bias = sdp_mvav(ortho.theta-90,ortho.error_norm,21,0,[]).m;
% Calculate observed difference
% Element-wise differences
diff_ori_bias = ortho_ori_bias - iso_ori_bias;
th=0:180;

% Observed mean difference
oribiaspeaks =  toplot(ismember(toplot.theta,[7:14 ]),:);
nooribiaspeaks =  toplot(ismember(toplot.theta,[45+ [7:14]]),:);

plot_mean_ci(oribiaspeaks.delta,oribiaspeaks.error_iqr_norm,3,21,'k');
plot_mean_ci(nooribiaspeaks.delta,nooribiaspeaks.error_iqr_norm,3,21,'r');


% Bootstrapping
n_bootstrap = 10000; % Number of bootstrap iterations
bootstrap_diffs = zeros(1, n_bootstrap);

% Bootstrapping
for i = 1:n_bootstrap
    % Resample with replacement
    resampled_diff = diff_ori_bias(randi(length(diff_ori_bias), 1, length(diff_ori_bias)));
    
    % Calculate the mean of the resampled differences
    bootstrap_diffs(i) = mean(resampled_diff);
end

% Calculate confidence intervals
ci_lower = prctile(bootstrap_diffs, 2.5); % Lower 2.5% percentile
ci_upper = prctile(bootstrap_diffs, 97.5); % Upper 97.5% percentile

% Determine significance
is_significant = observed_diff < ci_lower || observed_diff > ci_upper;

% Plot bootstrap distribution
figure;
histogram(bootstrap_diffs, 50, 'Normalization', 'probability');
hold on;
xline(observed_diff, 'r', 'LineWidth', 2); % Observed difference
xline(ci_lower, 'b--', 'LineWidth', 2); % Lower confidence bound
xline(ci_upper, 'b--', 'LineWidth', 2); % Upper confidence bound
hold off;

% Add labels and title
xlabel('Bootstrapped Differences');
ylabel('Probability');
title('Bootstrap Distribution of Differences');
legend('Bootstrap Distribution', 'Observed Difference', '95% CI Bounds');

% Output result
if is_significant
    disp('The observed difference is statistically significant.');
else
    disp('The observed difference is not statistically significant.');
end
%%
thrs = 0:45:90;%:10:20;
color = [1 0 0].*linspace(0,1,length(thrs))';
figure
subplot(121)
mv = movmean(repmat(grpstats(toplot1.error_norm,toplot1.delta,'mean'),3,1),35);
mv = mv(180:360); mv = mean([-mv(1:91) mv(end:-1:91)],2); mv = [-mv; mv(end-1:-1:1)];
plot(-90:90,mv,'Color','b','LineWidth',1.5);hold on;
subplot(122)
mv = movmean(repmat(grpstats(toplot1.error_norm,toplot1.delta,'std'),3,1),35);
mv = mv(180:360); mv = mean([mv(1:91) mv(end:-1:91)],2); mv = [mv; mv(end-1:-1:1)];
plot(-90:90,mv,'Color','b','LineWidth',1.5);hold on;
cc = 1;
for d=thrs
    toplot2=toplot(abs(toplot.delta)==d,:);
    subplot(121)
    
    % hold on
    mv = movmean(repmat(grpstats(toplot2.error_norm_1back,toplot2.delta_1back,'mean'),3,1),35);
    mv = mv(180:360); mv = mean([-mv(1:91) mv(end:-1:91)],2); mv = [-mv; mv(end-1:-1:1)];
    plot(-90:90,mv,'Color',color(cc,:));
    hold on
    subplot(122)
    % mv = movmean(repmat(grpstats(toplot1.error_norm,toplot1.delta,'std'),3,1),35);
    % mv = mv(180:360); mv = mean([mv(1:91) mv(end:-1:91)],2); mv = [mv; mv(end-1:-1:1)];
    % plot(-90:90,mv,'Color','k');
    % hold on
    mv = movmean(repmat(grpstats(toplot2.error_norm_1back,toplot2.delta_1back,'std'),3,1),35);
    mv = mv(180:360); mv = mean([mv(1:91) mv(end:-1:91)],2); mv = [mv; mv(end-1:-1:1)];
    plot(-90:90,mv,'Color',color(cc,:));
    hold on
    cc=cc+1;
end