close all;

% Load data
load(['..' filesep '..' filesep 'data' filesep 'datasets' filesep 'SD_ma_master_table.mat']);
tbl.theta(tbl.theta == 180) = 0;
tbl.delta(tbl.delta == -90) = 90;

% Filter and preprocess data
toplot = tbl(strcmp(tbl.stimulus, 'Orientation'), :);
toplot.theta = round(toplot.theta);
toplot.delta = round(toplot.delta);
toplot.theta90 = toplot.theta - 90;

% Line width for plots
lw = 2.5;

% Figure setup
figure('Units', 'normalized', 'Position', [0 0 0.5 1]);

% Plot stimulus-specific bias
subplot(3, 2, 1);
hold on;
plot([0 0], [-1 1], 'k--', 'LineWidth', 1);
plot([0 180], [0 0], 'k--', 'LineWidth', 1);

% Stimulus-specific bias computation and plot
dd = zeros(180, 1);
for i = 0:179
    [hh, pp] = ttest(toplot.error_iqr_norm(toplot.theta == i), nanmean(toplot.error_iqr_norm));
    dd(i + 1) = nanmean(toplot.error_iqr_norm(toplot.theta == i)) / nanstd(toplot.error_iqr_norm(toplot.theta == i));
    eb(i + 1) = errorbar(i, nanmean(toplot.error_iqr_norm(toplot.theta == i)), ...
        nanstd(toplot.error_iqr_norm(toplot.theta == i)) / sqrt(sum(toplot.theta == i)), ...
        'LineWidth', lw, 'Color', 'b');
end
xlabel('θ(°)');
ylabel('Stimulus-specific Bias (°)');
set(gca, 'FontSize', 20, 'LineWidth', 2, 'FontName', 'Times New Roman');
grid on;
box on;
xticks(0:45:180);

% Plot error scatter
subplot(3, 2, 2);
mv = nan(180, 1);
idx = unique(abs(toplot.theta)) + 1;
idx = idx(isfinite(idx));
mv(idx) = grpstats(toplot.error_iqr_norm, toplot.theta, 'std');
mv = movmean(repmat(mv, 3, 1), 31);
mv = mv(181:360);
plot(0:179, mv, 'LineWidth', lw, 'Color', 'k');
xlabel('θ (°)');
ylabel('Error Scatter (°)');
set(gca, 'FontSize', 20, 'LineWidth', 2, 'FontName', 'Times New Roman');
grid on;
box on;
xticks(0:45:180);

% Threshold processing and additional plots
abs_dd = abs(dd);
[abs_dd_sorted, idx] = sort(abs_dd);
num_theta_per_bin = 36;
thrs = abs_dd_sorted(num_theta_per_bin:num_theta_per_bin:min([length(abs_dd_sorted) num_theta_per_bin*ceil(length(abs_dd_sorted)/num_theta_per_bin+1) ]));
color = [1 0 0].*linspace(0,1,length(thrs))';
cc = 1;
out_theta = [];
for thr = thrs'
    % Identify thresholds and filter data
    if thr == thrs(1)
        card = find(abs(dd) <= thr) - 1;
        prev_card = card;
    else
        card = setdiff(find(abs(dd) <= thr) - 1, prev_card);
        prev_card = find(abs(dd) <= thr) - 1;
    end
    theseThetas = setdiff(card, out_theta);
    thresholded_toplot = toplot(ismember(toplot.theta, theseThetas), :);

    % Highlight selected data points
    for i = theseThetas'
        eb(i + 1).Color = color(cc, :);
    end

    % Serial dependence bias
    subplot(3, 2, 3);
    mv = nan(91, 1);
    idx = unique(abs(thresholded_toplot.delta)) + 1;
    idx = idx(isfinite(idx));
    error_iqr = thresholded_toplot.error_iqr_norm;
    error_iqr(thresholded_toplot.delta < 0) = -error_iqr(thresholded_toplot.delta < 0);
    mv(idx) = grpstats(error_iqr, abs(thresholded_toplot.delta), 'mean');
    mv = [-mv(end:-1:2); mv];
    mv = movmean(repmat(mv, 3, 1), 31);
    mv = mv(182:362);
    bias(cc) = max(mv);
    plot(-90:90, mv, 'LineWidth', lw, 'Color', color(cc, :));
    hold on;

    if cc == 1
        xlabel('Δ (°)');
        ylabel('Serial Dependence Bias (°)');
        set(gca, 'FontSize', 20, 'LineWidth', 2, 'FontName', 'Times New Roman');
        grid on;
        box on;
    end

    % Error scatter
    subplot(3, 2, 4);
    mv = nan(91, 1);
    mv(idx) = grpstats(thresholded_toplot.error_iqr_norm, abs(thresholded_toplot.delta), 'std');
    mv = [mv(end:-1:2); mv];
    mv = movmean(repmat(mv, 3, 1), 31);
    mv = mv(182:362);
    plot(-90:90, mv, 'LineWidth', lw, 'Color', color(cc, :));
    hold on;
    sup(cc) = mv(end)-mv(91);
    set(gca, 'FontSize', 20, 'LineWidth', 2, 'FontName', 'Times New Roman');
        grid on;
        box on;

    % Update metrics and correlations
    if cc > 1
        subplot(3, 2, 5);
        sc1.XData = [sc1.XData mean(abs_dd(abs_dd > thrs(cc - 1) & abs_dd <= thrs(cc)))];
        sc1.YData = sup(1:cc);
        ll1.XData = sc1.XData;
        ll1.YData = polyval(polyfit(sc1.XData, sup(1:cc), 1), sc1.XData);

        sc2.XData = sc1.XData;
        sc2.YData = bias(1:cc);
        ll2.XData = sc1.XData;
        ll2.YData = polyval(polyfit(sc1.XData, bias(1:cc), 1), sc1.XData);

    else
        % Initial scatter plot setup
        subplot(3, 2, 5);
        sc1 = scatter(mean(abs_dd(abs_dd <= thrs(cc))), sup(cc), 'filled', 'k');
        hold on;
        ll1 = plot(thrs(1:cc), polyval(polyfit(sc1.XData, sup(1:cc), 1), thrs(1:cc)), 'r', 'LineWidth', lw);
        xlabel('Stimulus-specific Bias Effect Size (d)');
        ylabel('Superiority (ortho-iso) (°)');
        set(gca, 'FontSize', 20, 'LineWidth', 2, 'FontName', 'Times New Roman');
        grid on;
        box on;

        % Update bias plot
        subplot(3, 2, 6);
        sc2 = scatter(sc1.XData, bias(cc), 'filled', 'k');
        hold on;
        ll2 = plot([sc2.XData], polyval(polyfit([sc2.XData], bias(1:cc), 1), [sc2.XData]), 'r', 'LineWidth', lw);
        xlabel('Stimulus-specific Bias Effect Size (d)');
        ylabel('Serial Dependence Bias Peak (°)');
        set(gca, 'FontSize', 20, 'LineWidth', 2, 'FontName', 'Times New Roman');
        grid on;
        box on;

    end
    cc = cc + 1;
    pause(0.0001);
end

subplot(3,2,5)
[rr,pp] = corr(sc2.XData',sup');
legend(ll1,['r=' num2str(round(rr,2)) ', p=' num2str(round(pp,3))],'Location','southeast')

subplot(3,2,6)
[rr,pp] = corr(sc2.XData',bias');
legend(ll2,['r=' num2str(round(rr,2)) ', p=' num2str(round(pp,3))],'Location','northeast')


% Save the figure
set(gcf, 'PaperOrientation', 'landscape', 'PaperUnits', 'normalized');
exportgraphics(gcf, ['figures' filesep 'SI_oriBiasVSscatter_stepByStep.pdf'], 'BackgroundColor', 'none', 'ContentType', 'vector');
