function figure3B()
load('results.mat')
FontName = 'Arial';
%% Figure 3B: Mean bias moving average and mean error scatter + polynomial fit
x = 0:90;
figure('Units','normalized','position',[.05 .4 .6 .4]);
% prepare bias for plotting
mvav_bias_aggregated_curve = mvav_bias_aggregated(:)';
% prepare scatter for plotting
mvav_scatter_aggregated_curve = mvav_scatter_aggregated(:)';
aggregated_fit_curve = best_fit_all(:)';

subplot(121);hold off;
hold on;
p=plot_mean_ci(repmat(0:90,1,length(aggregated_fit_curve)/91),mvav_bias_aggregated_curve,1,1,'k','mean',[]);
p.LineWidth = 6;
p.LineStyle = ':';

subplot(122)
hold on;
p=plot_mean_ci(repmat(0:90,1,length(aggregated_fit_curve)/91),mvav_scatter_aggregated_curve,1,1,'k','mean',[]);
p.LineWidth = 6;
p.LineStyle = ':';
p=plot_mean_ci(repmat(0:90,1,length(aggregated_fit_curve)/91),aggregated_fit_curve,1,1,'b','mean',[]);
p.LineWidth = 6;


subplot(121)
% Configure plot
xlabel('|Δ| (°)');
ylabel('Bias (°)');
title(['Data']);
xticks([0 45 90]);
set(gca, 'FontSize', 20, 'LineWidth', 1.9, 'GridLineWidth', 1.9, 'FontName', FontName);
xticklabels({'iso' 'mid' 'ortho'});

subplot(122)
% Configure plot
xlabel('|Δ| (°)');
ylabel('Error Scatter (°)');
title(['Data vs. Polynomial Fit']);
xticks([0 45 90]);
set(gca, 'FontSize', 20, 'LineWidth', 1.9, 'GridLineWidth', 1.9, 'FontName', FontName);
xticklabels({'iso' 'mid' 'ortho'});
set(gcf,'PaperOrientation','landscape')
set(gcf, 'PaperUnits', 'normalized');
exportgraphics(gcf,['figures' filesep 'main_results_polyfit.pdf'],'BackgroundColor','none')
savefig(gcf, ['figures' filesep 'main_results_polyfit.fig']);

end