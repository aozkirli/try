function figure2()
%% Figure 2: bias and scatter with boxplots (all datasets) 
load(['..' filesep '..' filesep 'data' filesep 'datasets' filesep 'SD_ma_master_table.mat'],'tbl')
load('results.mat','tbl_scatter')
FontName        = 'Arial';
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
        dataset_bias = medvalues;
    end
    dataset_scatter = medvalues;

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
savefig(gcf,['figures' filesep 'main_bias_scatter.fig'])
end