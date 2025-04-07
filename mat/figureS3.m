function figureS3()
load(['..' filesep '..' filesep 'data' filesep 'datasets' filesep 'SD_ma_master_table.mat'],'tbl')
FontName = 'Arial';
%% Figure S3.1: Show stimulus-specific bias cleaning effects on error scatter (as done in Fritsche et al. 2020)
bin_size_Fristche = 30;
figure('Units','normalized','position',[0 0 1 1]);
tmp = tbl(abs(tbl.delta)<=90,:);
tmp.theta= round(tmp.theta);
tmp.delta = round(tmp.delta);
datasets = unique(tmp.codenum);

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
        scatt_subj_iqr(:,oo)  = movmean(scatt,bin_size_Fristche,'omitnan')-mean(movmean(scatt,bin_size_Fristche,'omitnan'));

        tto.error_ori_deb(tto.delta<0) = -tto.error_ori_deb(tto.delta<0);
        scatt = nan(1,91);
        scatt(ismember([0:90],unique(abs(tto.delta)))) = grpstats(tto.error_ori_deb,abs(tto.delta),'std');scatt(scatt==0)=NaN;
        scatt_subj_ori_deb(:,oo)  = movmean(scatt,bin_size_Fristche,'omitnan')-mean(movmean(scatt,bin_size_Fristche,'omitnan'));
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
savefig(gcf,['figures' filesep 'SI_allDatasetsScatter.fig'])

%% Figure S3.2 and S3.3: Show how stimulus-specific introduces artifacts to error scatter
ori_bias_removed = {'with_oriBias','without_oriBias'};
whichVariable = {'error_iqr_norm','error_ori_deb_norm'};
rep_size = 3;
bin_size = 31;
n_bins = 10;
intersect_thetas = [0:181];
pastel_cmap = blue2pink(n_bins);

% filter the data and calculate the effect size of stim-specific bias for each theta
tmp    = tbl(abs(tbl.delta)<=90 & tbl.theta<=180 & ~isnan(tbl.delta),:);
dd = arrayfun(@(i) nanmean(tmp.error_iqr_norm(tmp.theta == i)) / ...
    nanstd(tmp.error_iqr_norm(tmp.theta == i)), 0:180);
[abs_dd_sorted,idx]=sort(abs(dd));

for removed = 1:2
    error_variable = whichVariable{removed};

    % Create a new figure
    figure('Units','normalized','position',[0 0 1 1]);
    % Define positions for 2x2 on the left-hand side
    left_positions = [0 .065 0 0]+[
        0.1, 0.5, 0.2, 0.35;  % Top-left
        0.1, 0.03, 0.2, 0.35;  % Bottom-left
        0.4, 0.5, 0.2, 0.35;  % Top-right
        0.4, 0.03, 0.2, 0.35   % Bottom-right
        ];

    superiority = nan(n_bins,1);
    sd_peak     = nan(n_bins,1);
    cohenD_mean = nan(n_bins,1);

    letters = {'A)', 'B)', 'C)', 'D)'};
    % Define positions for the 2x2 grid (left)
    for i = 1:4
        subplot('Position', left_positions(i,:));
        pos = get(gca, 'Position');
        annotation('textbox', [pos(1)-0.06, pos(2) + pos(4), 0.05, 0.05], ...
            'String', letters{i}, 'HorizontalAlignment', 'center', ...
            'FontName', FontName, 'FontSize', 22, 'EdgeColor', 'none', 'Units', 'normalized');
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
                        pos = get(gca, 'Position');
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
                        pos = get(gca, 'Position');
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
                        pos = get(gca, 'Position');
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
        0.675, 0.725, 0.12, 0.22;  % Top-right
        0.675, 0.4, 0.12, 0.22;  % Middle-right
        0.675, 0.075, 0.12, 0.22   % Bottom-right
        ];

    letters = {'E)', 'F)', 'G)'};
    % % Define positions for the 3x1 grid (right)
    for i = 1:3
        % Adjust position manually for the 3x1 grid
        subplot('Position', right_positions(i,:));
        pos = get(gca, 'Position');
        annotation('textbox', [pos(1)-0.06, pos(2) + pos(4), 0.05, 0.05], ...
            'String', letters{i}, 'HorizontalAlignment', 'center', ...
            'FontName', FontName, 'FontSize', 22, 'EdgeColor', 'none', 'Units', 'normalized');
                    
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
    exportgraphics(gcf, ['figures' filesep 'SI_confound_ ' ori_bias_removed{removed} '.pdf'], 'BackgroundColor', 'none', 'ContentType', 'vector');
    savefig(gcf,['figures' filesep 'SI_confound_ ' ori_bias_removed{removed} '.fig'])

end

%% Figure S3.4: Mean stim-specific bias amplitude vs. mean superiority across datasets
tmp = tbl(abs(tbl.delta)<=90,:);
datasets = unique(tmp.code);
colors          = linspecer(length(datasets));
ss_bias = cell(length(datasets),1);

ssfreethr = [10:25 100:115; 65:80 155:170];
figure('Units','normalized','position',[0.9 0.3 .3 .4]);
for i = 1:length(datasets)
    tbl_i               = tmp(tmp.codenum==i,:);
    observers = unique(tbl_i.obs);
    ssbias = nan(length(observers),1);
    sup    = nan(length(observers),1);
    for o = 1:length(observers)
        tbl_io               = tbl_i(tbl_i.obs==o,:);
        ssbias(o) = mean(tbl_io.error_iqr(ismember(tbl_io.theta,ssfreethr(1,:))),'omitnan')-...
            mean(tbl_io.error_iqr(ismember(tbl_io.theta,ssfreethr(2,:))),'omitnan');
        sup(o) = std(tbl_io.error_iqr(tbl_io.bin==3),'omitnan')-std(tbl_io.error_iqr(tbl_io.bin==1),'omitnan');
    end
    idx = ~isnan(sup) & ~isnan(ssbias);
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
tbl_datasets_mean.ss_bias = cellfun(@mean, ss_bias);

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


set(gcf,'PaperOrientation','landscape')
set(gcf, 'PaperUnits', 'normalized');
exportgraphics(gcf, ['figures' filesep 'SI_oriBias_scatter_dataset.pdf'], 'BackgroundColor', 'none', 'ContentType', 'vector');
savefig(gcf,['figures' filesep 'SI_oriBias_scatter_dataset.fig'])

end
