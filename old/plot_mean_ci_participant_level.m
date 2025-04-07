function p = plot_mean_ci_participant_level(x,y,g,rep_size,mv_av,color,varargin)
if length(varargin)<1
    stat_type = 'mean';
else
    stat_type = varargin{1};
end
assert((strcmp(stat_type,'mean') |  strcmp(stat_type,'std')),'stat_type should be either mean or std!')
% Combine data into a table
data = table(x, y, g);

% Use 'groupsummary' to calculate the standard deviation
result = grpstats(data, {'g', 'x'}, stat_type);
if strcmp(stat_type,'std')
    result(result.GroupCount == 1,:)=[];
    result.stat = result.std_y;
else
    result.stat = result.mean_y;
end

groups = unique(result.g);
for gg = 1:length(groups)
    temp = result(strcmp(result.g,groups(gg)),:);
    participant_stat = nan(181,1);
    if any(temp.x<0)
        participant_stat(temp.x+91) = temp.stat;
    else
        participant_stat(temp.x+1) =  temp.stat;
    end
    nan_idx = isnan(participant_stat);
    mv_avg  = (movmean(repmat(participant_stat,rep_size,1),mv_av,'omitnan')');
    mv_avg(nan_idx) = NaN;
    all_participants(gg,:) = mv_avg;
end
ci = 1.96*std(all_participants,'omitnan')./sqrt(size(all_participants,1));
y_plot = mean(all_participants,'omitnan');

% get unique x for the plot
x_plot = unique(x)';x_plot = x_plot(:)';
if rep_size ~= 1
    y_plot   = y_plot(length(x_plot)+1:2*length(x_plot));
    if strcmp(stat_type,'mean') & length(y_plot)>180
        y_plot   = y_plot-y_plot(x_plot==0);
    end
    ci_lower = y_plot - ci(length(x_plot)+1:2*length(x_plot));
    ci_upper = y_plot + ci(length(x_plot)+1:2*length(x_plot));
else
    ci_lower = y_plot - ci;
    ci_upper = y_plot + ci;
end

% Plot confidence intervals as a shaded area
fill([x_plot, fliplr(x_plot)], [ci_upper, fliplr(ci_lower)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'FaceColor', color);
hold on

% Plot the mean values line
p = plot(x_plot, y_plot, 'Color', color, 'LineWidth', 1.5);
xlim([x_plot(1)-5 x_plot(end)+5])
xticks(x_plot(1):30:x_plot(end))
set(gca,'FontSize',15)
end