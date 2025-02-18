function p = plot_mean_ci_bootstrap(x, y, rep_size, mv_av, color, varargin)
if length(varargin) < 1
    stat_type = 'mean';
    n_bootstrap = [];
else
    stat_type = varargin{1};
    n_bootstrap = varargin{2};
end

% Get unique x values
x_unique = unique(x);

if isempty(n_bootstrap)
    y_plot = nan(181,1);
    n      = nan(181,1);
    if any(x_unique<0)
        y_plot(x_unique+91) = grpstats(y,x,stat_type);
        n(x_unique+91) = grpstats(y,x,'numel');
        x_plot = -90:90; x_plot = x_plot(:)';
    else
        y_plot(x_unique+1)  = grpstats(y,x,stat_type);
        n(x_unique+1) = grpstats(y,x,'numel');
        x_plot = 0:180; x_plot = x_plot(:)';
    end
    y_plot = movmean(repmat(y_plot,rep_size,1),mv_av,'omitnan')';
    n      = movsum(repmat(n,rep_size,1),mv_av,'omitnan')';

    alpha = 0.05;
    % Chi-squared values for CI
    chi2_upper = chi2inv(1 - alpha / 2, n - 1);
    chi2_lower = chi2inv(alpha / 2, n - 1);
    
    % Confidence intervals for standard deviation
    ci_lower = y_plot .* sqrt((n - 1)  ./ chi2_upper);
    ci_upper = y_plot .* sqrt((n - 1)  ./ chi2_lower);

    % get unique x for the plot
    % x_plot = unique(x)';x_plot = x_plot(:)';
    if rep_size ~= 1
        y_plot   = y_plot(length(x_plot)+1:2*length(x_plot));
        if strcmp(stat_type,'mean') & length(y_plot)>180
            y_plot   = y_plot-y_plot(x_plot==0);
        end
        ci_lower = ci_lower(length(x_plot)+1:2*length(x_plot));
        ci_upper = ci_upper(length(x_plot)+1:2*length(x_plot));
    else
        ci_lower = y_plot - ci;
        ci_upper = y_plot + ci;
    end
else
    % Parameters for bootstrapping
    ci_percentile = 0.975; % Upper percentile for 95% CI (two-tailed)

    % Generate bootstrap resamples
    bootstrap_samples = zeros(n_bootstrap, length(unique(x)));
    parfor b = 1:n_bootstrap
        [resampled_y,resampling_idx] = datasample(y, length(y));
        resampled_x = x(resampling_idx);
        if strcmp(stat_type,'mean')
            bootstrap_samples(b,:) = grpstats(resampled_y,resampled_x);
        elseif strcmp(stat_type,'std')
            bootstrap_samples(b,:) = grpstats(resampled_y,resampled_x,'std');
        end
    end

    % Smooth the bootstrap samples
    bootstrap_samples_mv_av = movmean(repmat(bootstrap_samples, 1, rep_size)', mv_av)';
    % Compute mean and CI from bootstrap samples
    y_boot_mean = mean(bootstrap_samples_mv_av,'omitnan');
    y_boot_ci_lower = prctile(bootstrap_samples_mv_av, (1 - ci_percentile) * 100);
    y_boot_ci_upper = prctile(bootstrap_samples_mv_av, ci_percentile * 100);

    % Adjust for rep_size if needed
    if rep_size ~= 1
        y_plot = y_boot_mean(length(x_unique) + 1:2 * length(x_unique));
        ci_lower = y_boot_ci_lower(length(x_unique) + 1:2 * length(x_unique));
        ci_upper = y_boot_ci_upper(length(x_unique) + 1:2 * length(x_unique));
    end
end
% Plot confidence intervals as a shaded area
fill([x_plot fliplr(x_plot)], [ci_upper fliplr(ci_lower)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'FaceColor', color);
hold on;

% Plot the mean values line
p = plot(x_plot, y_plot, 'Color', color, 'LineWidth', 1.5);
xlim([x_unique(1) - 5, x_unique(end) + 5]);
xticks(x_unique(1):30:x_unique(end));
set(gca, 'FontSize', 15);
