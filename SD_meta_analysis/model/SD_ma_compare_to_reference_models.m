clc  % Clear command window

% Get the full path of the current script and navigate to the directory
script_path = mfilename('fullpath');
[cur_dir, ~, ~] = fileparts(script_path);
cd(cur_dir)

% Add paths to in-house toolboxes (replace with actual shared folder paths)
addpath(genpath('C:\Users\chare\Google Drive\Work\09_Code\BEIM_toolbox\matlabv\'))
addpath(genpath('/Users/ayberkozkirli/Documents/GitHub/BEIM_toolbox'))

% Change to 'datasets' directory
cd(['..' filesep 'data' filesep 'datasets'])
load('SD_ma_master_table.mat')

%% Simulation Parameters
sigma_enc       = 4:2:30;              % Standard deviations to test
n_iter          = 100;                 % Number of bootstrap iterations
n_trials        = 1000;                % Number of trials per iteration
sup_norm_impr   = nan(n_iter, numel(sigma_enc));
sup_norm_impr_opt = sup_norm_impr;     % optimal cue integration
bins_bay        = nan(n_iter, numel(sigma_enc), 3);
bins_opt        = nan(n_iter, numel(sigma_enc), 3);
hamp_bay        = nan(n_iter, numel(sigma_enc));
hamp_opt        = nan(n_iter, numel(sigma_enc));

% Simulation loop
for iter = 1:n_iter
    fprintf('Iteration = %d of %d\n', iter, n_iter);
    
    for k = 1:numel(sigma_enc)
        % Simulate model for current standard deviation (sigmaEncoding)
        theta_stim  = datasample(0:179, n_trials)';  % Random orientations
        enc_std     = sigma_enc(k);
        dec_std     = datasample(4:30,1);
        trans_prob  = datasample(0.25:.05:.75,1);
        decay_tau   = datasample(.5:.05:1.5,1);
        meas_sigma  = enc_std;%datasample(4:30,1); % 0 to hold fixed
        
        % Run the model for the current parameters
        model_out   = SD_ma_model_bayesian(enc_std, dec_std, trans_prob, decay_tau, theta_stim, [], [], meas_sigma);
        
        % Store SD bias half-amp
        [par,f]     = f_dog(model_out.delta,model_out.error);
        hamp_bay(iter,k) = par(1);
        % Remove the bias from errors
        model_out.error = model_out.error-f(par,model_out.delta);
        
        % Bin model errors based on delta (orientation difference)
        iso_cond    = std(model_out.error(ismember(abs(model_out.delta), 0:10)));
        mid_cond    = std(model_out.error(ismember(abs(model_out.delta), 40:50)));
        ortho_cond  = std(model_out.error(ismember(abs(model_out.delta), 80:90)));
        bins_bay(iter,k,:) = [iso_cond mid_cond ortho_cond];

        % Calculate normalized improvement (Superiority Index)
        sup_norm_impr(iter, k) = (ortho_cond - iso_cond) / iso_cond;

        %------------------------------------------------------------------
        % Predict with simple optimal cue integration
        delta       = sdp_acute_ang(nbk(theta_stim)-theta_stim);
        wp          = 1./(2+(delta/sigma_enc(k)).^2);
        sigma       = sqrt(wp.^2.*sigma_enc(k).^2+(1-wp).^2.*sigma_enc(k)^2);
        errors      = wp.*delta;
         
        % Store SD bias half-amp
        par         = f_dog(delta,errors);
        hamp_opt(iter,k) = par(1);

        % Bin model errors based on delta (orientation difference)
        iso_cond    = mean(sigma(ismember(abs(delta),0:10)));
        mid_cond    = mean(sigma(ismember(abs(delta),40:50)));
        ortho_cond  = mean(sigma(ismember(abs(delta),80:90)));
        bins_opt(iter,k,:) = [iso_cond mid_cond ortho_cond];
        
        % Calculate normalized improvement (Superiority Index)
        sup_norm_impr_opt(iter,k) = (ortho_cond - iso_cond) / iso_cond;
    end
end

%% Calculate Mean and 95% CI of Normalized Improvement
% For Bayesian inference
mean_impr       = mean(sup_norm_impr);        % Mean across bootstrap iterations

% Calculate standard deviation and standard error
std_impr        = std(sup_norm_impr);         % Standard deviation of bootstrap sample
n               = length(sup_norm_impr);      % Number of bootstrap samples
se_impr         = std_impr / sqrt(n);         % Standard error

% Calculate the t-value for 95% CI with (n - 1) degrees of freedom
alpha           = 0.05;
t_crit          = tinv(1 - alpha/2, n - 1);

% Calculate 95% confidence intervals
ci_lower        = mean_impr - t_crit * se_impr;
ci_upper        = mean_impr + t_crit * se_impr;

% For Optimal linear integration
mean_impr_opt   = mean(sup_norm_impr_opt);        % Mean across bootstrap iterations

% Calculate standard deviation and standard error
std_impr        = std(sup_norm_impr_opt);         % Standard deviation of bootstrap sample
n               = length(sup_norm_impr_opt);      % Number of bootstrap samples
se_impr         = std_impr / sqrt(n);         % Standard error

% Calculate the t-value for 95% CI with (n - 1) degrees of freedom
alpha           = 0.05;
t_crit          = tinv(1 - alpha/2, n - 1);

% Calculate 95% confidence intervals
ci_lower_opt    = mean_impr_opt - t_crit * se_impr;
ci_upper_opt    = mean_impr_opt + t_crit * se_impr;

%% Plot Mean and 95% CI of Normalized Improvement

figure;
p               = plot(sigma_enc, mean_impr, 'LineWidth', 1.5, 'Color', 'b');
hold on;
fill([sigma_enc fliplr(sigma_enc)], [ci_upper fliplr(ci_lower)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% optimal
p2              = plot(sigma_enc, mean_impr_opt, 'LineWidth', 1.5, 'Color', 'k');
hold on;
fill([sigma_enc fliplr(sigma_enc)], [ci_upper_opt fliplr(ci_lower_opt)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

xlabel('Encoding Standard Deviation');
ylabel('Normalized Improvement');

% Additional Analysis on `tbl_scatter` Data

% Filter table columns to include only numeric and logical types
var_types       = varfun(@class, tbl, 'OutputFormat', 'cell');
num_log_vars    = ismember(var_types, {'double', 'logical'});
filtered_tbl    = tbl(:, num_log_vars);

% Aggregate and group data by specific columns
tbl_scatter     = grpstats(filtered_tbl, {'bin', 'studynum', 'obsid', 'codenum'}, 'std');
tbl_scatter_count = grpstats(filtered_tbl, {'bin', 'studynum', 'obsid', 'codenum'}, 'numel');

% Filter out rows with insufficient observations
low_count       = tbl_scatter_count.numel_error_ori_deb_sd_deb < 15;
tbl_scatter.numel_error_ori_deb_sd_deb(low_count, :) = nan;

% Calculate group statistics for different bins
group_stats     = grpstats(tbl_scatter, {'codenum', 'bin'}, {'mean', 'meanci'});

% Calculate normalized differences for visual comparison
iso_mean_std    = group_stats.mean_std_error_ori_deb_sd_deb(group_stats.bin == 1);
ortho_mean_std  = group_stats.mean_std_error_ori_deb_sd_deb(group_stats.bin == 3);
norm_diff       = (ortho_mean_std - iso_mean_std) ./ iso_mean_std;

% Calculate overall mean scatter per dataset
group_stats_scat = grpstats(filtered_tbl , {'codenum', 'obsid'}, {'std'});
group_stats_scat = grpstats(group_stats_scat, {'codenum'}, {'mean'});
overall_scat = group_stats_scat.mean_std_error_ori_deb_sd_deb;

% Plot comparison of group statistics with scatter
sc1 = scatter(overall_scat, norm_diff, 'black', 'filled', 'MarkerFaceAlpha', .2, 'MarkerEdgeColor', 'flat');

xlabel('Overall Scatter');
ylabel('Normalized Difference (Ortho - Iso) / Iso');
title('Comparison of datasets scatter to reference models');

xlim([2 32])
hline(0,'k--')

legend([p p2 sc1 ],{'Bayesian','Optimal(linear)','Data'})

%% Plot predictions and CI at the 3 bins
bin_bayesian_matrix     = reshape(bins_bay,n_iter*numel(sigma_enc),3);
bin_bayesian_matrix     = (bin_bayesian_matrix-bin_bayesian_matrix(:,1))./bin_bayesian_matrix(:,1);

bin_optimal_matrix      = reshape(bins_opt,n_iter*numel(sigma_enc),3);
bin_optimal_matrix      = (bin_optimal_matrix-bin_optimal_matrix(:,1))./bin_optimal_matrix(:,1);

colors                  = lines(2);
figure
pl1     = plot_line(bin_bayesian_matrix*100,'x_labels',{'iso','mid','ortho'},'err_plot','dots','color',colors(1,:));
pl2     = plot_line(bin_optimal_matrix*100 ,'x_labels',{'iso','mid','ortho'},'err_plot','dots','color',colors(2,:));
xlim([0 4]);
ylim([-5 35]);
set(gca,'xtick',1:3,'xticklabel',{'iso','mid','ortho'});
format_figure(0,nan,'bin','% change wrt iso');
format_legend([pl1 pl2],{'Bayesian inference','Optimal cue integration'},2,15)