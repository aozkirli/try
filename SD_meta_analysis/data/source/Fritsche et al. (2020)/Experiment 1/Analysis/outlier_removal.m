function [data_clean, data_removed] = outlier_removal(cfg)

data = cfg.data;

data_clean      = cell(size(cfg.subjects,1),1);
data_removed    = cell(size(cfg.subjects,1),1);

for i = 1:size(cfg.subjects,2)

    % compute the error for each trial
    trial_errors    = mod(data{i,1}(1:end, 5) - data{i,1}(1:end, 4) + 90, 180) - 90;
    
    mean_error      = circ_rad2ang(circ_mean(circ_ang2rad(trial_errors)));
    
    %{
    s = sum(sin(circ_ang2rad(trial_errors)));
    c = sum(cos(circ_ang2rad(trial_errors)));
    r_bar = sqrt(s^2 + c^2) / size(trial_errors,1);
    delta_hat = (1-sum(cos(2*circ_ang2rad(trial_errors - mean_error))/size(trial_errors,1))) / (2*r_bar^2);
    constant = 6;
    
    ci = circ_ang2rad(mean_error) + [-1 1] * constant * delta_hat;
    
    figure;
    plot(cos(circ_ang2rad(trial_errors)),sin(circ_ang2rad(trial_errors)),'.');
    xlim([-1 1])
    ylim([-1 1])
    line([0 cos(ci(1))],[0 sin(ci(1))])
    line([0 cos(ci(2))],[0 sin(ci(2))])
    %}
    
    [~,s0]  = circ_std(circ_ang2rad(trial_errors)); % compute the circular standard deviation
    s0      = circ_rad2ang(s0);
    cutoff  = 3*s0;
    
    % remove trials with errors larger than 3 sds of error distribution
    data_clean{i,1} = data{i,1};
    data_clean{i,1}(abs(trial_errors) > mean_error + cutoff,4:end) = nan;
    
    data_removed{i,1} = data{i,1}(abs(trial_errors) > mean_error + cutoff,:);
    
end