function [data_clean, data_removed] = outlier_removal(cfg)

data = cfg.data;

data_clean      = cell(size(cfg.subjects,1),1);
data_removed    = cell(size(cfg.subjects,1),1);

for iSubject = 1:size(cfg.subjects,2)

    % compute the error for each trial
    trial_errors    = mod(data{iSubject,1}(:, 5) - data{iSubject,1}(:, 4) + 90, 180) - 90;
    
    mean_error      = circ_rad2ang(circ_mean(circ_ang2rad(trial_errors)));
    
    [~,s0]  = circ_std(circ_ang2rad(trial_errors)); % compute the circular standard deviation
    s0      = circ_rad2ang(s0);
    
    cutoff  = 3*s0;
    
    % remove trials with errors larger than 3 sds of error distribution
    data_clean{iSubject,1} = data{iSubject,1};
    data_clean{iSubject,1}((trial_errors > mean_error + cutoff) | (trial_errors < mean_error - cutoff),4:end) = nan;
    
    data_removed{iSubject,1} = data{iSubject,1}(abs(trial_errors) > mean_error + cutoff,:);
    
    
    if cfg.plot
        
        figure; hold on;
        plot(data_clean{iSubject,1}(:, 4),data_clean{iSubject,1}(:, 5),'.b');
        plot(data_removed{iSubject,1}(:, 4),data_removed{iSubject,1}(:, 5),'.r');
        plot(-90:90,-90:90,'-k')
        xlabel('Stimulus orientation')
        ylabel('Response orientation')

        title(['Responses - Subject ' num2str(cfg.subjects(iSubject))])
        
    end
    
    
end