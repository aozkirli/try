function [data_clean, data_removed] = outlier_removal(cfg)

data = cfg.data;

data_clean      = cell(size(cfg.subjects,1),1);

for iSubject = 1:size(cfg.subjects,2)

    idx_short   = find(data{iSubject,1}(:,7)==1);
    idx_long    = find(data{iSubject,1}(:,7)==2);
    
    % compute the error for each trial
    trial_errors_short  = mod(data{iSubject,1}(idx_short, 5) - data{iSubject,1}(idx_short, 4) + 90, 180) - 90;
    trial_errors_long   = mod(data{iSubject,1}(idx_long, 5) - data{iSubject,1}(idx_long, 4) + 90, 180) - 90;
    
    mean_error_short    = circ_rad2ang(circ_mean(circ_ang2rad(trial_errors_short)));
    mean_error_long     = circ_rad2ang(circ_mean(circ_ang2rad(trial_errors_long)));
    
    [~,s0_short]    = circ_std(circ_ang2rad(trial_errors_short)); % compute the circular standard deviation
    s0_short        = circ_rad2ang(s0_short);
    
    [~,s0_long]     = circ_std(circ_ang2rad(trial_errors_long));
    s0_long         = circ_rad2ang(s0_long);
        
    cutoff_short    = 3*s0_short;
    cutoff_long     = 3*s0_long;
    
    
    % remove trials with errors larger than 3 sds of error distribution
    data_clean{iSubject,1} = data{iSubject,1};
    
    data_clean{iSubject,1}(idx_short((trial_errors_short > mean_error_short + cutoff_short) | (trial_errors_short < mean_error_short - cutoff_short)),4:7) = nan;
    data_clean{iSubject,1}(idx_long((trial_errors_long > mean_error_long + cutoff_long) | (trial_errors_long < mean_error_long - cutoff_short)),4:7) = nan;
    
    %data_clean{iSubject,1}((trial_errors > mean_error + cutoff) | (trial_errors < mean_error - cutoff),:) = nan;
    
    
    data_removed.short{iSubject,1} = data{iSubject,1}(idx_short((trial_errors_short > mean_error_short + cutoff_short) | (trial_errors_short < mean_error_short - cutoff_short)),:);
    data_removed.long{iSubject,1} = data{iSubject,1}(idx_long((trial_errors_long > mean_error_long + cutoff_long) | (trial_errors_long < mean_error_long - cutoff_short)),:);
    
    %data_removed{iSubject,1} = data{iSubject,1}(abs(trial_errors) > mean_error + cutoff,:);
    
    
    if cfg.plot
        
        figure; hold on;
        plot(data_clean{iSubject,1}(:, 4),data_clean{iSubject,1}(:, 5),'.b');
        plot(data_removed.short{iSubject,1}(:, 4),data_removed.short{iSubject,1}(:, 5),'.r');
        plot(data_removed.long{iSubject,1}(:, 4),data_removed.long{iSubject,1}(:, 5),'.g');
        plot(-90:90,-90:90,'-k')
        xlabel('Stimulus orientation')
        ylabel('Response orientation')

        title(['Responses - Subject ' num2str(cfg.subjects(iSubject))])
        
    end
    
    
end