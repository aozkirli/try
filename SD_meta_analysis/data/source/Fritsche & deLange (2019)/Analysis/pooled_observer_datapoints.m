function [mean_datapoints,sem_datapoints,mean_smoothed_datapoints,sem_smoothed_datapoints] = pooled_observer_datapoints(cfg)

subjects = cfg.subjects;
sd_data_ori = cfg.sd_data_ori;
sd_data_size = cfg.sd_data_size;

ifc_ref_stim_relative_orientation = -80:10:90;

%% Orientation condition

for iSubject = 1:length(subjects)
    
    % unsmoothed
    for iStimLevel = 1:length(ifc_ref_stim_relative_orientation)
        mean_subject_datapoints{iSubject,1}(iStimLevel) = ...
            mean(sd_data_ori{iSubject}(sd_data_ori{iSubject}(:,1) == ifc_ref_stim_relative_orientation(iStimLevel),2));
    end
    
    % smoothed
    for iStimLevel = 2:length(ifc_ref_stim_relative_orientation)-1
        mean_smoothed_subject_datapoints{iSubject,1}(iStimLevel) = ...
            (mean_subject_datapoints{iSubject,1}(iStimLevel-1)+ ...
             mean_subject_datapoints{iSubject,1}(iStimLevel)+...
             mean_subject_datapoints{iSubject,1}(iStimLevel+1))/3;
    end
    
    mean_smoothed_subject_datapoints{iSubject,1}(1) = ...
        (mean_subject_datapoints{iSubject,1}(length(ifc_ref_stim_relative_orientation)) + ...
         mean_subject_datapoints{iSubject,1}(1) +  ...
         mean_subject_datapoints{iSubject,1}(2))/3;
    
    mean_smoothed_subject_datapoints{iSubject,1}(length(ifc_ref_stim_relative_orientation)) = ...
        (mean_subject_datapoints{iSubject,1}(length(ifc_ref_stim_relative_orientation)-1) + ...
         mean_subject_datapoints{iSubject,1}(length(ifc_ref_stim_relative_orientation)) + ...
         mean_subject_datapoints{iSubject,1}(1))/3;
    
end

mean_subject_datapoints = cell2mat(mean_subject_datapoints);
mean_smoothed_subject_datapoints = cell2mat(mean_smoothed_subject_datapoints);

mean_datapoints(1,:) = mean(mean_subject_datapoints);
sem_datapoints(1,:) = std(mean_subject_datapoints,1,1) / sqrt(length(subjects));

mean_smoothed_datapoints(1,:) = mean(mean_smoothed_subject_datapoints);
sem_smoothed_datapoints(1,:) = std(mean_smoothed_subject_datapoints,1,1) / sqrt(length(subjects));

clear mean_subject_datapoints mean_smoothed_subject_datapoints

%% Size condition

for iSubject = 1:length(subjects)
    
    % unsmoothed
    for iStimLevel = 1:length(ifc_ref_stim_relative_orientation)
        mean_subject_datapoints{iSubject,1}(iStimLevel) = ...
            mean(sd_data_size{iSubject}(sd_data_size{iSubject}(:,1) == ifc_ref_stim_relative_orientation(iStimLevel),2));
    end
    
    % smoothed
    for iStimLevel = 2:length(ifc_ref_stim_relative_orientation)-1
        mean_smoothed_subject_datapoints{iSubject,1}(iStimLevel) = ...
            (mean_subject_datapoints{iSubject,1}(iStimLevel-1)+ ...
             mean_subject_datapoints{iSubject,1}(iStimLevel)+...
             mean_subject_datapoints{iSubject,1}(iStimLevel+1))/3;
    end
    
    mean_smoothed_subject_datapoints{iSubject,1}(1) = ...
        (mean_subject_datapoints{iSubject,1}(length(ifc_ref_stim_relative_orientation)) + ...
         mean_subject_datapoints{iSubject,1}(1) +  ...
         mean_subject_datapoints{iSubject,1}(2))/3;
    
    mean_smoothed_subject_datapoints{iSubject,1}(length(ifc_ref_stim_relative_orientation)) = ...
        (mean_subject_datapoints{iSubject,1}(length(ifc_ref_stim_relative_orientation)-1) + ...
         mean_subject_datapoints{iSubject,1}(length(ifc_ref_stim_relative_orientation)) + ...
         mean_subject_datapoints{iSubject,1}(1))/3;
    
end

mean_subject_datapoints = cell2mat(mean_subject_datapoints);
mean_smoothed_subject_datapoints = cell2mat(mean_smoothed_subject_datapoints);

mean_datapoints(2,:) = mean(mean_subject_datapoints);
sem_datapoints(2,:) = std(mean_subject_datapoints,1,1) / sqrt(length(subjects));

mean_smoothed_datapoints(2,:) = mean(mean_smoothed_subject_datapoints);
sem_smoothed_datapoints(2,:) = std(mean_smoothed_subject_datapoints,1,1) / sqrt(length(subjects));

end

