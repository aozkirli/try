function [rt_stats] = rt_analysis(cfg)

data = cfg.data;

mean_rt = nan(size(cfg.subjects,2),1);

for i = 1:size(cfg.subjects,2)
    
    % remove nans from data
    data{i,1} = data{i,1}(~isnan(data{i,1}(:,1)),:);
    
    mean_rt(i,1) = mean(data{i,1}(:,6));
    
    % Condition on location change wrt to previous stimulus
    rt_prev_same_loc = [];
    rt_prev_diff_loc = [];
    
    for iTrial = 2:length(data{i,1})
        if data{i,1}(iTrial,7) == data{i,1}(iTrial-1,7)
            rt_prev_same_loc = [rt_prev_same_loc; data{i,1}(iTrial,6)];
        elseif data{i,1}(iTrial,7) ~= data{i,1}(iTrial-1,7)
            rt_prev_diff_loc = [rt_prev_diff_loc; data{i,1}(iTrial,6)];
        end
    end
    
    rt_stats.mean_rt_prev_same_loc(i,1) = mean(rt_prev_same_loc);
    rt_stats.mean_rt_prev_diff_loc(i,1) = mean(rt_prev_diff_loc);
    
    % Condition on location change wrt to next stimulus
    rt_next_same_loc = [];
    rt_next_diff_loc = [];
    
    for iTrial = 1:length(data{i,1})-1
        if data{i,1}(iTrial,7) == data{i,1}(iTrial+1,7)
            rt_next_same_loc = [rt_next_same_loc; data{i,1}(iTrial,6)];
        elseif data{i,1}(iTrial,7) ~= data{i,1}(iTrial+1,7)
            rt_next_diff_loc = [rt_next_diff_loc; data{i,1}(iTrial,6)];
        end
    end
    
    rt_stats.mean_rt_next_same_loc(i,1) = mean(rt_next_same_loc);
    rt_stats.mean_rt_next_diff_loc(i,1) = mean(rt_next_diff_loc);
    
    
end

rt_stats.mean_rt     = mean_rt;
rt_stats.group_mean_rt = mean(mean_rt);

rt_stats.group_mean_rt_prev_same_loc = mean(rt_stats.mean_rt_prev_same_loc);
rt_stats.group_mean_rt_prev_diff_loc = mean(rt_stats.mean_rt_prev_diff_loc);

rt_stats.group_mean_rt_next_same_loc = mean(rt_stats.mean_rt_next_same_loc);
rt_stats.group_mean_rt_next_diff_loc = mean(rt_stats.mean_rt_next_diff_loc);

end