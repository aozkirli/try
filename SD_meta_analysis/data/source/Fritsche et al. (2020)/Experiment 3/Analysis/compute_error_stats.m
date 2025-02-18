function errorstats = compute_error_stats(cfg)

data = cfg.data;

for iSubject = 1:length(data)
   
    error = mod(data{iSubject}(:,5) - data{iSubject}(:,4) + 90, 180) - 90;
    sd_error(iSubject,1) = nanstd(error);
    
end

group_sd_error = mean(sd_error);

errorstats.sd_error = sd_error;
errorstats.group_sd_error = group_sd_error;

end

