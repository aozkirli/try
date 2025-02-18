function group_fit = group_fitting(cfg)

group_fit = fit_dog(cfg);
save(cfg.savePath,'group_fit');

end

