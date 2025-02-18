function bootstrap_params = bootstrap_dog(cfg)

data                = cfg.data;
nBootstrapSamples   = cfg.nBootstrapSamples;

bootstrap_params    = NaN(nBootstrapSamples,2);

for iBootstrap = 1:nBootstrapSamples
   
    bootstrap_sample_idx = randi(length(cfg.subjects),[length(cfg.subjects) 1]);
    bootstrap_sample = cell2mat(data(bootstrap_sample_idx));
    
    funcfg               = [];
    funcfg.fittingsteps  = cfg.fittingsteps;
    funcfg.fixedwidth    = cfg.fixedwidth;
    
    funcfg.data          = bootstrap_sample;
    bootstrap_sample_fit = fit_dog(funcfg);
  
    bootstrap_params(iBootstrap,:) = bootstrap_sample_fit.coeffs;
    
end

save(cfg.savePath,'bootstrap_params')

end

