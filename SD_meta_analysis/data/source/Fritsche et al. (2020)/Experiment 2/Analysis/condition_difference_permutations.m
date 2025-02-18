function condition_difference_permutations(cfg)

data_cond1 = cfg.data_cond1;
data_cond2 = cfg.data_cond2;

% matrix for storing the model parameters for each permutation
permutation_distribution_cond1 = zeros(cfg.nperms,2);
permutation_distribution_cond2 = zeros(cfg.nperms,2);

for i = 1:cfg.nperms
    
    disp(['Permutation ' num2str(i) ' of ' num2str(cfg.nperms) '...'])
    
    %% Make a random partition
    rnd_partition = randi([1 2],size(cfg.subjects,2),1);
 
    perm_data_cond1 = cell(size(cfg.subjects,2),1);
    perm_data_cond2 = cell(size(cfg.subjects,2),1);
    
    for iSubject = 1:size(cfg.subjects,2)
       
        if rnd_partition(iSubject) == 1
            perm_data_cond1{iSubject} = data_cond1{iSubject};
            perm_data_cond2{iSubject} = data_cond2{iSubject};
        elseif rnd_partition(iSubject) == 2
            perm_data_cond1{iSubject} = data_cond2{iSubject};
            perm_data_cond2{iSubject} = data_cond1{iSubject};
        end
        
    end
    
    perm_data_cond1 = cell2mat(perm_data_cond1);
    perm_data_cond2 = cell2mat(perm_data_cond2);
    
    %% Condition 1  
    funcfg                              = [];
    funcfg.data                         = perm_data_cond1(:,[3 4]);
    funcfg.fittingsteps                 = cfg.fittingsteps;
    funcfg.fixedwidth                   = cfg.fixedwidth;
    if cfg.fixedwidth == true
        funcfg.width                    = cfg.width_cond1;
    end
    permutation_fit                     = fit_dog(funcfg);
       
    permutation_distribution_cond1(i,:) = permutation_fit.coeffs;
    
    %% Condition 2
    funcfg                              = [];
    funcfg.data                         = perm_data_cond2(:,[3 4]);
    funcfg.fittingsteps                 = cfg.fittingsteps;
    funcfg.fixedwidth                   = cfg.fixedwidth;
    if cfg.fixedwidth == true
        funcfg.width                    = cfg.width_cond2;
    end
    permutation_fit                     = fit_dog(funcfg);
       
    permutation_distribution_cond2(i,:) = permutation_fit.coeffs;
     
end

permutation_distribution_amplitude_difference   = permutation_distribution_cond1(:,1) - permutation_distribution_cond2(:,1);
permutation_distribution_width_difference       = permutation_distribution_cond1(:,2) - permutation_distribution_cond2(:,2);

save(cfg.savePath, 'permutation_distribution_cond1','permutation_distribution_cond2','permutation_distribution_amplitude_difference','permutation_distribution_width_difference')

end

