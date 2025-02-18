function group_permutations(cfg)

data = cfg.data;

% matrix for storing the model parameters for each permutation
permutation_distribution = zeros(cfg.nperms,2);

for i = 1:cfg.nperms
    
    disp(['Permutation ' num2str(i) ' of ' num2str(cfg.nperms) '...'])
    
    % randomly flip signs for each participant's responses
    A = [-1 1];
    tmp_permsign = randi([1 2],size(cfg.subjects,2),1);
    permsign = A(tmp_permsign);
    
    perm_data = data;
    
    for iSubject = 1:length(data)   
        perm_data{iSubject}(:,4) = permsign(iSubject) * perm_data{iSubject}(:,4);      
    end
    
    perm_data = cell2mat(perm_data);

    funcfg                              = [];
    funcfg.data                         = perm_data(:,[3 4]);
    funcfg.fittingsteps                 = cfg.fittingsteps;
    funcfg.fixedwidth                   = cfg.fixedwidth;
    if cfg.fixedwidth == true
        funcfg.width                    = cfg.width_same_loc;
    end
    permutation_fit                     = fit_dog(funcfg);
       
    permutation_distribution(i,:)       = permutation_fit.coeffs;
   
end

save(cfg.savePath, 'permutation_distribution')


end

