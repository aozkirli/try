addpath /home/common/matlab/fieldtrip/qsub/
addpath /home/predatt/matfri/Toolboxes/CircStat2012a

clear all; close all;

%% Main analysis of the estimation task experiment

%% Load data
subjects = [1:9,11:25];

cfg = [];
cfg.subjects = subjects;

[data, log] = load_data(cfg);

%% Outlier removal
cfg                     = [];
cfg.data                = data;
cfg.subjects            = subjects;
cfg.plot                = false;

[data_outlier_rm, data_removed]    = outlier_removal(cfg);

%% Compute error stats
cfg = [];
cfg.data = data_outlier_rm;

errorstats = compute_error_stats(cfg);


%% RT analysis
cfg             = [];
cfg.subjects    = subjects;
cfg.data        = data;

[rt_stats] = rt_analysis(cfg);


%% Sort SD data
cfg             = [];
cfg.subjects    = subjects;
cfg.data        = data_outlier_rm;
cfg.nBack       = 10;
cfg.demean      = true; % demean response errors to remove overall bias

[sd_data] = sort_data(cfg);


%% Compute moving averages
cfg             = [];
cfg.data        = sd_data;
cfg.bin_width   = 20; % possible to specify multiple window sizes
        
[moving_averages] = compute_moving_averages(cfg);


%% ===================================================================== %%
%% Fitting to  Dog Models to group data (pooled datapoints)
compiled_group_fitting  = qsubcompile('group_fitting','toolbox', {'curvefit'});

time            = 1 * 3600; % hours to seconds
space           = 1 * 1024^3; % Gb to bytes

cfg                 = [];
cfg.fittingsteps    = 100;
cfg.fixedwidth      = false;

for iBack = 1:10
    
    data_all            = cell2mat(sd_data.all{iBack});
    data_same_loc       = cell2mat(sd_data.same_loc{iBack});
    data_diff_loc       = cell2mat(sd_data.diff_loc{iBack});
    
    % all trials
    cfg.data        = data_all(:,[3 4]);
    cfg.savePath    = ['Matlab datafiles/group fits/group_fits_all_' num2str(iBack) '-back.mat'];
    jobidarray{1}   = qsubfeval(compiled_group_fitting,cfg,'timreq',time,'memreq',space);
       
    % same location
    cfg.data        = data_same_loc(:,[3 4]);
    cfg.savePath    = ['Matlab datafiles/group fits/group_fits_same_loc_' num2str(iBack) '-back.mat'];
    jobidarray{2}   = qsubfeval(compiled_group_fitting,cfg,'timreq',time,'memreq',space);
    
    % different location
    cfg.data        = data_diff_loc(:,[3 4]);
    cfg.savePath    = ['Matlab datafiles/group fits/group_fits_diff_loc_' num2str(iBack) '-back.mat'];
    jobidarray{3}   = qsubfeval(compiled_group_fitting,cfg,'timreq',time,'memreq',space);
       
end



%% ===================================================================== %%
%% ================= Bootstrap DoG parameter estimates ================= %%

% precompile permutation function for cluster execution
compiled_bootstrap_dog  = qsubcompile('bootstrap_dog','toolbox', {'curvefit'});
time = 4 * 3600; % hours to seconds
space = 1 * 1024^3; % Gb to bytes

cfg                     = [];
cfg.subjects            = subjects;
cfg.fittingsteps        = 1;
cfg.nBootstrapSamples   = 1000;
cfg.fixedwidth          = false;

for iBack = 1:10
    
    % all trials
    cfg.data                = sd_data.all{iBack};
    cfg.savePath            = ['Matlab datafiles/pooled_observer_bootstrap_all_' num2str(iBack) '-back.mat'];
    jobidarray{1}           = qsubfeval(compiled_bootstrap_dog,cfg,'timreq',time,'memreq',space);  
    
    % same location
    cfg.data                = sd_data.same_loc{iBack};
    cfg.savePath            = ['Matlab datafiles/pooled_observer_bootstrap_same_loc_' num2str(iBack) '-back.mat'];
    jobidarray{1}           = qsubfeval(compiled_bootstrap_dog,cfg,'timreq',time,'memreq',space);    
    
    % different location
    cfg.data                = sd_data.diff_loc{iBack};
    cfg.savePath            = ['Matlab datafiles/pooled_observer_bootstrap_diff_loc_' num2str(iBack) '-back.mat'];
    jobidarray{2}           = qsubfeval(compiled_bootstrap_dog,cfg,'timreq',time,'memreq',space);
    
end


%% ===================================================================== %%
%% Permutation test for testing SD against zero

% precompile permutation function for cluster execution
compiled_group_permutations  = qsubcompile('group_permutations','toolbox', {'curvefit'});
time = 36 * 3600; % hours to seconds
space = 1 * 1024^3; % Gb to bytes


cfg                 = [];
cfg.subjects        = subjects;
cfg.nperms          = 10000;
cfg.fittingsteps    = 1;
cfg.fixedwidth      = false;

for iBack = 1:10
    
    % all trials
    cfg.data        = sd_data.all{iBack};
    cfg.savePath    = ['Matlab datafiles/group fits/all_permutation_fits_' num2str(iBack) '-back.mat'];
    jobidarray{1}   = qsubfeval(compiled_group_permutations,cfg,'timreq',time,'memreq',space);   
    
    % same location
    cfg.data        = sd_data.same_loc{iBack};
    cfg.savePath    = ['Matlab datafiles/group fits/same_loc_permutation_fits_' num2str(iBack) '-back.mat'];
    jobidarray{2}   = qsubfeval(compiled_group_permutations,cfg,'timreq',time,'memreq',space);   
    
    % different location
    cfg.data        = sd_data.diff_loc{iBack};
    cfg.savePath    = ['Matlab datafiles/group fits/diff_loc_permutation_fits_' num2str(iBack) '-back.mat'];
    jobidarray{3}   = qsubfeval(compiled_group_permutations,cfg,'timreq',time,'memreq',space);
    
end


%% ===================================================================== %%
%% Results of permutation test of SD against zero

for iBack = 1:10
   
    %% load fitting parameters
    
    load(['Matlab datafiles/group fits/group_fits_all_' num2str(iBack) '-back.mat']);
    params_all(iBack,:) = group_fit.coeffs;
    
    %% load permutation distribution
    load(['Matlab datafiles/group fits/all_permutation_fits_' num2str(iBack) '-back.mat']);
    permutation_distribution_all = permutation_distribution;
    
    %% compute p-values
    if params_all(iBack,1) <= 0
        p_val_all(iBack,1)      = sum(permutation_distribution_all(:,1) <= params_all(iBack,1)) / length(permutation_distribution_all);
    else
        p_val_all(iBack,1)      = sum(permutation_distribution_all(:,1) >= params_all(iBack,1)) / length(permutation_distribution_all);
    end
    

end

p_val_all = p_val_all * 2;


%% ===================================================================== %%
%% Permutation test for testing difference of SD between location switches

% precompile permutation function for cluster execution
compiled_condition_difference_permutations  = qsubcompile('condition_difference_permutations','toolbox', {'curvefit'});


cfg                 = [];
cfg.subjects        = subjects;

for iBack = 1:10
    
    cfg.data_cond1   = sd_data.same_loc{iBack};
    cfg.data_cond2   = sd_data.diff_loc{iBack};
    
    cfg.fixedwidth      = false;
    
    cfg.nperms          = 10000;
    cfg.fittingsteps    = 1;
    cfg.savePath        = ['Matlab datafiles/group fits/permutation_fits_location_diff_' num2str(iBack) '-back.mat'];
    
    time                = 36 * 3600; % hours to seconds
    space               = 1 * 1024^3; % Gb to bytes
    
    jobidarray{iBack}    = qsubfeval(compiled_condition_difference_permutations,cfg,'timreq',time,'memreq',space);   
    
end


%% ===================================================================== %%
%% Results of permutation test of differences between location switches

for iBack = 1:10
   
    %% load fitting parameters
    load(['Matlab datafiles/group fits/group_fits_same_loc_' num2str(iBack) '-back.mat']);
    params_same_loc(iBack,:) = group_fit.coeffs;
    
    load(['Matlab datafiles/group fits/group_fits_diff_loc_' num2str(iBack) '-back.mat']);
    params_diff_loc(iBack,:) = group_fit.coeffs;
    
    amp_diff(iBack,1) = params_same_loc(iBack,1) - params_diff_loc(iBack,1);
    
    
    %% load permutation distribution    
    load(['Matlab datafiles/group fits/permutation_fits_location_diff_' num2str(iBack) '-back.mat']);

    %% compute p-values
    p_val_location_diff(iBack,1) = sum(permutation_distribution_amplitude_difference <= amp_diff(iBack,1)) / length(permutation_distribution_amplitude_difference); 
    
    
end

%% ========================================================================
%% ========================= Bar plots ====================================

%% compute error bars
for iBack = 1:10
   
    load(['Matlab datafiles/pooled_observer_bootstrap_all_' num2str(iBack) '-back.mat']);
    sem_all(iBack,:) = std(bootstrap_params(:,1));
    
end


figure; hold on;
bar(1:10,params_all(:,1)');
errorbar(1:10,params_all(:,1)',sem_all,'.k','LineWidth',2,'Capsize',0)
xlabel('n-back')
ylabel('Bias amplitude (deg)')
ylim([-1.5 2])




%% ===================================================================== %%
%% === Permutation test for location differences of avg n-back trials == %%

%% pool the n-back data into subject cells
for iSubject = 1:length(subjects)
    
    pooled_sddata_same_loc_49{iSubject,1} = [sd_data.same_loc{4}{iSubject};...
        sd_data.same_loc{5}{iSubject};...
        sd_data.same_loc{6}{iSubject};...
        sd_data.same_loc{7}{iSubject};...
        sd_data.same_loc{8}{iSubject};...
        sd_data.same_loc{9}{iSubject}];
    
    
    pooled_sddata_diff_loc_49{iSubject,1} = [sd_data.diff_loc{4}{iSubject};...
        sd_data.diff_loc{5}{iSubject};...
        sd_data.diff_loc{6}{iSubject};...
        sd_data.diff_loc{7}{iSubject};...
        sd_data.diff_loc{8}{iSubject};...
        sd_data.diff_loc{9}{iSubject}];
    
    % also pool the moving averages
    pooled_moving_average_same_loc_49{iSubject,1} = mean([moving_averages.same_loc{4}{iSubject};...
        moving_averages.same_loc{5}{iSubject};...
        moving_averages.same_loc{6}{iSubject};...
        moving_averages.same_loc{7}{iSubject};...
        moving_averages.same_loc{8}{iSubject};...
        moving_averages.same_loc{9}{iSubject}],1);
    
    pooled_moving_average_diff_loc_49{iSubject,1} = mean([moving_averages.diff_loc{4}{iSubject};...
        moving_averages.diff_loc{5}{iSubject};...
        moving_averages.diff_loc{6}{iSubject};...
        moving_averages.diff_loc{7}{iSubject};...
        moving_averages.diff_loc{8}{iSubject};...
        moving_averages.diff_loc{9}{iSubject}],1);
        
end


%% ======================= Fit DoG model to data ==========================

compiled_group_fitting  = qsubcompile('group_fitting','toolbox', {'curvefit'});

time            = 1 * 3600; % hours to seconds
space           = 1 * 1024^3; % Gb to bytes


%% 4-9 back trials
cfg                 = [];
cfg.fittingsteps    = 100;
cfg.fixedwidth      = false;

% same location
pooled_sddata_same_loc_49_mat = cell2mat(pooled_sddata_same_loc_49);
cfg.data        = pooled_sddata_same_loc_49_mat(:,[3 4]);
cfg.savePath    = ['Matlab datafiles/group fits/group_fits_same_loc_4-9-back.mat'];
jobidarray{2}   = qsubfeval(compiled_group_fitting,cfg,'timreq',time,'memreq',space);

% different location
pooled_sddata_diff_loc_49_mat = cell2mat(pooled_sddata_diff_loc_49);
cfg.data        = pooled_sddata_diff_loc_49_mat(:,[3 4]);
cfg.savePath    = ['Matlab datafiles/group fits/group_fits_diff_loc_4-9-back.mat'];
jobidarray{3}   = qsubfeval(compiled_group_fitting,cfg,'timreq',time,'memreq',space);


%% ======================= Run permutation test ===========================
%% =============== Difference between spatial locations ===================

% precompile permutation function for cluster execution
compiled_condition_difference_permutations  = qsubcompile('condition_difference_permutations','toolbox', {'curvefit'});

cfg                 = [];
cfg.subjects        = subjects;
cfg.fixedwidth      = false;
cfg.nperms          = 100;
cfg.fittingsteps    = 1;

time                = 1 * 3600; % hours to seconds
space               = 1 * 1024^3; % Gb to bytes

for i = 1:100
    
    % 4-9 back serial dependence
    cfg.data_cond1   = pooled_sddata_same_loc_49;
    cfg.data_cond2   = pooled_sddata_diff_loc_49;
    cfg.savePath     = ['Matlab datafiles/group fits/pooled permutations/permutation_fits_location_diff_4-9-back_' num2str(i) '.mat'];
    
    jobidarray{1}    = qsubfeval(compiled_condition_difference_permutations,cfg,'timreq',time,'memreq',space);
      
end


%% ===================== Evaluate permutation tests =======================
pooled_fit_same_loc = load(['Matlab datafiles/group fits/group_fits_same_loc_4-9-back.mat']);
pooled_fit_diff_loc = load(['Matlab datafiles/group fits/group_fits_diff_loc_4-9-back.mat']);

amp_diff = pooled_fit_same_loc.group_fit.coeffs(1) - pooled_fit_diff_loc.group_fit.coeffs(1);

perm_amp_diff = [];

for i = 1:100
   
    load(['Matlab datafiles/group fits/pooled permutations/permutation_fits_location_diff_4-9-back_' num2str(i) '.mat'])
    perm_amp_diff = [perm_amp_diff; permutation_distribution_amplitude_difference];
    
    
end

pval = sum(perm_amp_diff <= amp_diff) / length(perm_amp_diff) * 2;



%% ======================= Run permutation test ===========================
%% ======================= Difference from zero  ==========================
%% ===================================================================== %%

% precompile permutation function for cluster execution
compiled_group_permutations  = qsubcompile('group_permutations','toolbox', {'curvefit'});
time = 36 * 3600; % hours to seconds
space = 1 * 1024^3; % Gb to bytes


cfg                 = [];
cfg.subjects        = subjects;
cfg.nperms          = 10000;
cfg.fittingsteps    = 1;
cfg.fixedwidth      = false;

% same loc
cfg.data        = pooled_sddata_same_loc_49;
cfg.savePath    = ['Matlab datafiles/group fits/pooled permutations/same_loc_permutation_fits_4-9-back.mat'];
jobidarray{1}   = qsubfeval(compiled_group_permutations,cfg,'timreq',time,'memreq',space);
    
% diff loc
cfg.data        = pooled_sddata_diff_loc_49;
cfg.savePath    = ['Matlab datafiles/group fits/pooled permutations/diff_loc_permutation_fits_4-9-back.mat'];
jobidarray{2}   = qsubfeval(compiled_group_permutations,cfg,'timreq',time,'memreq',space);


%% Evaluation
pooled_fit_same_loc = load(['Matlab datafiles/group fits/group_fits_same_loc_4-9-back.mat']);
pooled_fit_diff_loc = load(['Matlab datafiles/group fits/group_fits_diff_loc_4-9-back.mat']);

load(['Matlab datafiles/group fits/pooled permutations/same_loc_permutation_fits_4-9-back.mat']);
permutation_distribution_same_loc = permutation_distribution;

load(['Matlab datafiles/group fits/pooled permutations/diff_loc_permutation_fits_4-9-back.mat']);
permutation_distribution_diff_loc = permutation_distribution;

%% compute p-values
p_val_same_loc = sum(permutation_distribution_same_loc(:,1) <= pooled_fit_same_loc.group_fit.coeffs(1)) / length(permutation_distribution_same_loc) * 2;
p_val_diff_loc = sum(permutation_distribution_diff_loc(:,1) <= pooled_fit_diff_loc.group_fit.coeffs(1)) / length(permutation_distribution_diff_loc) * 2;


%% ========================================================================
%% Plot pooled n-back model fits

dog_fun = @(a,w,x) x.*a.*w.*(sqrt(2)./exp(-0.5)).*exp(-(w.*x).^2);

pooled_fit_same_loc = load(['Matlab datafiles/group fits/group_fits_same_loc_4-9-back.mat']);
pooled_fit_diff_loc = load(['Matlab datafiles/group fits/group_fits_diff_loc_4-9-back.mat']);


figure; hold on;
h1 = plot(-90:90,mean(cell2mat(pooled_moving_average_same_loc_49)),'-b');
h2 = plot(-90:90,mean(cell2mat(pooled_moving_average_diff_loc_49)),'-r');
shadedErrorBar(-90:90,mean(cell2mat(pooled_moving_average_same_loc_49)),std(cell2mat(pooled_moving_average_same_loc_49))/sqrt(24),'-b',0);
shadedErrorBar(-90:90,mean(cell2mat(pooled_moving_average_diff_loc_49)),std(cell2mat(pooled_moving_average_diff_loc_49))/sqrt(24),'-r',0);


plot([-90:0.01:90], feval(dog_fun, ...
        pooled_fit_same_loc.group_fit.coeffs(1), pooled_fit_same_loc.group_fit.coeffs(2),...
        [-90:0.01:90]), '-b','LineWidth',3);

plot([-90:0.01:90], feval(dog_fun, ...
        pooled_fit_diff_loc.group_fit.coeffs(1), pooled_fit_diff_loc.group_fit.coeffs(2),...
        [-90:0.01:90]), '-r','LineWidth',3);

ylim([-1.5 1.5]);
xlim([-90 90])
grid on;
    
xlabel('prev. ori.')
ylabel('response error')

legend([h1 h2],'same location','different location')