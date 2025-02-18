clear all; close all;
addpath /home/common/matlab/fieldtrip/qsub/
addpath CircStat2012a/

%% Main analysis of the serial dependence experiment

%% Main parameters of the analysis
subjects            = 1:23;
remove_outliers     = true;

%% =================== load the raw data ================================ %
cfg                 = [];
cfg.dataPath        = 'C:\Users\chare\Google Drive\Work\11_WorkInProgress\EPFL_Ambizione\Datasets\DSC_3018012.11_596_v1\Experiment 1\Data\';%'../Data/';
cfg.subjects        = subjects; % Subjects to be included in the analysis

[data, log]         = load_data(cfg);

%% ===================== Outlier removal   ============================== %
if remove_outliers
    cfg                     = [];
    cfg.data                = data;
    cfg.subjects            = subjects;
    
    [data, data_removed]    = outlier_removal(cfg);
end

%% ================== Compute variance of response errors =============== %
cfg         = [];
cfg.data    = data;

errorstats = compute_error_stats(cfg);


%% ===================== create SD sorted data ========================== %
cfg         = [];
cfg.data    = data;
cfg.nBack   = 40;
cfg.demean  = true;

sd_data = sort_data(cfg);


%% ===================================================================== %%
%% Fitting to Dog Models to group data (pooled datapoints)
compiled_group_fitting  = qsubcompile('group_fitting','toolbox', {'curvefit'});

time            = 1 * 3600; % hours to seconds
space           = 1 * 1024^3; % Gb to bytes

cfg                 = [];
cfg.fittingsteps    = 100;
cfg.fixedwidth      = false;

for iBack = 1:40
    
    data_all            = cell2mat(sd_data.all{iBack});
    
    % all trials
    cfg.data        = data_all(:,[3 4]);
    cfg.savePath    = ['Matlab datafiles/group fits/group_fits_all_' num2str(iBack) '-back.mat'];
    jobidarray{1}   = qsubfeval(compiled_group_fitting,cfg,'timreq',time,'memreq',space);
    
    pause(1)
    
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

for iBack = 1:40
    
    cfg.data                = sd_data.all{iBack};
    cfg.savePath            = ['Matlab datafiles/bootstrap/pooled_observer_bootstrap_all_' num2str(iBack) '-back.mat'];
    jobidarray{1}           = qsubfeval(compiled_bootstrap_dog,cfg,'timreq',time,'memreq',space);  
    
    pause(1)
    
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

for iBack = 1:40
    
    cfg.data        = sd_data.all{iBack};
    cfg.savePath    = ['Matlab datafiles/permutation/all_permutation_fits_' num2str(iBack) '-back.mat'];
    jobidarray{1}   = qsubfeval(compiled_group_permutations,cfg,'timreq',time,'memreq',space); 
    
    pause(1)
    
end


%% ===================================================================== %%
%% Results of permutation test of SD against zero

for iBack = 1:40
   
    %% load fitting parameters
    
    load(['Matlab datafiles/group fits/group_fits_all_' num2str(iBack) '-back.mat']);
    params_all(iBack,:) = group_fit.coeffs;
        
    %% load permutation distribution
    load(['Matlab datafiles/permutation/all_permutation_fits_' num2str(iBack) '-back.mat']);
    permutation_distribution_all = permutation_distribution;
    
    %% compute p-values
    if params_all(iBack,1) <= 0
        p_val_all(iBack,1)      = sum(permutation_distribution_all(:,1) <= params_all(iBack,1)) / length(permutation_distribution_all); 
    else
        p_val_all(iBack,1)      = sum(permutation_distribution_all(:,1) >= params_all(iBack,1)) / length(permutation_distribution_all); 
    end
end

p_val_all = p_val_all * 2;

significance = [p_val_all < 0.05, ...
                p_val_all < 0.05 / 40];


%% ===================================================================== %%
%% ========================== Visualization ============================ %%
%% ===================================================================== %%

%% Compute moving averages
cfg             = [];
cfg.data        = sd_data;
cfg.bin_width   = 20; % possible to specify multiple window sizes
        
[moving_averages] = compute_moving_averages(cfg);


%% load fitting parameters
for iBack = 1:40
   
    load(['Matlab datafiles/group fits/group_fits_all_' num2str(iBack) '-back.mat']);
    params_all(iBack,:) = group_fit.coeffs;
  
end

%% compute standard errors of moving averages
for iBack = 1:40
   
    this_mov_avg = cell2mat(moving_averages.all{iBack});
    sem_mov_avg(iBack,:) = std(this_mov_avg) / sqrt(length(subjects));
  
end

dog_fun = @(a,w,x) x.*a.*w.*(sqrt(2)./exp(-0.5)).*exp(-(w.*x).^2);

figure; hold on;
shadedErrorBar(-90:90,moving_averages.grand_all{1},sem_mov_avg(1,:),'-k',0);
plot([-90:0.01:90], feval(dog_fun, ...
        params_all(1,1), params_all(1,2),...
        [-90:0.01:90]), '-k','LineWidth',3);
ylim([-3 3])
xlim([-90 90])


figure; hold on;
shadedErrorBar(-90:90,moving_averages.grand_all{4},sem_mov_avg(4,:),'-k',0);
plot([-90:0.01:90], feval(dog_fun, ...
        params_all(4,1), params_all(4,2),...
        [-90:0.01:90]), '-k','LineWidth',3);
ylim([-3 3])
xlim([-90 90])



%% ========================================================================
%% ========================= Bar plots ====================================




%% load fits and compute error bars
for iBack = 1:40
   
    load(['Matlab datafiles/group fits/group_fits_all_' num2str(iBack) '-back.mat']);
    params_all(iBack,:) = group_fit.coeffs;
    
    load(['Matlab datafiles/bootstrap/pooled_observer_bootstrap_all_' num2str(iBack) '-back.mat']);
    sem_all(iBack,:) = std(bootstrap_params(:,1));
     
end

figure; hold on;
bar(1:40,params_all(:,1)');
errorbar(1:40,params_all(:,1)',sem_all,'.k','LineWidth',2,'Capsize',0)
xlabel('n-back')
ylabel('Bias amplitude (deg)')


