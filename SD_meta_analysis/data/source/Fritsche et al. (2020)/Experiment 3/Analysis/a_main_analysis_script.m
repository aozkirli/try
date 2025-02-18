addpath /home/common/matlab/fieldtrip/qsub/
addpath /home/predatt/matfri/Toolboxes/CircStat2012a

clear all; close all;

subjects = [1:7,9:25];

%% Load data
cfg = [];
cfg.subjects = subjects;

[data, log] = load_data(cfg);

%{
% compute number of missed trials
for iSubject = 1:24
    missed_trials(iSubject) = sum(isnan(data{iSubject}(:,6)));   
end
%}

%% Outlier removal
cfg                     = [];
cfg.data                = data;
cfg.subjects            = subjects;
cfg.plot                = false;

[data, data_removed]    = outlier_removal(cfg);

%% Compute error stats
cfg = [];
cfg.data = data;

errorstats = compute_error_stats(cfg);


%% Sort SD data
cfg             = [];
cfg.subjects    = subjects;
cfg.data        = data;
cfg.demean      = true; % demean response errors to remove overall bias
cfg.nBack       = 10;

[sd_data] = sort_data(cfg);


%% Compute moving averages
cfg             = [];
cfg.data        = sd_data;
cfg.bin_width   = 20;% possible to specify multiple window sizes

moving_averages = compute_moving_averages(cfg);


%% ===================================================================== %%
%% Fitting to  Dog Models to group data (pooled datapoints)
compiled_group_fitting  = qsubcompile('group_fitting','toolbox', {'curvefit'});

time            = 1 * 3600; % hours to seconds
space           = 1 * 1024^3; % Gb to bytes

cfg                 = [];
cfg.fittingsteps    = 100;
cfg.fixedwidth      = false;

for iBack = 1:10
    
    data_all            = cell2mat(sd_data.data_nback{iBack});
    data_current_short  = cell2mat(sd_data.data_nback_short{iBack});
    data_current_long   = cell2mat(sd_data.data_nback_long{iBack});
    
    % all trials
    cfg.data        = data_all(:,[3 4]);
    cfg.savePath    = ['Matlab datafiles/group fits/group_fits_all_' num2str(iBack) '-back.mat'];
    jobidarray{1}   = qsubfeval(compiled_group_fitting,cfg,'timreq',time,'memreq',space);
       
    % current response delay short
    cfg.data        = data_current_short(:,[3 4]);
    cfg.savePath    = ['Matlab datafiles/group fits/group_fits_current_short_' num2str(iBack) '-back.mat'];
    jobidarray{2}   = qsubfeval(compiled_group_fitting,cfg,'timreq',time,'memreq',space);
    
    % current response delay long
    cfg.data        = data_current_long(:,[3 4]);
    cfg.savePath    = ['Matlab datafiles/group fits/group_fits_current_long_' num2str(iBack) '-back.mat'];
    jobidarray{3}   = qsubfeval(compiled_group_fitting,cfg,'timreq',time,'memreq',space);
       
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
    cfg.data        = sd_data.data_nback{iBack};
    cfg.savePath    = ['Matlab datafiles/group fits/all_permutation_fits_' num2str(iBack) '-back.mat'];
    jobidarray{1}   = qsubfeval(compiled_group_permutations,cfg,'timreq',time,'memreq',space);   
    
    % current short
    cfg.data        = sd_data.data_nback_short{iBack};
    cfg.savePath    = ['Matlab datafiles/group fits/current_short_permutation_fits_' num2str(iBack) '-back.mat'];
    jobidarray{2}   = qsubfeval(compiled_group_permutations,cfg,'timreq',time,'memreq',space);   
    
    % current long
    cfg.data        = sd_data.data_nback_long{iBack};
    cfg.savePath    = ['Matlab datafiles/group fits/current_long_permutation_fits_' num2str(iBack) '-back.mat'];
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
        p_val_all(iBack,1) = sum(permutation_distribution_all(:,1) <= params_all(iBack,1)) / length(permutation_distribution_all); 
    else
        p_val_all(iBack,1) = sum(permutation_distribution_all(:,1) >= params_all(iBack,1)) / length(permutation_distribution_all);
    end
     
end

p_val_all = p_val_all * 2;


%% ===================================================================== %%
%% Permutation test for testing difference of SD between current response delay levels

% precompile permutation function for cluster execution
compiled_condition_difference_permutations  = qsubcompile('condition_difference_permutations','toolbox', {'curvefit'});


cfg                 = [];
cfg.subjects        = subjects;

for iBack = 1:10
    
    cfg.data_cond1   = sd_data.data_nback_short{iBack};
    cfg.data_cond2   = sd_data.data_nback_long{iBack};
    
    cfg.fixedwidth      = false;
    
    cfg.nperms          = 10000;
    cfg.fittingsteps    = 1;
    cfg.savePath        = ['Matlab datafiles/group fits/permutation_fits_current_delay_' num2str(iBack) '-back.mat'];
    
    time                = 36 * 3600; % hours to seconds
    space               = 1 * 1024^3; % Gb to bytes
    
    jobidarray{iBack}    = qsubfeval(compiled_condition_difference_permutations,cfg,'timreq',time,'memreq',space);   
    
end


%% ===================================================================== %%
%% Results of permutation test of differences between response delays

for iBack = 1:10
   
    %% load fitting parameters
    load(['Matlab datafiles/group fits/group_fits_current_long_' num2str(iBack) '-back.mat']);
    params_current_long(iBack,:) = group_fit.coeffs;
    
    load(['Matlab datafiles/group fits/group_fits_current_short_' num2str(iBack) '-back.mat']);
    params_current_short(iBack,:) = group_fit.coeffs;
    
    amp_diff(iBack,1) = params_current_long(iBack,1) - params_current_short(iBack,1);
    
    
    %% load permutation distribution    
    load(['Matlab datafiles/group fits/permutation_fits_current_delay_' num2str(iBack) '-back.mat']);

    %% compute p-values
    p_val_delay_diff(iBack,1) = sum(permutation_distribution_amplitude_difference >= amp_diff(iBack,1)) / length(permutation_distribution_amplitude_difference); 
    
    
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
    cfg.data                = sd_data.data_nback{iBack};
    cfg.savePath            = ['Matlab datafiles/pooled_observer_bootstrap_all_' num2str(iBack) '-back.mat'];
    jobidarray{1}           = qsubfeval(compiled_bootstrap_dog,cfg,'timreq',time,'memreq',space);  
    
    % Short delay
    cfg.data                = sd_data.data_nback_short{iBack};
    cfg.savePath            = ['Matlab datafiles/pooled_observer_bootstrap_short_' num2str(iBack) '-back.mat'];
    jobidarray{1}           = qsubfeval(compiled_bootstrap_dog,cfg,'timreq',time,'memreq',space);    
    
    % Long delay
    cfg.data                = sd_data.data_nback_long{iBack};
    cfg.savePath            = ['Matlab datafiles/pooled_observer_bootstrap_long_' num2str(iBack) '-back.mat'];
    jobidarray{2}           = qsubfeval(compiled_bootstrap_dog,cfg,'timreq',time,'memreq',space);
    
end


%% ========================================================================
%% ========================= Bar plot =====================================

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
%% ==== Permutation test for delay differences of avg n-back trials ==== %%

%% ======================= Fit DoG model to data ==========================

compiled_group_fitting  = qsubcompile('group_fitting','toolbox', {'curvefit'});

time            = 1 * 3600; % hours to seconds
space           = 1 * 1024^3; % Gb to bytes


%% 2-6 back trials

data_pooled_short = [cell2mat(sd_data.data_nback_short{2});...
    cell2mat(sd_data.data_nback_short{3});...
    cell2mat(sd_data.data_nback_short{4});...
    cell2mat(sd_data.data_nback_short{5});...
    cell2mat(sd_data.data_nback_short{6})];

data_pooled_long = [cell2mat(sd_data.data_nback_long{2});...
    cell2mat(sd_data.data_nback_long{3});...
    cell2mat(sd_data.data_nback_long{4});...
    cell2mat(sd_data.data_nback_long{5});...
    cell2mat(sd_data.data_nback_long{6})];

cfg                 = [];
cfg.fittingsteps    = 100;
cfg.fixedwidth      = false;

% short response delay
cfg.data        = data_pooled_short(:,[3 4]);
cfg.savePath    = ['Matlab datafiles/group fits/group_fits_current_short_2-6-back.mat'];
jobidarray{2}   = qsubfeval(compiled_group_fitting,cfg,'timreq',time,'memreq',space);

% long response delay
cfg.data        = data_pooled_long(:,[3 4]);
cfg.savePath    = ['Matlab datafiles/group fits/group_fits_current_long_2-6-back.mat'];
jobidarray{3}   = qsubfeval(compiled_group_fitting,cfg,'timreq',time,'memreq',space);


%% ======================= Run permutation test= ==========================

%% pool the n-back data into subject cells
for iSubject = 1:length(subjects)
    
    pooled_sddata_short_26{iSubject,1} = [sd_data.data_nback_short{2}{iSubject};...
        sd_data.data_nback_short{3}{iSubject};...
        sd_data.data_nback_short{4}{iSubject};...
        sd_data.data_nback_short{5}{iSubject};...
        sd_data.data_nback_short{6}{iSubject}];
    
    
    pooled_sddata_long_26{iSubject,1} = [sd_data.data_nback_long{2}{iSubject};...
        sd_data.data_nback_long{3}{iSubject};...
        sd_data.data_nback_long{4}{iSubject};...
        sd_data.data_nback_long{5}{iSubject};...
        sd_data.data_nback_long{6}{iSubject}]; 
    
    % also pool the moving averages
    pooled_moving_average_short_26{iSubject,1} = mean([moving_averages.short{2}{iSubject};...
        moving_averages.short{3}{iSubject};...
        moving_averages.short{4}{iSubject};...
        moving_averages.short{5}{iSubject};...
        moving_averages.short{6}{iSubject}],1);
    
    pooled_moving_average_long_26{iSubject,1} = mean([moving_averages.long{2}{iSubject};...
        moving_averages.long{3}{iSubject};...
        moving_averages.long{4}{iSubject};...
        moving_averages.long{5}{iSubject};...
        moving_averages.long{6}{iSubject}],1); 
    
end


% precompile permutation function for cluster execution
compiled_condition_difference_permutations  = qsubcompile('condition_difference_permutations','toolbox', {'curvefit'});

cfg                 = [];
cfg.subjects        = subjects;
cfg.fixedwidth      = false;
cfg.nperms          = 100;
cfg.fittingsteps    = 1;

time                = 2 * 3600; % hours to seconds
space               = 1 * 1024^3; % Gb to bytes


for i = 1:100
    
    % 2-6 back serial dependence
    cfg.data_cond1   = pooled_sddata_short_26;
    cfg.data_cond2   = pooled_sddata_long_26;
    cfg.savePath     = ['Matlab datafiles/group fits/pooled permutations/permutation_fits_current_delay_2-6-back_' num2str(i) '.mat'];
    
    jobidarray{i}    = qsubfeval(compiled_condition_difference_permutations,cfg,'timreq',time,'memreq',space);
    
    
end


%% ========================================================================
%% Results of permutation test

%% 2-6 back

load('Matlab datafiles/group fits/group_fits_current_short_2-6-back.mat');
params_short = group_fit.coeffs;

load('Matlab datafiles/group fits/group_fits_current_long_2-6-back.mat');
params_long = group_fit.coeffs;

amp_diff = params_short(1) - params_long(1);

% load permutation distribution
perm_amp_diff = [];

for i = 1:100

    load(['Matlab datafiles/group fits/pooled permutations/permutation_fits_current_delay_2-6-back_' num2str(i) '.mat']);
    perm_amp_diff = [perm_amp_diff; permutation_distribution_amplitude_difference];
    
end
    
pval = sum(perm_amp_diff >= amp_diff) / length(perm_amp_diff);

%% ===================================================================== %%
%% ================= Permutation test against zero ===================== %%
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
cfg.data        = pooled_sddata_short_26;
cfg.savePath    = ['Matlab datafiles/group fits/pooled permutations/short_delay_permutation_fits_2-6-back.mat'];
jobidarray{1}   = qsubfeval(compiled_group_permutations,cfg,'timreq',time,'memreq',space);
    
% diff loc
cfg.data        = pooled_sddata_long_26;
cfg.savePath    = ['Matlab datafiles/group fits/pooled permutations/long_delay_permutation_fits_2-6-back.mat'];
jobidarray{2}   = qsubfeval(compiled_group_permutations,cfg,'timreq',time,'memreq',space);


%% Evaluation
pooled_fit_short    = load('Matlab datafiles/group fits/group_fits_current_short_2-6-back.mat');
pooled_fit_long     = load('Matlab datafiles/group fits/group_fits_current_long_2-6-back.mat');

load(['Matlab datafiles/group fits/pooled permutations/short_delay_permutation_fits_2-6-back.mat']);
permutation_distribution_short = permutation_distribution;

load(['Matlab datafiles/group fits/pooled permutations/long_delay_permutation_fits_2-6-back.mat']);
permutation_distribution_long = permutation_distribution;

%% compute p-values
p_val_same_loc = sum(permutation_distribution_short(:,1) <= pooled_fit_short.group_fit.coeffs(1)) / length(permutation_distribution_short);
p_val_diff_loc = sum(permutation_distribution_long(:,1) <= pooled_fit_long.group_fit.coeffs(1)) / length(permutation_distribution_long);


%% ===================================================================== %%
%% == Non-parametric quantifcation of pooled n-back delay dependence === %%
%% ===================================================================== %%
bins = [-90 0; 0 90];


for iSubject = 1:length(subjects)
    
    thisdata = pooled_sddata_short_26{iSubject};
    %thisdata = sd_data.data_nback_short{1, 1}{iSubject};
    
    bias(iSubject,1) = (mean(thisdata(thisdata(:,3) > 0, 4)) - mean(thisdata(thisdata(:,3) < 0, 4)))/2;
    
    thisdata = sd_data.data_nback_long{1, 1}{iSubject};
    
    bias(iSubject,2) = (mean(thisdata(thisdata(:,3) > 0, 4)) - mean(thisdata(thisdata(:,3) < 0, 4)))/2;
    
   
end
