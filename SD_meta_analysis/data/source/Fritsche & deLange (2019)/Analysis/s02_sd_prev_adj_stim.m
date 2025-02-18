addpath /home/common/matlab/fieldtrip/qsub/

clear all; close all;

%% ===================================================================== %%
%% =================== Subject information ============================= %%
%% ===================================================================== %%
subjects = [2,3,4,5,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,25,26,...
    27,28,29,30,31,32,33,35,36,37,38];


%% ========================================================================
%% Load data
%% ========================================================================
cfg = [];
cfg.subjects = subjects;

[data_ori,~,data_size,~] = load_data(cfg);


%% ===================================================================== %%
%% ========= Compute SD data (for previous adjusmtent stimulus)========= %%
%% ===================================================================== %%

for iSubject = 1:length(subjects)
   
    %% Attention to orientation

    % compute adjustment errors
    adj_errors = mod(data_ori{iSubject}(:,11) - data_ori{iSubject}(:,10) + 90, 180) - 90;

    % compute upper and lower bounds for outliers
    sd = 3 * std(adj_errors,1);
    lower_bound = mean(adj_errors, 1) - sd;
    upper_bound = mean(adj_errors, 1) + sd;

    % reject outliers
    idx_reject = find(adj_errors < lower_bound | adj_errors > upper_bound);

    adj_errors(idx_reject) = nan;
        
    % demean adjustment errors
    adj_errors = bsxfun(@minus, adj_errors, nanmean(adj_errors));
    
    % Sort SD data
    tmp_data = [];
    
    for j = 2:size(data_ori{iSubject,1},1)
        
        if (data_ori{iSubject,1}(j,2) == data_ori{iSubject,1}(j-1,2))  && (~isnan(adj_errors(j)) && ~isnan(adj_errors(j-1)))
            
            prev_ori   = data_ori{iSubject,1}(j-1,10);
            current_ori = data_ori{iSubject,1}(j,10);
            error = adj_errors(j);
            
            stim_ori_diff = mod(prev_ori - current_ori + 90, 180) - 90;
            
            tmp_data = [tmp_data; stim_ori_diff error];
            
        end
        
    end
    
    sd_data_ori{iSubject,1} = tmp_data;
    
    clear tmp_data
    
    
    %% Attention to size

    % compute adjustment errors
    adj_errors = mod(data_size{iSubject}(:,11) - data_size{iSubject}(:,10) + 90, 180) - 90;

    % compute upper and lower bounds for outliers
    sd = 3 * std(adj_errors,1);
    lower_bound = mean(adj_errors, 1) - sd;
    upper_bound = mean(adj_errors, 1) + sd;

    % reject outliers
    idx_reject = find(adj_errors < lower_bound | adj_errors > upper_bound);

    adj_errors(idx_reject) = nan;
        
    % demean adjustment errors
    adj_errors = bsxfun(@minus, adj_errors, nanmean(adj_errors));
    
    % Sort SD data
    tmp_data = [];
    
    for j = 2:size(data_size{iSubject,1},1)
        
        if (data_size{iSubject,1}(j,2) == data_size{iSubject,1}(j-1,2))  && (~isnan(adj_errors(j)) && ~isnan(adj_errors(j-1)))
            
            prev_ori   = data_size{iSubject,1}(j-1,10);
            current_ori = data_size{iSubject,1}(j,10);
            error = adj_errors(j);
            
            stim_ori_diff = mod(prev_ori - current_ori + 90, 180) - 90;
            
            tmp_data = [tmp_data; stim_ori_diff error];
            
        end
        
    end
    
    sd_data_size{iSubject,1} = tmp_data;
    
    clear tmp_data
    
end


%% ===================================================================== %%
%% ===================== Compute moving averages ======================= %%
%% ===================================================================== %%
bin_width = 20;

for iSubject = 1:length(subjects)
    
    %% Attention to orientation
    x = sd_data_ori{iSubject}(:,1);
    y = sd_data_ori{iSubject}(:,2);
    
    x_padded = [x - 180; x; x + 180];
    y_padded = [y; y; y];
    
    moving_averages.orientation{iSubject,1} = zeros(181,1);
    
    for b = 0:180
        moving_averages.orientation{iSubject,1}(b+1) = ...
            circ_rad2ang(circ_mean(circ_ang2rad(y_padded(inrange(x_padded, [1, bin_width] - floor(bin_width/2) - 1 - 90 + b)))));
    end
    
    %% Attention to size
    x = sd_data_size{iSubject}(:,1);
    y = sd_data_size{iSubject}(:,2);
    
    x_padded = [x - 180; x; x + 180];
    y_padded = [y; y; y];
    
    moving_averages.size{iSubject,1} = zeros(181,1);
    
    for b = 0:180
        moving_averages.size{iSubject,1}(b+1) = ...
            circ_rad2ang(circ_mean(circ_ang2rad(y_padded(inrange(x_padded, [1, bin_width] - floor(bin_width/2) - 1 - 90 + b)))));
    end
    
end

grand_mov_avg_ori = mean(cell2mat(moving_averages.orientation'),2);
grand_mov_avg_size = mean(cell2mat(moving_averages.size'),2);

sem_mov_avg_ori = std(cell2mat(moving_averages.orientation'),1,2) /sqrt(34);
sem_mov_avg_size = std(cell2mat(moving_averages.size'),1,2) /sqrt(34);


%% ===================================================================== %%
%% ==================== Pooled Observer Analysis ======================= %%
%% ===================================================================== %%

%% ================ Fit DoG models to pooled observer data ============= %%
 
% orientation
cfg                 = [];
cfg.data            = cell2mat(sd_data_ori);
cfg.fittingsteps    = 100;
cfg.fixedwidth      = false;

pooled_observer_dog_fit_ori = fit_dog(cfg);

save('Matlab datafiles/prev_adj_stimulus/prev_adj_stim_pooled_observer_dog_fit_ori.mat','pooled_observer_dog_fit_ori');

% size
cfg                 = [];
cfg.data            = cell2mat(sd_data_size);
cfg.fittingsteps    = 100;
cfg.fixedwidth      = false;

pooled_observer_dog_fit_size = fit_dog(cfg);

save('Matlab datafiles/prev_adj_stimulus/prev_adj_stim_pooled_observer_dog_fit_size.mat','pooled_observer_dog_fit_size');


%% ================= Bootstrap DoG parameter estimates ================= %%

% precompile permutation function for cluster execution
compiled_bootstrap_dog  = qsubcompile('bootstrap_dog','toolbox', {'curvefit'});
time = 1 * 3600; % hours to seconds
space = 1 * 1024^3; % Gb to bytes

cfg                     = [];
cfg.subjects            = subjects;
cfg.fittingsteps        = 1;
cfg.nBootstrapSamples   = 1000;
cfg.fixedwidth          = false;

% Orientation
cfg.data                = sd_data_ori;
cfg.savePath            = 'Matlab datafiles/prev_adj_stimulus/pooled_observer_bootstrap_ori.mat';
jobidarray{1}           = qsubfeval(compiled_bootstrap_dog,cfg,'timreq',time,'memreq',space);

% Size
cfg.data                = sd_data_size;
cfg.savePath            = 'Matlab datafiles/prev_adj_stimulus/pooled_observer_bootstrap_size.mat';
jobidarray{2}           = qsubfeval(compiled_bootstrap_dog,cfg,'timreq',time,'memreq',space);


%% ======================== Visualization ============================== %%
load('Matlab datafiles/prev_adj_stimulus/prev_adj_stim_pooled_observer_dog_fit_ori.mat');
load('Matlab datafiles/prev_adj_stimulus/prev_adj_stim_pooled_observer_dog_fit_size.mat');

dog_fun = @(a,w,x) x.*a.*w.*(sqrt(2)./exp(-0.5)).*exp(-(w.*x).^2);

figure; hold on;
shadedErrorBar(-90:90,grand_mov_avg_ori,sem_mov_avg_ori,'b');
shadedErrorBar(-90:90,grand_mov_avg_size,sem_mov_avg_size,'r');

h1 = plot([-90:0.01:90], feval(dog_fun, ...
    pooled_observer_dog_fit_ori.coeffs(1), pooled_observer_dog_fit_ori.coeffs(2),...
    [-90:0.01:90]), '-b','LineWidth',3);

h2 = plot([-90:0.01:90], feval(dog_fun, ...
    pooled_observer_dog_fit_size.coeffs(1), pooled_observer_dog_fit_size.coeffs(2),...
    [-90:0.01:90]), '-r','LineWidth',3);

ylim([-3 3])
xlim([-90 90])

xlabel('Rel. orientation of previous adjustment stimulus')
ylabel('Response Error')

%legend([h1,h2],'attention to orientation','attention to size')

%% Bar plot
load('Matlab datafiles/prev_adj_stimulus/pooled_observer_bootstrap_ori.mat');
bootstrap_sd_amp_ori = std(bootstrap_params(:,1),1,1);
bootstrap_sd_width_ori = std(bootstrap_params(:,2),1,1);

load('Matlab datafiles/prev_adj_stimulus/pooled_observer_bootstrap_size.mat');
bootstrap_sd_amp_size = std(bootstrap_params(:,1),1,1);
bootstrap_sd_width_size = std(bootstrap_params(:,2),1,1);

amps = [pooled_observer_dog_fit_ori.coeffs(1),pooled_observer_dog_fit_size.coeffs(1)];
sem_amps = [bootstrap_sd_amp_ori bootstrap_sd_amp_size];

widths = [pooled_observer_dog_fit_ori.coeffs(2),pooled_observer_dog_fit_size.coeffs(2)];
sem_widths = [bootstrap_sd_width_ori bootstrap_sd_width_size];

figure

subplot(1,2,1); hold on;
bar(1:2,amps,0.5);
errorbar(1:2,amps,sem_amps,'.r','LineWidth',2)
xlim([0.5 2.5])
ylim([0 3])
ylabel('DoG Amplitude (deg)')
title('SD Amplitude')

ax = gca;
ax.XTick = [1 2];
ax.XTickLabel = {'Orientation','Size'};

subplot(1,2,2); hold on;
bar(1:2,widths,0.5);
errorbar(1:2,widths,sem_widths,'.r','LineWidth',2)
xlim([0.5 2.5])
ylabel('DoG Width')
title('SD Width')

ax = gca;
ax.XTick = [1 2];
ax.XTickLabel = {'Orientation','Size'};




%% ===== Compute permutation distributions for test against zero  ====== %%

% precompile permutation function for cluster execution
compiled_group_permutations  = qsubcompile('group_permutations','toolbox', {'curvefit'});
time = 36 * 3600; % hours to seconds
space = 1 * 1024^3; % Gb to bytes

cfg                 = [];
cfg.subjects        = subjects;
cfg.nperms          = 10000;
cfg.fittingsteps    = 1;
cfg.fixedwidth      = false;

% orientation condition
cfg.data        = sd_data_ori;
cfg.savePath    = 'Matlab datafiles/prev_adj_stimulus/prev_adj_stim_pooled_observer_permutations_ori.mat';
jobidarray{1}   = qsubfeval(compiled_group_permutations,cfg,'timreq',time,'memreq',space);

% size condition
cfg.data        = sd_data_size;
cfg.savePath    = 'Matlab datafiles/prev_adj_stimulus/prev_adj_stim_pooled_observer_permutations_size.mat';
jobidarray{2}   = qsubfeval(compiled_group_permutations,cfg,'timreq',time,'memreq',space);

%% ===== Evaluate permutation distributions for test against zero  ===== %%
load('Matlab datafiles/prev_adj_stimulus/prev_adj_stim_pooled_observer_dog_fit_ori.mat');
load('Matlab datafiles/prev_adj_stimulus/prev_adj_stim_pooled_observer_dog_fit_size.mat');

perm_dist_ori = load('Matlab datafiles/prev_adj_stimulus/prev_adj_stim_pooled_observer_permutations_ori.mat');
perm_dist_size = load('Matlab datafiles/prev_adj_stimulus/prev_adj_stim_pooled_observer_permutations_size.mat');

p_val_ori = sum(perm_dist_ori.permutation_distribution(:,1) >= ...
    pooled_observer_dog_fit_ori.coeffs(1)) / length(perm_dist_ori.permutation_distribution(:,1));

p_val_size = sum(perm_dist_size.permutation_distribution(:,1) >= ...
    pooled_observer_dog_fit_size.coeffs(1)) / length(perm_dist_size.permutation_distribution(:,1));


%% === Compute permutation distributions for test between conditions  === %%

% precompile permutation function for cluster execution
compiled_condition_difference_permutations  = qsubcompile('condition_difference_permutations','toolbox', {'curvefit'});

cfg                 = [];
cfg.subjects        = subjects;

cfg.data_cond1   = sd_data_ori;
cfg.data_cond2   = sd_data_size;

cfg.fixedwidth      = false;

cfg.nperms          = 10000;
cfg.fittingsteps    = 1;
cfg.savePath        = 'Matlab datafiles/prev_adj_stimulus/prev_adj_stim_pooled_observer_permutations_condition_difference.mat';

time                = 36 * 3600; % hours to seconds
space               = 1 * 1024^3; % Gb to bytes

jobidarray{3}    = qsubfeval(compiled_condition_difference_permutations,cfg,'timreq',time,'memreq',space);


%% === Evaluate permutation distributions for test between conditions  === %%
load('Matlab datafiles/prev_adj_stimulus/prev_adj_stim_pooled_observer_dog_fit_ori.mat');
load('Matlab datafiles/prev_adj_stimulus/prev_adj_stim_pooled_observer_dog_fit_size.mat');

amp_diff = pooled_observer_dog_fit_ori.coeffs(1) - pooled_observer_dog_fit_size.coeffs(1);
width_diff = pooled_observer_dog_fit_ori.coeffs(2) - pooled_observer_dog_fit_size.coeffs(2);

load('Matlab datafiles/prev_adj_stimulus/prev_adj_stim_pooled_observer_permutations_condition_difference.mat');

p_value_amp_diff = sum(permutation_distribution_amplitude_difference <= amp_diff) / length(permutation_distribution_amplitude_difference);
p_value_width_diff = sum(permutation_distribution_width_difference <= width_diff) / length(permutation_distribution_width_difference);
