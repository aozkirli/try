addpath /home/common/matlab/fieldtrip/qsub/

clear all; close all;

%% ===================================================================== %%
%% =================== Subject information ============================= %%
%% ===================================================================== %%
subjects = [2,3,4,5,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,25,26,...
    27,28,29,30,31,32,33,35,36,37,38];

subjects_ori_first = [2,3,8,9,10,13,14,17,18,21,22,25,26,30,31,35,38];
subjects_ori_first_idx = arrayfun(@(x)find(subjects==x,1),subjects_ori_first);

subjects_size_first = [4,5,7,11,12,15,16,19,20,23,27,28,29,32,33,36,37];
subjects_size_first_idx = arrayfun(@(x)find(subjects==x,1),subjects_size_first);

%% ========================================================================
%% Load data
%% ========================================================================
cfg = [];
cfg.subjects = subjects;

[data_ori,~,data_size,~] = load_data(cfg);

%% ========================================================================
%% Analyze the 2IFC accuracies
%% ========================================================================
ifc_orientation_accuracy    = nan(1,length(subjects));
ifc_size_accuracy           = nan(1,length(subjects));

for iSubject = 1:length(subjects)
    ifc_orientation_accuracy(iSubject) = sum(data_ori{iSubject}(:,13))/length(data_ori{iSubject}(:,13));
    ifc_size_accuracy(iSubject) = sum(data_size{iSubject}(:,13))/length(data_size{iSubject}(:,13));
end

mean_ifc_orientation_accuracy = mean(ifc_orientation_accuracy);
mean_ifc_size_accuracy = mean(ifc_size_accuracy);

sem_ifc_accuracy = withinstde([ifc_orientation_accuracy' ifc_size_accuracy']);

figure; hold on;
plot(repmat([0.5; 1],1,length(mean_ifc_orientation_accuracy)),[ifc_orientation_accuracy; ifc_size_accuracy], '-s')
errorbar([0.5 1],[mean_ifc_orientation_accuracy mean_ifc_size_accuracy],sem_ifc_accuracy,'-sk','LineWidth',3);
ylim([0.5 1])
xlim([0.25 1.25])

ax = gca;
ax.XTick = [0.5 1];
ax.XTickLabel = {'Orientation','Size'};

ylabel('Accuracy')

grid on;
grid minor;


%% ========================================================================
%% Analyze the 2IFC thresholds
%% ========================================================================
for iSubject = 1:length(subjects)
    
    tmp = load(['../Data/S' num2str(subjects(iSubject)) '/orientation_threshold.mat'],'orientation_threshold');
    orientation_threshold(iSubject,1) = tmp.orientation_threshold;
    
    tmp = load(['../Data/S' num2str(subjects(iSubject)) '/size_threshold.mat']);
    size_threshold(iSubject,1) = tmp.size_threshold;
        
end

cutoff_orientation_threshold = mean(orientation_threshold) + 3*std(orientation_threshold,1,1);
cutoff_size_threshold = mean(size_threshold) + 3*std(size_threshold,1,1);



%% ========================================================================
%% Analyze mean adjustment errors
%% ========================================================================
for iSubject = 1:length(subjects)
  
    adj_errors_ori = mod(data_ori{iSubject}(:,11) - data_ori{iSubject}(:,10) + 90, 180) - 90;
    adj_errors_size = mod(data_size{iSubject}(:,11) - data_size{iSubject}(:,10) + 90, 180) - 90;
    
    overall_mean_adj_error_ori(iSubject,1) = mean(abs(adj_errors_ori));
    overall_mean_adj_error_size(iSubject,1) = mean(abs(adj_errors_size));
    
    overall_sd_adjustment_error_ori(iSubject,1) = std(adj_errors_ori,1,1);
    overall_sd_adjustment_error_size(iSubject,1) = std(adj_errors_size,1,1);
    
    overall_mean_adj_error(iSubject,1) = mean([overall_mean_adj_error_ori(iSubject,1),overall_mean_adj_error_size(iSubject,1)]);
    overall_sd_adj_error(iSubject,1) = std([adj_errors_ori; adj_errors_size],1,1);
    
end

cutoff_adjustment_error = mean(overall_mean_adj_error) + 3*std(overall_mean_adj_error,1,1);

sem_mean_adj_error = withinstde([overall_mean_adj_error_ori,overall_mean_adj_error_size]);

figure; hold on;
plot(repmat([0.5; 1],1,length(overall_mean_adj_error_ori)),[overall_mean_adj_error_ori'; overall_mean_adj_error_size'], '-s')
errorbar([0.5 1],[mean(overall_mean_adj_error_ori) mean(overall_mean_adj_error_size)],sem_mean_adj_error,'-sk','LineWidth',3);
%ylim([0.5 1])
xlim([0.25 1.25])

ax = gca;
ax.XTick = [0.5 1];
ax.XTickLabel = {'Orientation','Size'};

ylabel('Accuracy')

grid on;
grid minor;



%% ========================================================================
%% Analyze response times
%% ========================================================================
for iSubject = 1:length(subjects)
  
    subject_rt_2ifc(iSubject,1) = nanmean(data_ori{iSubject}(:,14));
    subject_rt_2ifc(iSubject,2) = nanmean(data_size{iSubject}(:,14));
    
    subject_rt_adjustment(iSubject,1) = mean(data_ori{iSubject}(:,12));
    subject_rt_adjustment(iSubject,2) = mean(data_size{iSubject}(:,12));
   
end

%% 2IFC RT plot
sem_rt_2ifc = withinstde(subject_rt_2ifc);

figure; hold on;
plot(repmat([0.5; 1],1,length(subject_rt_2ifc)),[subject_rt_2ifc(:,1)'; subject_rt_2ifc(:,2)'], '-s')
errorbar([0.5 1],[mean(subject_rt_2ifc(:,1)) mean(subject_rt_2ifc(:,2))],sem_rt_2ifc,'-sk','LineWidth',3);
ylim([0 1.1])
xlim([0.25 1.25])

ax = gca;
ax.XTick = [0.5 1];
ax.XTickLabel = {'Orientation','Size'};

ylabel('2IFC RT')

grid on;
grid minor;


%% Adjustment RT plot
sem_rt_adjustment = withinstde(subject_rt_adjustment);

figure; hold on;
plot(repmat([0.5; 1],1,length(subject_rt_adjustment)),[subject_rt_adjustment(:,1)'; subject_rt_adjustment(:,2)'], '-s')
errorbar([0.5 1],[mean(subject_rt_adjustment(:,1)) mean(subject_rt_adjustment(:,2))],sem_rt_adjustment,'-sk','LineWidth',3);
%ylim([0 1.1])
xlim([0.25 1.25])

ax = gca;
ax.XTick = [0.5 1];
ax.XTickLabel = {'Orientation','Size'};

ylabel('Adjustment RT')

grid on;
grid minor;


%% ===================================================================== %%
%% =================== Participant Exclusion Check ===================== %%
%% ===================================================================== %%
exclude_ifc_accuracy_ori = find(ifc_orientation_accuracy < 0.6);
exclude_ifc_accuracy_size = find(ifc_size_accuracy < 0.6);

exclude_ifc_threshold_ori = find(orientation_threshold > cutoff_orientation_threshold);
exclude_ifc_threshold_size = find(size_threshold > cutoff_size_threshold);

exclude_adjustment_error = find(overall_mean_adj_error > cutoff_adjustment_error);

exclude_subjects = unique([exclude_ifc_accuracy_ori, exclude_ifc_accuracy_size...
    exclude_ifc_threshold_ori, exclude_ifc_threshold_size, exclude_adjustment_error]);

subjects(exclude_subjects) = [];



%% ===================================================================== %%
%% ======== Compute single subject SD data and exclude outliers ======== %%
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
    idx_reject_outliers = find(adj_errors < lower_bound | adj_errors > upper_bound);
    idx_reject_noresponse = find(isnan(data_ori{iSubject}(:,14)));
    
    num_reject_outliers(iSubject,1) = length(idx_reject_outliers);
    num_reject_noresponse(iSubject,1) = length(idx_reject_noresponse);

    tmp = [data_ori{iSubject}(:,9), adj_errors];
    tmp([idx_reject_outliers;idx_reject_noresponse],:) = [];
        
    % demean adjustment errors
    tmp(:,2) = bsxfun(@minus, tmp(:,2), mean(tmp(:,2), 1));
    
    sd_data_ori{iSubject,1} = tmp;
    
    clear tmp
  
    
    %% Attention to Size    
    
    % compute adjustment errors
    adj_errors = mod(data_size{iSubject}(:,11) - data_size{iSubject}(:,10) + 90, 180) - 90;
 
    % compute upper and lower bounds for outliers
    sd = 3 * std(adj_errors,1);
    lower_bound = mean(adj_errors, 1) - sd;
    upper_bound = mean(adj_errors, 1) + sd;

    % reject outliers
    idx_reject_outliers = find(adj_errors < lower_bound | adj_errors > upper_bound);
    idx_reject_noresponse = find(isnan(data_size{iSubject}(:,14)));
    
    num_reject_outliers(iSubject,2) = length(idx_reject_outliers);
    num_reject_noresponse(iSubject,2) = length(idx_reject_noresponse);

    tmp = [data_size{iSubject}(:,9), adj_errors];
    tmp([idx_reject_outliers;idx_reject_noresponse],:) = [];
        
    % demean adjustment errors
    tmp(:,2) = bsxfun(@minus, tmp(:,2), mean(tmp(:,2), 1));
    
    sd_data_size{iSubject,1} = tmp;
    
    clear tmp
    
end

%% Compute trial exclusion numbers
mean_num_reject_outliers = mean(num_reject_outliers);
std_num_reject_outliers = std(num_reject_outliers,1,1);

mean_num_reject_noresponse = mean(num_reject_noresponse);
std_num_reject_noresponse = std(num_reject_noresponse,1,1);

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

save('Matlab datafiles/pooled_observer_dog_fit_ori.mat','pooled_observer_dog_fit_ori');

% size
cfg                 = [];
cfg.data            = cell2mat(sd_data_size);
cfg.fittingsteps    = 100;
cfg.fixedwidth      = false;

pooled_observer_dog_fit_size = fit_dog(cfg);

save('Matlab datafiles/pooled_observer_dog_fit_size.mat','pooled_observer_dog_fit_size');


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
cfg.savePath            = 'Matlab datafiles/pooled_observer_bootstrap_ori.mat';
jobidarray{1}           = qsubfeval(compiled_bootstrap_dog,cfg,'timreq',time,'memreq',space);

% Size
cfg.data                = sd_data_size;
cfg.savePath            = 'Matlab datafiles/pooled_observer_bootstrap_size.mat';
jobidarray{2}           = qsubfeval(compiled_bootstrap_dog,cfg,'timreq',time,'memreq',space);


%% ======================== Visualization ============================== %%
load('Matlab datafiles/pooled_observer_dog_fit_ori.mat');
load('Matlab datafiles/pooled_observer_dog_fit_size.mat');

cfg = [];
cfg.subjects = subjects;
cfg.sd_data_ori = sd_data_ori;
cfg.sd_data_size = sd_data_size;

[mean_datapoints,sem_datapoints,mean_smoothed_datapoints,sem_smoothed_datapoints] = pooled_observer_datapoints(cfg);

dog_fun = @(a,w,x) x.*a.*w.*(sqrt(2)./exp(-0.5)).*exp(-(w.*x).^2);

figure; hold on;
errorbar(-90:10:90,[mean_smoothed_datapoints(1,end), mean_smoothed_datapoints(1,:)],[sem_smoothed_datapoints(1,end), sem_smoothed_datapoints(1,:)],'-b','LineWidth',2);
errorbar(-90:10:90,[mean_smoothed_datapoints(2,end), mean_smoothed_datapoints(2,:)],[sem_smoothed_datapoints(2,end), sem_smoothed_datapoints(2,:)],'-r','LineWidth',2);

plot([-90:0.01:90], feval(dog_fun, ...
    pooled_observer_dog_fit_ori.coeffs(1), pooled_observer_dog_fit_ori.coeffs(2),...
    [-90:0.01:90]), '-b','LineWidth',3);

plot([-90:0.01:90], feval(dog_fun, ...
    pooled_observer_dog_fit_size.coeffs(1), pooled_observer_dog_fit_size.coeffs(2),...
    [-90:0.01:90]), '-r','LineWidth',3);

ylim([-3.5 3.5])
xlim([-90 90])

xlabel('Rel. orientation of inducer')
ylabel('Response Error')

legend('attention to orientation','attention to size')

%% Bar plot
load('Matlab datafiles/pooled_observer_bootstrap_ori.mat');
bootstrap_sd_amp_ori = std(bootstrap_params(:,1),1,1);
bootstrap_sd_width_ori = std(bootstrap_params(:,2),1,1);

load('Matlab datafiles/pooled_observer_bootstrap_size.mat');
bootstrap_sd_amp_size = std(bootstrap_params(:,1),1,1);
bootstrap_sd_width_size = std(bootstrap_params(:,2),1,1);

amps = [pooled_observer_dog_fit_ori.coeffs(1),pooled_observer_dog_fit_size.coeffs(1)];
sem_amps = [bootstrap_sd_amp_ori bootstrap_sd_amp_size];

widths = [pooled_observer_dog_fit_ori.coeffs(2),pooled_observer_dog_fit_size.coeffs(2)];
sem_widths = [bootstrap_sd_width_ori bootstrap_sd_width_size];

figure

%subplot(1,2,1); 
hold on;
bar(1:2,amps,0.5);
errorbar(1:2,amps,sem_amps,'.r','LineWidth',2)
xlim([0.5 2.5])
ylabel('DoG Amplitude (deg)')
title('SD Amplitude')

ax = gca;
ax.XTick = [1 2];
ax.XTickLabel = {'Orientation','Size'};

%subplot(1,2,2);
figure
hold on;
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
cfg.savePath    = 'Matlab datafiles/pooled_observer_permutations_ori.mat';
jobidarray{1}   = qsubfeval(compiled_group_permutations,cfg,'timreq',time,'memreq',space);

% size condition
cfg.data        = sd_data_size;
cfg.savePath    = 'Matlab datafiles/pooled_observer_permutations_size.mat';
jobidarray{2}   = qsubfeval(compiled_group_permutations,cfg,'timreq',time,'memreq',space);


%% ===== Evaluate permutation distributions for test against zero  ===== %%
load('Matlab datafiles/pooled_observer_dog_fit_ori.mat');
load('Matlab datafiles/pooled_observer_dog_fit_size.mat');

perm_dist_ori = load('Matlab datafiles/pooled_observer_permutations_ori.mat');
perm_dist_size = load('Matlab datafiles/pooled_observer_permutations_size.mat');

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
cfg.savePath        = 'Matlab datafiles/pooled_observer_permutations_condition_difference.mat';

time                = 36 * 3600; % hours to seconds
space               = 1 * 1024^3; % Gb to bytes

jobidarray{3}    = qsubfeval(compiled_condition_difference_permutations,cfg,'timreq',time,'memreq',space);


%% === Evaluate permutation distributions for test between conditions  === %%
load('Matlab datafiles/pooled_observer_dog_fit_ori.mat');
load('Matlab datafiles/pooled_observer_dog_fit_size.mat');

amp_diff = pooled_observer_dog_fit_ori.coeffs(1) - pooled_observer_dog_fit_size.coeffs(1);
width_diff = pooled_observer_dog_fit_ori.coeffs(2) - pooled_observer_dog_fit_size.coeffs(2);

load('Matlab datafiles/pooled_observer_permutations_condition_difference.mat');

p_value_amp_diff = sum(permutation_distribution_amplitude_difference >= amp_diff) / length(permutation_distribution_amplitude_difference);
p_value_width_diff = sum(permutation_distribution_width_difference <= width_diff) / length(permutation_distribution_width_difference);




%% ===================================================================== %%
%% ================ Individual Observer Analysis ======================= %%
%% ===================================================================== %%
showPlot = false;

for iSubject = 1:length(subjects)  
    
    disp(['Fitting subject ' num2str(iSubject)]);
    
    cfg                 = [];
    cfg.data            = sd_data_ori{iSubject};
    cfg.fittingsteps    = 10;
    cfg.fixedwidth      = false;
    
    fit_ori = fit_dog(cfg);
    
    cfg                 = [];
    cfg.data            = sd_data_size{iSubject};
    cfg.fittingsteps    = 10;
    cfg.fixedwidth      = false;
    
    fit_size = fit_dog(cfg);
    
    save(['Matlab datafiles/subject fits/subj-' num2str(iSubject) '_fit.mat'],'fit_ori','fit_size');
    
  
    %% Plot
    if showPlot
        
        figure;
        
        % Attention to orientation
        
        % compute mean adjustment errors for plot
        ifc_ref_stim_relative_orientation = -80:10:90;
        
        for iStimLevel = 1:length(ifc_ref_stim_relative_orientation)
            mean_adj_error(iStimLevel) = mean(sd_data_ori{iSubject}(sd_data_ori{iSubject}(:,1) == ifc_ref_stim_relative_orientation(iStimLevel),2));
        end
        
        for iStimLevel = 2:length(ifc_ref_stim_relative_orientation)-1
            mov_avg_mean_adj_error(iStimLevel) = (mean_adj_error(iStimLevel-1)+ mean_adj_error(iStimLevel) + mean_adj_error(iStimLevel+1))/3;
        end
        
        mov_avg_mean_adj_error(1) = (mean_adj_error(length(ifc_ref_stim_relative_orientation))+ mean_adj_error(1) + mean_adj_error(2))/3;
        mov_avg_mean_adj_error(length(ifc_ref_stim_relative_orientation)) = (mean_adj_error(length(ifc_ref_stim_relative_orientation)-1) + mean_adj_error(length(ifc_ref_stim_relative_orientation)) + mean_adj_error(1))/3;
        
        dog_fun = @(a,w,x) x.*a.*w.*(sqrt(2)./exp(-0.5)).*exp(-(w.*x).^2);
        
        suptitle(['Subject ' num2str(iSubject)])
        
        subplot(1,2,1);
        hold on;
        plot(sd_data_ori{iSubject}(:,1),sd_data_ori{iSubject}(:,2),'.','MarkerSize',10);
        %plot(ifc_ref_stim_relative_orientation,mean_adj_error,'.r','MarkerSize',20);
        plot(ifc_ref_stim_relative_orientation,mov_avg_mean_adj_error,'.k','MarkerSize',20);
        plot([-90:0.01:90], feval(dog_fun, ...
            fit_ori.coeffs(1),fit_ori.coeffs(2),[-90:0.01:90]), '-k','LineWidth',2);
        xlim([-90 90]);
        ylim([-10 10])
        grid on;
        
        ylabel('Response error (deg)')
        xlabel('Relative orientation of first 2IFC stimulus (deg)')
        
        title('Attention to orientation')
        
        
        % Attention to size
        
        % compute mean adjustment errors for plot
        ifc_ref_stim_relative_orientation = -80:10:90;
        
        for iStimLevel = 1:length(ifc_ref_stim_relative_orientation)
            mean_adj_error(iStimLevel) = mean(sd_data_size{iSubject}(sd_data_size{iSubject}(:,1) == ifc_ref_stim_relative_orientation(iStimLevel),2));
        end
        
        for iStimLevel = 2:length(ifc_ref_stim_relative_orientation)-1
            mov_avg_mean_adj_error(iStimLevel) = (mean_adj_error(iStimLevel-1)+ mean_adj_error(iStimLevel) + mean_adj_error(iStimLevel+1))/3;
        end
        
        mov_avg_mean_adj_error(1) = (mean_adj_error(length(ifc_ref_stim_relative_orientation))+ mean_adj_error(1) + mean_adj_error(2))/3;
        mov_avg_mean_adj_error(length(ifc_ref_stim_relative_orientation)) = (mean_adj_error(length(ifc_ref_stim_relative_orientation)-1) + mean_adj_error(length(ifc_ref_stim_relative_orientation)) + mean_adj_error(1))/3;
        
        dog_fun = @(a,w,x) x.*a.*w.*(sqrt(2)./exp(-0.5)).*exp(-(w.*x).^2);
        
        subplot(1,2,2);
        hold on;
        plot(sd_data_size{iSubject}(:,1),sd_data_size{iSubject}(:,2),'.','MarkerSize',10);
        %plot(ifc_ref_stim_relative_orientation,mean_adj_error,'.r','MarkerSize',20);
        plot(ifc_ref_stim_relative_orientation,mov_avg_mean_adj_error,'.k','MarkerSize',20);
        plot([-90:0.01:90], feval(dog_fun, ...
            fit_size.coeffs(1),fit_size.coeffs(2),[-90:0.01:90]), '-k','LineWidth',2);
        xlim([-90 90]);
        ylim([-10 10])
        grid on;
        
        ylabel('Response error (deg)')
        xlabel('Relative orientation of first 2IFC stimulus (deg)')
        
        title('Attention to size')
    end

end

%% Load single-subject DoG parameters and do statistics
for iSubject = 1:length(subjects)
   
    tmp = load(['Matlab datafiles/subject fits/subj-' num2str(iSubject) '_fit.mat']);
    
    ss_params_ori(iSubject,:) = tmp.fit_ori.coeffs;
    ss_params_size(iSubject,:) = tmp.fit_size.coeffs;
    
end


%% Plot

% Amplitude
mean_amp_ori = mean(ss_params_ori(:,1));
mean_amp_size = mean(ss_params_size(:,1));

sem_amp_ori = std(ss_params_ori(:,1),1,1) / sqrt(length(subjects));
sem_amp_size = std(ss_params_size(:,1),1,1) / sqrt(length(subjects));

sem_amp = withinstde([ss_params_ori(:,1) ss_params_size(:,1)]);

figure; hold on;
plot(repmat([0.5; 1],1,length(ss_params_ori)),[ss_params_ori(:,1)'; ss_params_size(:,1)'], '-s')
errorbar([0.5 1],[mean_amp_ori mean_amp_size],sem_amp,'-sk','LineWidth',3);
ylim([-3 8])
xlim([0.25 1.25])

ax = gca;
ax.XTick = [0.5 1];
ax.XTickLabel = {'Orientation','Size'};

ylabel('DoG Amplitude')

grid on;
grid minor;


% Width
mean_width_ori = mean(ss_params_ori(:,2));
mean_width_size = mean(ss_params_size(:,2));

sem_width = withinstde([ss_params_ori(:,2) ss_params_size(:,2)]);

figure; hold on;
plot(repmat([0.5; 1],1,length(ss_params_ori)),[ss_params_ori(:,2)'; ss_params_size(:,2)'], '-s')
errorbar([0.5 1],[mean_width_ori mean_width_size],sem_width,'-sk','LineWidth',3);
ylim([0.01 0.08])
xlim([0.25 1.25])

ax = gca;
ax.XTick = [0.5 1];
ax.XTickLabel = {'Orientation','Size'};

ylabel('DoG Width')

grid on;
grid minor;

    
[h,p,ci,stats] = ttest(ss_params_ori(:,1),0,'tail','right');
[h,p,ci,stats] = ttest(ss_params_size(:,1),0,'tail','both');

[h,p,ci,stats] = ttest(ss_params_ori(:,1),ss_params_size(:,1));
[h,p,ci,stats] = ttest(ss_params_ori(:,2),ss_params_size(:,2));


%% ===================================================================== %%
%% ==================== Peripheral bump analysis ======================= %%
%% ===================================================================== %%

[ori_bin,size_bin] = meshgrid(40:10:80);
ori_bin = reshape(ori_bin,1,numel(ori_bin));
size_bin = reshape(size_bin,1,numel(size_bin));

for iBin = 1:length(ori_bin)
    
    % specify bins to be evaluated
    pos_bin_ori = ori_bin(iBin):10:80;
    neg_bin_ori = pos_bin_ori * -1;
    
    pos_bin_size = size_bin(iBin):10:80;
    neg_bin_size = pos_bin_size * -1;
    
    % compute the mean response errors in the specified bins
    for iSubject = 1:length(subjects)
        
        pos_bin_bias_ori = mean(sd_data_ori{iSubject}(ismember(sd_data_ori{iSubject}(:,1),pos_bin_ori),2));
        neg_bin_bias_ori = mean(sd_data_ori{iSubject}(ismember(sd_data_ori{iSubject}(:,1),neg_bin_ori),2));
        
        side_bump(iSubject,iBin,1) = (pos_bin_bias_ori - neg_bin_bias_ori) / 2;
               
        pos_bin_bias_size = mean(sd_data_size{iSubject}(ismember(sd_data_size{iSubject}(:,1),pos_bin_size),2));
        neg_bin_bias_size = mean(sd_data_size{iSubject}(ismember(sd_data_size{iSubject}(:,1),neg_bin_size),2));
        
        side_bump(iSubject,iBin,2) = (pos_bin_bias_size - neg_bin_bias_size) / 2;
              
    end
    
end

% one-sample t-tests
[h_ori,p_ori,ci_ori,stats_ori] = ttest(squeeze(side_bump(:,:,1)));
[h_size,p_size,ci_size,stats_size] = ttest(squeeze(side_bump(:,:,2)));

p_ori = reshape(p_ori,5,5);
p_size = reshape(p_size,5,5);

% paired t-test
[h,p,ci,stats] = ttest(squeeze(side_bump(:,:,1)),squeeze(side_bump(:,:,2)));

p = reshape(p,5,5);


tmp_mean_side_bump = squeeze(mean(side_bump,1));

mean_side_bump(:,:,1) = reshape(tmp_mean_side_bump(:,1),5,5);
mean_side_bump(:,:,2) = reshape(tmp_mean_side_bump(:,2),5,5);


%% Plot of side bump summary statistics
figure

subplot(2,3,1);
imagesc(mean_side_bump(:,:,1),[-0.7 0.3])
set(gca,'XTick',[1:5]);
set(gca,'XTickLabel',str2mat('40','50','60','70','80'))
set(gca,'YTick',[1:5]);
set(gca,'YTickLabel',str2mat('40','50','60','70','80'))
ylabel('Size bin')
xlabel('Orientation bin')
title('Avg. error - Orientation condition')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
caxis([-0.7 0.4]);
colorbar

subplot(2,3,2);
imagesc(mean_side_bump(:,:,2),[-0.7 0.3])
set(gca,'XTick',[1:5]);
set(gca,'XTickLabel',str2mat('40','50','60','70','80'))
set(gca,'YTick',[1:5]);
set(gca,'YTickLabel',str2mat('40','50','60','70','80'))
ylabel('Size bin')
xlabel('Orientation bin')
colorbar
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
caxis([-0.7 0.4]);
title('Avg. error - Size condition')

subplot(2,3,3);
figure
imagesc((mean_side_bump(:,:,1) - mean_side_bump(:,:,2)))
set(gca,'XTick',[1:5]);
set(gca,'XTickLabel',str2mat('40','50','60','70','80'))
set(gca,'YTick',[1:5]);
set(gca,'YTickLabel',str2mat('40','50','60','70','80'))
ylabel('Size bin')
xlabel('Orientation bin')
colorbar
caxis([-0.7 0.4]);
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
title('Avg. error - Condition difference')


%% P-values

subplot(2,3,4);
imagesc(p_ori,[0 1])
set(gca,'XTick',[1:5]);
set(gca,'XTickLabel',str2mat('40','50','60','70','80'))
set(gca,'YTick',[1:5]);
set(gca,'YTickLabel',str2mat('40','50','60','70','80'))
ylabel('Size bin')
xlabel('Orientation bin')
colorbar
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
title('p-value - Orientation condition')


subplot(2,3,5);
imagesc(p_size,[0 1])
set(gca,'XTick',[1:5]);
set(gca,'XTickLabel',str2mat('40','50','60','70','80'))
set(gca,'YTick',[1:5]);
set(gca,'YTickLabel',str2mat('40','50','60','70','80'))
ylabel('Size bin')
xlabel('Orientation bin')
colorbar
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
title('p-value - Size condition')


subplot(2,3,6);
imagesc(p,[0 1])
set(gca,'XTick',[1:5]);
set(gca,'XTickLabel',str2mat('40','50','60','70','80'))
set(gca,'YTick',[1:5]);
set(gca,'YTickLabel',str2mat('40','50','60','70','80'))
ylabel('Size bin')
xlabel('Orientation bin')
colorbar
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
title('p-value - Condition difference')


% Bayes factor analysis
for i = 1:25
    bf_ori(i)   = t1smpbf(stats_ori.tstat(i),34);
    bf_size(i)  = t1smpbf(stats_size.tstat(i),34);
    bf_diff(i)  = 1 ./ t1smpbf(stats.tstat(i),34);
    
end

bf_ori  = reshape(bf_ori,5,5);
bf_size = reshape(bf_size,5,5);
bf_diff = reshape(bf_diff,5,5);

figure;

subplot(1,3,1);
imagesc(bf_ori,[0 10])
set(gca,'XTick',[1:5]);
set(gca,'XTickLabel',str2mat('40','50','60','70','80'))
set(gca,'YTick',[1:5]);
set(gca,'YTickLabel',str2mat('40','50','60','70','80'))
ylabel('Size bin')
xlabel('Orientation bin')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
title('BF_{10} - Orientation condition')

subplot(1,3,2);
imagesc(bf_size,[0 10])
set(gca,'XTick',[1:5]);
set(gca,'XTickLabel',str2mat('40','50','60','70','80'))
set(gca,'YTick',[1:5]);
set(gca,'YTickLabel',str2mat('40','50','60','70','80'))
ylabel('Size bin')
xlabel('Orientation bin')
colorbar
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
title('BF_{10} - Size condition')

subplot(1,3,3);
imagesc(bf_diff)
set(gca,'XTick',[1:5]);
set(gca,'XTickLabel',str2mat('40','50','60','70','80'))
set(gca,'YTick',[1:5]);
set(gca,'YTickLabel',str2mat('40','50','60','70','80'))
ylabel('Size bin')
xlabel('Orientation bin')
colorbar
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
title('BF_{01} - Condition difference')