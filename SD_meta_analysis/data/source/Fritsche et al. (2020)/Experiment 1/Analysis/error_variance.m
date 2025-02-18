clear all; close all;
addpath /home/common/matlab/fieldtrip/qsub/
addpath CircStat2012a/

%% Main analysis of the serial dependence experiment

%% Main parameters of the analysis
subjects            = 1:23;
remove_outliers     = true;

%% =================== load the raw data ================================ %
cfg                 = [];
cfg.dataPath        = '../Data/';
cfg.subjects        = subjects; % Subjects to be included in the general analysis

[data, log]         = load_data(cfg);

%% ===================== Outlier removal   ============================== %
if remove_outliers
    cfg                     = [];
    cfg.data                = data;
    cfg.subjects            = subjects;
    
    [data, data_removed]    = outlier_removal(cfg);
end

%% ===================== create SD sorted data ========================== %
cfg         = [];
cfg.data    = data;
cfg.nBack   = 20;
cfg.demean  = true;

sd_data = sort_data(cfg);


%% Compute moving averages of sdev
cfg             = [];
cfg.data        = sd_data;
cfg.bin_width   = 30; % possible to specify multiple window sizes
cfg.demean      = true;
        
[moving_averages] = compute_moving_averages_std(cfg);


%% Compute n-back means and SEMs
for iBack = 1:20
    
    % all (pooled over same and different location trials)
    thisdata_all = cell2mat(moving_averages.all{iBack});  
    grand_mean_all(iBack,:) = mean(thisdata_all);
    grand_sem_all(iBack,:)  = std(thisdata_all)/sqrt(length(subjects));
      
end


%% Fit second derivative of a gaussian to mean variance profile
coeffs = zeros(20,3);

for iBack = 1:20
    
    disp(['Fitting ' num2str(iBack) '-back variance...'])
   
    cfg                 = [];
    cfg.data            = grand_mean_all(iBack,:)';
    cfg.fixedwidth      = false;
    cfg.fittingsteps    = 1;
    thisfit             = fit_dog2(cfg);
    
    coeffs(iBack,:) = thisfit.coeffs;
    
end

save('Matlab datafiles/group fits/variance_demeaned.mat','coeffs');

%% Plot variances with model fit
dog2_fun = @(a,b,w,x) a * (2*w.^2.*x.^2-1).*exp(-(w.*x).^2) + b;


figure;
for iBack = 1:10
    subplot(2,5,iBack); hold on;
    
    shadedErrorBar(-90:90,grand_mean_all(iBack,:),grand_sem_all(iBack,:));
    
    
    plot([-90:0.01:90], feval(dog2_fun, ...
        coeffs(iBack,1), coeffs(iBack,2), coeffs(iBack,3), ...
        [-90:0.01:90]), '-k','LineWidth',3);
    
    
    title([num2str(iBack) "-back"])
    xlim([-90 90])
    %ylim([9.5 13])
    ylim([-1.5 1.5])
    xlabel('Rel. n-back ori')
    ylabel('SD of response distribution')
    set(gca,'FontSize',14);
end
suptitle("Exp. 1 - All trials")



%% Bootstrap variance modulation

% precompile permutation function for cluster execution
compiled_bootstrap  = qsubcompile('bootstrap_dog2','toolbox', {'curvefit'});
time    = 10 * 3600; % hours to seconds
space   = 1 * 1024^3; % Gb to bytes

for iBack = 1:20
   
    disp(['Bootstraping ' num2str(iBack) '-back variance...'])
    
    cfg                     = [];
    cfg.subjects            = subjects;
    cfg.data                = moving_averages.all{iBack};
    cfg.nBootstrapSamples   = 10000;
    cfg.fittingsteps        = 1;
    cfg.fixedwidth          = false;
    cfg.savePath            = ['Matlab datafiles/bootstrap/variance_bootstrap_' num2str(iBack) '-back_demeaned.mat'];
    
    %variance_bootstrap_params{iBack,1} = bootstrap_dog2(cfg);
    jobidarray{iBack} = qsubfeval(compiled_bootstrap,cfg,'timreq',time,'memreq',space); 
 
end


%% ================= Bar plot of variance estimates   =================== %
load('Matlab datafiles/group fits/variance_demeaned.mat','coeffs');
amplitudes = coeffs(:,1) * -1;

nback = 20;


% load bootstrap estimates and determine CIs and SDs
for iBack = 1:nback
    
    load(['Matlab datafiles/bootstrap/variance_bootstrap_' num2str(iBack) '-back_demeaned.mat'])
    bootstrap_amplitudes = sort(bootstrap_params(:,1)*-1);
        
    bootstrap_sd(iBack,:) = std(bootstrap_amplitudes);
    
end


figure; hold on;
bar(1:nback,amplitudes(1:nback));
errorbar(1:nback,amplitudes(1:nback),bootstrap_sd,'.k','CapSize',0)
ylim([-1.5 1])
xticks(1:nback)


%% =============== Individual n-back variance plots   =====================
dog2_fun = @(a,b,w,x) a * (2*w.^2.*x.^2-1).*exp(-(w.*x).^2) + b;

iBack = 1;

figure; hold on;

shadedErrorBar(-90:90,grand_mean_all(iBack,:),grand_sem_all(iBack,:));

plot([-90:0.01:90], feval(dog2_fun, ...
    coeffs(iBack,1), coeffs(iBack,2), coeffs(iBack,3), ...
    [-90:0.01:90]), '-k','LineWidth',3);

xlim([-90 90])
ylim([-1.5 1.5])




