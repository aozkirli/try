addpath /home/common/matlab/fieldtrip/qsub/
addpath /home/predatt/matfri/Toolboxes/CircStat2012a

clear all; close all;

subjects = [1:7,9:25];

%% Load data
cfg = [];
cfg.subjects = subjects;

[data, log] = load_data(cfg);


%% Outlier removal
cfg                     = [];
cfg.data                = data;
cfg.subjects            = subjects;
cfg.plot                = false;

[data, data_removed]    = outlier_removal(cfg);


%% Sort SD data
cfg             = [];
cfg.subjects    = subjects;
cfg.data        = data;
cfg.demean      = true; % demean response errors to remove overall bias
cfg.nBack       = 20;

[sd_data] = sort_data(cfg);


%% Compute moving averages
cfg             = [];
cfg.data        = sd_data;
cfg.bin_width   = 30;% possible to specify multiple window sizes
cfg.demean      = true;

moving_averages = compute_moving_averages_variance(cfg);


%% Compute n-back means and variances
for iBack = 1:20
    
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


%% Plot variances
dog2_fun = @(a,b,w,x) a * (2*w.^2.*x.^2-1).*exp(-(w.*x).^2) + b;

figure;
for iBack = 1:10
    subplot(2,5,iBack); hold on;
    shadedErrorBar(-90:90,grand_mean_all(iBack,:),grand_sem_all(iBack,:));
    plot([-90:0.01:90], feval(dog2_fun, ...
        coeffs(iBack,1), coeffs(iBack,2), coeffs(iBack,3), ...
        [-90:0.01:90]), '-k','LineWidth',3);
    title([num2str(iBack) "-back"])
    ylim([-1 1])
    xlim([-90 90])
    xlabel('Rel. n-back ori')
    ylabel("SD of response distribution")
    set(gca,'FontSize',14);
end
suptitle("Exp. 3")


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
ylim([-1 0.5])
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
