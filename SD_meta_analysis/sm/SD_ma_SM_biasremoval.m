clc  % Clear command window

% Get the full path of the current script
scriptPath = mfilename('fullpath');

% Extract directory path of the current script
[currentDir, ~, ~] = fileparts(scriptPath);

% Navigate to the current directory
cd(currentDir)

% Add in-house toolbox (must be replaced with the shared folder path containing necessary functions)
addpath(genpath('C:\Users\chare\Google Drive\Work\09_Code\BEIM_toolbox\matlabv\'))
addpath(genpath('/Users/ayberkozkirli/Documents/GitHub/BEIM_toolbox'))

% Change to 'datasets' directory
cd(['..' filesep 'data' filesep 'datasets'])
load('SD_ma_master_table.mat')

%% Plot stimulus-specific bias in orientation and motion
plt                 = [];
plt.x_name          = 'theta_cent';
plt.y_name          = 'errorsd';
plt.nboot           = 10;
plt.space           = 180;

figure('Position',[200 280 1000 500])
subplot(121)
sdp_plot(tbl_subset(tbl,'stimulus','Orientation'),plt);
plt.color           = [1 0 0];
plt.y_name          = 'error_ori_deb_sd_deb';
sdp_plot(tbl_subset(tbl,'stimulus','Orientation'),plt);
ylim([-8 8])
format_figure(0,0,'\Theta(°)','Error(°)');
tl                  = title('Orientation');tl.FontWeight = 'normal';

plt.color           = [0 0 0];
plt.space           = 360;
plt.y_name          = 'errorsd';
subplot(122)
sdp_plot(tbl_subset(tbl,'stimulus','Motion'),plt);
plt.color           = [1 0 0];
plt.y_name          = 'error_ori_deb_sd_deb';
sdp_plot(tbl_subset(tbl,'stimulus','Motion'),plt);
ylim([-8 8])
format_figure(0,0,'\Theta(°)','Error(°)');
tl                  = title('Motion');tl.FontWeight = 'normal';

%% Separate ISO vs. ORTHO
colors              = lines(2);
plt.space           = 180;
plt.color           = colors(1,:);
plt.y_name          = 'errorsd';

figure('Position',[200 280 1000 500])
subplot(121)
[pl1,bt1]           = sdp_plot(tbl_subset(tbl,'stimulus','Orientation','bin',1),plt);
plt.color           = colors(2,:);
[pl2,bt2]           = sdp_plot(tbl_subset(tbl,'stimulus','Orientation','bin',3),plt);
ylim([-8 8])
format_figure(0,0,'\Theta(°)','Error(°)');
tl                  = title('Orientation');tl.FontWeight = 'normal';
format_legend([pl1 pl2],{'iso','ortho'});

std_ori_before      = [bootstrp(1000, @std, bt1.m) bootstrp(1000, @std, bt2.m)];

plt.space           = 360;
plt.color           = colors(1,:);
subplot(122)
[pl1,bt1]           = sdp_plot(tbl_subset(tbl,'stimulus','Motion','bin',1),plt);
plt.color           = colors(2,:);
[pl2,bt2]           = sdp_plot(tbl_subset(tbl,'stimulus','Motion','bin',3),plt);
ylim([-8 8])
format_figure(0,0,'\Theta(°)','Error(°)');
tl                  = title('Motion');tl.FontWeight = 'normal';

std_mot_before      = [bootstrp(1000, @std, bt1.m) bootstrp(1000, @std, bt2.m)];

% Plot also after cleaning
colors              = lines(2);
plt.space           = 180;
plt.color           = colors(1,:);
plt.y_name          = 'error_ori_deb_sd_deb';

figure('Position',[200 280 1000 500])
subplot(121)
[pl1,bt1]           = sdp_plot(tbl_subset(tbl,'stimulus','Orientation','bin',1),plt);
plt.color           = colors(2,:);
[pl2,bt2]           = sdp_plot(tbl_subset(tbl,'stimulus','Orientation','bin',3),plt);
ylim([-8 8])
format_figure(0,0,'\Theta(°)','Error(°)');
tl                  = title('Orientation');tl.FontWeight = 'normal';
format_legend([pl1 pl2],{'iso','ortho'});

std_ori_after       = [bootstrp(1000, @std, bt1.m) bootstrp(1000, @std, bt2.m)];

plt.space           = 360;
plt.color           = colors(1,:);
subplot(122)
[pl1,bt1]           = sdp_plot(tbl_subset(tbl,'stimulus','Motion','bin',1),plt);
plt.color           = colors(2,:);
[pl2,bt2]           = sdp_plot(tbl_subset(tbl,'stimulus','Motion','bin',3),plt);
ylim([-8 8])
format_figure(0,0,'\Theta(°)','Error(°)');
tl                  = title('Motion');tl.FontWeight = 'normal';

std_mot_after       = [bootstrp(1000, @std, bt1.m) bootstrp(1000, @std, bt2.m)];

% Create bar plot with error bars for bootstrap std
figure;
subplot(1, 2, 1) % Orientation
barData = [std_ori_before std_ori_after];
bar(1:4,mean(barData));
hold on;
errorbar(1:4, mean(barData), std(barData), 'k', 'LineStyle', 'none');
set(gca, 'XTickLabel', {'ISO(3SD)', 'ORTHO(3SD)','ISO(brm)','ORTHO(brm)'})
title('Orientation')
xlabel('Condition')
ylabel('Bootstrap Standard Deviation')

subplot(1, 2, 2) % Motion
barData = [std_mot_before std_mot_after];
bar(1:4,mean(barData));
hold on;
errorbar(1:4, mean(barData), std(barData), 'k', 'LineStyle', 'none');
set(gca, 'XTickLabel', {'ISO(3SD)', 'ORTHO(3SD)','ISO(brm)','ORTHO(brm)'})
title('Motion')
xlabel('Condition')
ylabel('Bootstrap Standard Deviation')

%% 2D plots for scatter and bias as a function of theta (using norm errors)
% orientation
tbl_ori             = tbl_subset(tbl,'stimulus','Orientation');
x                   = 0:180;
sldwin              = 20;
scatter_ori         = nan(numel(x),numel(x));
bias_ori            = scatter_ori;
for k = 1:numel(x)
    k
    idx             = mod(x(k)+[-sldwin:sldwin],180);
    tmp             = tbl_ori(ismember(tbl_ori.theta,idx),:);
    mv              = sdp_mvav(tmp.delta,tmp.error_norm);
    scatter_ori(k,:)= mv.std;
    bias_ori(k,:)   = mv.m;
end

%% plot 
figure
subplot(121)
imagesc(x-90,x,scatter_ori-mean(scatter_ori,2));
format_figure(nan,nan,'\Delta(°)','\Theta(°)');
tl                  = title('Scatter'); tl.FontWeight = 'normal';
subplot(122)
imagesc(x-90,x,bias_ori-mean(bias_ori,2))
format_figure(nan,nan,'\Delta(°)','\Theta(°)');
tl                  = title('Bias'); tl.FontWeight = 'normal';

iso_idx             = ismember(abs(x-90),0:10);
ortho_idx           = ismember(abs(x-90),80:90);
sup_ori             = -(mean(scatter_ori(:,iso_idx),2)-mean(scatter_ori(:,ortho_idx),2)); % advantage of iso over ortho, as positive value
% sup_ori             = sup_ori./mean(scatter_ori(:,iso_idx),2);

pos_idx             = ismember(x-90,0:20);
neg_idx             = ismember(x-90,-(0:20));
sdp_ori             = mean(bias_ori(:,pos_idx),2)-mean(bias_ori(:,neg_idx),2);
% sdp_ori             = sdp_ori./mean(bias_ori(:,pos_idx),2);

% Plot scatter and regression line with prediction confidence interval
figure
scatter(sdp_ori, sup_ori, 'k', 'filled', 'MarkerFaceAlpha', 0.2);
xlabel('SD bias');
ylabel('Superiority');
hold on;

% Calculate the model
m                   = fitlm(sdp_ori, sup_ori );

% Get fitted values and prediction intervals
x_fit               = linspace(min(sdp_ori), max(sdp_ori), 100)';
[y_fit, y_ci]       = predict(m, x_fit);

% Plot regression line
plot(x_fit, y_fit, 'k', 'LineWidth', 1.5);

% Plot confidence intervals
plot(x_fit, y_ci(:,1), 'k--', 'LineWidth', 1); % Lower CI
plot(x_fit, y_ci(:,2), 'k--', 'LineWidth', 1); % Upper CI

legend('Data', 'Regression Line', '95% CI', 'Location', 'best');
title('Regression with Prediction Confidence Interval');
hold off;



% %%  Temporary part (problematic, please check...)
% % Create sine and cosine basis sets
% tmp         = tbl_ori(tbl_ori.bin==1 | tbl_ori.bin==3,:);
% nbasis      = 6;
% theta360    = tmp.theta * 2;  % Scale by 2 for orientation data
% theta_fit   = linspace(min(tmp.theta ),max(tmp.theta ),numel(tmp.theta ))'; % for model predictions
% theta360_fit= theta_fit * 2;
% 
% Xsin        = nan(size(theta360, 1), nbasis);
% Xcos        = nan(size(theta360, 1), nbasis);
% Xsin_fit    = Xsin; % for model predictions
% Xcos_fit    = Xsin; % for model predictions
% 
% pred_sin    = [];
% pred_cos    = [];   
% for b = 1:nbasis
%     Xsin(:, b) = sin(deg2rad(theta360 * b));
%     Xcos(:, b) = cos(deg2rad(theta360 * b));
%     pred_sin   = [pred_sin 'Xsin' num2str(b) '*bin' '+'];
%     pred_cos   = [pred_cos 'Xcos' num2str(b) '*bin' '+'];
% 
%     % for model predictions
%     Xsin_fit(:, b) = sin(deg2rad(theta360_fit * b));
%     Xcos_fit(:, b) = cos(deg2rad(theta360_fit * b));
% 
% end
% 
% tmp                 = horzcat(tmp,array2table(Xsin),array2table(Xcos));
% 
% % for model predictions
% Xsin                = Xsin_fit;
% Xcos                = Xcos_fit;
% tmp_iso             = table;
% tmp_iso.bin         = ones(size(Xsin,1),1);
% tmp_iso             = horzcat(tmp_iso,array2table(Xsin),array2table(Xcos));
% 
% tmp_ortho           = table;
% tmp_ortho.bin       = 3*ones(size(Xsin,1),1);
% tmp_ortho           = horzcat(tmp_ortho,array2table(Xsin),array2table(Xcos));
% 
% 
% m                   = fitlm(tmp,['error_norm ~ ' pred_sin pred_cos(1:end-1)]);
% tmp.res             = m.Residuals.Raw;
% 
% 
% plt.space           = 180;
% plt.color           = colors(1,:);
% plt.y_name          = 'error_norm';
% 
% figure('Position',[200 280 1000 500])
% subplot(121)
% [pl1,bt1]           = sdp_plot(tbl_subset(tmp,'bin',1),plt);
% plot(theta_fit-90,predict(m,tmp_iso),'color',plt.color,'linewidth',2);
% plt.color           = colors(2,:);
% [pl2,bt2]           = sdp_plot(tbl_subset(tmp,'bin',3),plt);
% plot(theta_fit-90,predict(m,tmp_ortho),'color',plt.color,'linewidth',2);
% ylim([-1 1])
% format_figure(0,0,'\Theta(°)','Error(z)');
% tl                  = title('Orientation');tl.FontWeight = 'normal';
% format_legend([pl1 pl2],{'iso','ortho'});
% 
% 
% plt.y_name          = 'res';
% plt.color           = colors(1,:);
% 
% subplot(122)
% [pl1,bt1]           = sdp_plot(tbl_subset(tmp,'bin',1),plt);
% plt.color           = colors(2,:);
% [pl2,bt2]           = sdp_plot(tbl_subset(tmp,'bin',3),plt);
% ylim([-1 1])
% format_figure(0,0,'\Theta(°)','Error(z)');
% tl                  = title('Orientation');tl.FontWeight = 'normal';
% format_legend([pl1 pl2],{'iso','ortho'});
% 
% %%
% % Select only numeric columns from the table 'tmp'
% numericColumns = varfun(@isnumeric, tmp, 'OutputFormat', 'uniform'); 
% 
% % Keep only the numeric columns
% tmp = tmp(:, numericColumns);
% 
% scatter_tbl_cleaned = grpstats(tmp,{'bin','studynum','obsid','codenum'},'std');
% 
% fitlm(scatter_tbl_cleaned,'std_res ~ bin') % std_error_ori_deb_sd_deb
% 
