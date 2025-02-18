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

%% Plot the scatter with 3SD or the orideb cleaning
study_list          = unique(tbl.study);
nstudies            = numel(study_list);

figure
for i = 1:nstudies
    tmp             = tbl_subset(tbl,'study',study_list{i});
    tmp(abs(tmp.delta)>90,:) = [];
    subplot(5,5,i);    
    mv              = sdp_mvav(tmp.delta,tmp.errorsd);
    plot(mv.x,mv.std-mean(mv.std),'k-','linewidth',2)
    hold on;
    mv              = sdp_mvav(tmp.delta,tmp.error_ori_deb);
    plot(mv.x,mv.std-mean(mv.std),'r-','linewidth',2)
    xlim([-90 90]);ylim([-1 1].*max(abs(ylim)));
    format_figure(0,0,'\Delta(°)',{'Scatter','(centered)'},[],[],5);
    tl              = title(study_list{i}); tl.FontWeight = 'normal';
    drawnow
end
sgtitle('Scatter')
legend('3SD','Stim-bias removed')

figure % Bias
for i = 1:nstudies
    tmp             = tbl_subset(tbl,'study',study_list{i});
    tmp(abs(tmp.delta)>90,:) = [];
    subplot(5,5,i);    
    mv              = sdp_mvav(tmp.delta,tmp.errorsd);
    plot(mv.x,mv.m-mean(mv.m),'k-','linewidth',2)
    hold on;
    mv              = sdp_mvav(tmp.delta,tmp.error_ori_deb);
    plot(mv.x,mv.m-mean(mv.m),'r-','linewidth',2)
    xlim([-90 90]);ylim([-1 1].*max(abs(ylim)));
    format_figure(0,0,'\Delta(°)',{'Error(°)'},[],[],5);
    tl              = title(study_list{i}); tl.FontWeight = 'normal';
    drawnow
end
sgtitle('Serial Dependence Bias')
legend('3SD','Stim-bias removed')

