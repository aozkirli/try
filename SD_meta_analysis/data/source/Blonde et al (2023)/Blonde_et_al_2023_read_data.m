clc
% full path of the current script
scriptPath          = mfilename('fullpath');

% directory path
[currentDir, ~, ~]  = fileparts(scriptPath);

cd(currentDir)
fileID              = fopen('source.txt', 'r');

source              = fread(fileID, '*char')';
fclose(fileID);

fprintf(source)


% add inhouse toolbox
% this path must be replaced with the shared folder containing the
% collection of functions required
addpath(genpath('C:\Users\chare\Google Drive\Work\09_Code\BEIM_toolbox\matlabv\'))


% % read experiment 1 % EXCLUDED BECAUSE no 90 deg
% tbl             = readtable('Datatable_EXP1.csv');
% 
% % check and remove outlier subjects (as in the manuscript)
% tbl             = tbl_subset(tbl,'ptrialout','<.2','errorsigm','<45');
% actual_n        = numel(unique(tbl.obs));
% 
% fprintf('Converting to standard table format...\n\n');
% 
% tmp             = table;
% tmp.obs         = tbl.obs;
% tmp.theta       = tbl.theta;
% tmp.resp        = tbl.resp;
% tmp.error       = tbl.error;
% tmp.rt          = tbl.rt;
% tmp.block       = tbl.block;
% tmp.delta       = tbl.delta_btw; % the delta of interest: previously reported orientations
% tmp.trial       = tbl.trial;
% % tmp.cond        = tbl.condition; % the condition manipulation did not
% % affect delta_btw effects
% tbl             = tmp;
% 
% tbl_name      = 'Blonde et al. (2023)';
% tbl.study     = repmat({tbl_name},numel(tbl.theta),1);
% tbl.experiment= repmat({'Experiment 1'},numel(tbl.theta),1);
% tbl.stimulus  = repmat({'Orientation'},numel(tbl.theta),1);
% 
% fprintf('Storing the standard table \n\n');
% 
% cd('..\..\tables\single_experiments')
% save Blonde_et_al_2023_Exp1.mat tbl tbl_name
% writetable(tbl,'Blonde_et_al_2023_Exp1.csv','delimiter',';')


cd(currentDir)

% read experiment 2
tbl             = readtable('Datatable_EXP2.csv');

% keep only the uniform condition
tbl             = tbl_subset(tbl,'condition',0);

% check and remove outlier subjects (as in the manuscript)
tbl             = tbl_subset(tbl,'ptrialout','<.2','errorsigm','<45');
actual_n        = numel(unique(tbl.obs));

fprintf('Converting to standard table format...\n\n');

tmp             = table;
tmp.obs         = tbl.obs;
tmp.theta       = tbl.theta;
tmp.resp        = tbl.resp;
tmp.error       = tbl.error;
tmp.rt          = tbl.rt;
tmp.block       = tbl.block;
tmp.delta       = tbl.delta_btw; % the delta of interest: previously reported orientations
tmp.trial       = tbl.trial;
% tmp.cond        = tbl.condition;
tbl             = tmp;

tbl_name      = 'Blonde et al. (2023)';
tbl.study     = repmat({tbl_name},numel(tbl.theta),1);
tbl.experiment= repmat({'Experiment 2'},numel(tbl.theta),1);
tbl.stimulus  = repmat({'Orientation'},numel(tbl.theta),1);

fprintf('Storing the standard table \n\n');

cd('..\..\tables\single_experiments')
save Blonde_et_al_2023_Exp2.mat tbl tbl_name
writetable(tbl,'Blonde_et_al_2023_Exp2.csv','delimiter',';')
