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


%-------------------------------------------------------% read experiment 1
fprintf('Reading Experiment 1...excluding the mixed block\n\n');
load('Experiment_01_dataset.mat')
tbl             = tabledata;
% remove string variables
tbl(:,[2 12])   = [];
tbl             = tbl_subset(tbl,'block','~=3'); % keep only blocked conditions (exclude the mixed block)


fprintf('Converting to standard table format...\n\n');

tbl_name      = 'Ceylan et al. (2021)';
tbl.study     = repmat({tbl_name},numel(tbl.theta),1);
tbl.experiment= repmat({'Experiment 1'},numel(tbl.theta),1);
tbl.stimulus  = repmat({'Orientation'},numel(tbl.theta),1);
tbl.cond      = double(tbl.stimtype==2);

fprintf('Storing the standard table \n\n');

cd('..\..\tables\single_experiments')
save Ceylan_et_al_2021_Exp1.mat tbl tbl_name
writetable(tbl,'Ceylan_et_al_2021_Exp1.csv','delimiter',';')

%-------------------------------------------------------% read experiment 2
cd(currentDir)
fprintf('Reading Experiment 2...excluding the mixed block\n\n');
load('Experiment_02_dataset.mat')
tbl             = tabledata;
% remove string variables
tbl(:,[2 12])   = [];
tbl             = tbl_subset(tbl,'block','~=3'); % keep only blocked conditions (exclude the mixed block)

fprintf('Converting to standard table format...\n\n');

tbl_name      = 'Ceylan et al. (2021)';
tbl.study     = repmat({tbl_name},numel(tbl.theta),1);
tbl.experiment= repmat({'Experiment 2'},numel(tbl.theta),1);
tbl.stimulus  = repmat({'Orientation'},numel(tbl.theta),1);
tbl.cond      = double(tbl.stimtype==2);

fprintf('Storing the standard table \n\n');

cd('..\..\tables\single_experiments')
save Ceylan_et_al_2021_Exp2.mat tbl tbl_name
writetable(tbl,'Ceylan_et_al_2021_Exp2.csv','delimiter',';')
