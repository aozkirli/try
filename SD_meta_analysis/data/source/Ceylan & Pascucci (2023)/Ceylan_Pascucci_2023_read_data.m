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


% read experiment 2
tbl             = readtable('Ceylan_Pascucci_Exp2_2023.csv');


fprintf('Converting to standard table format...\n\n');

tmp             = table;
tmp.obs         = tbl.obs;
tmp.theta       = tbl.target;
tmp.resp        = tbl.resp;
tmp.error       = tbl.error;
tmp.rt          = tbl.rt;
tmp.block       = tbl.block;
tmp.delta       = tbl.delta_targ; % relevant delta (previous target inducing attractive SD)
tmp.trial       = tbl.nTrial;
tbl             = tmp;


tbl_name      = 'Ceylan & Pascucci (2023)';
tbl.study     = repmat({tbl_name},numel(tbl.theta),1);
tbl.experiment= repmat({'Experiment 2'},numel(tbl.theta),1);
tbl.stimulus  = repmat({'Orientation'},numel(tbl.theta),1);

fprintf('Storing the standard table \n\n');

cd('..\..\tables\single_experiments')
save Ceylan_Pascucci_2023_Exp2.mat tbl tbl_name
writetable(tbl,'Ceylan_Pascucci_2023_Exp2.csv','delimiter',';')