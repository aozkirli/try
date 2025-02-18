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
tmp                 = readtable('all_adjustment_task_data.csv');


fprintf('Converting to standard table format...\n\n');

% recode
tbl                 = table;
tbl.obs             = tmp.ParticipantNo_;
tbl.theta           = tmp.StimulusOrientation;
tbl.resp            = tmp.ParticipantResponse;
tbl.cond            = tmp.NoiseLevel_ConcentrationOfVonMises_; % noise
tbl.rt              = 9*ones(size(tmp.ParticipantNo_)); %% missing RT just put 1
% add variables 
nobs            = unique(tbl.obs);
n               = numel(nobs);
retbl           = table();
for i = 1:n
    i_tbl       = tbl_subset(tbl,'obs',nobs(i));
    i_tbl.delta = sdp_acute_ang(nbk(i_tbl.theta)-i_tbl.theta);
    i_tbl.error = sdp_acute_ang(i_tbl.resp-i_tbl.theta);
    retbl       = vertcat(retbl,i_tbl);
end
tbl             = retbl;

tbl_name        = 'Gallagher & Benton (2022)';
tbl.study       = repmat({tbl_name},numel(tbl.theta),1);
tbl.experiment  = repmat({'Experiment 1'},numel(tbl.theta),1);
tbl.stimulus    = repmat({'Orientation'},numel(tbl.theta),1);

fprintf('Storing the standard table \n\n');

cd('..\..\tables\single_experiments')
save Gallagher_Benton_2022.mat tbl tbl_name
writetable(tbl,'Gallagher_Benton_2022.csv','delimiter',';')
