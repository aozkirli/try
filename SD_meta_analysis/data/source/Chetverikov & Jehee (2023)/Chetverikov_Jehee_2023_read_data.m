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

fprintf('Loading experiment 1 (low coherence)...\n\n');

% read experiment 1
tmp                 = readtable('low_coherence.csv');


fprintf('Converting to standard table format...\n\n');

% recode
tbl                 = table;
tbl.obs             = tmp.subject_id;
tbl.theta           = tmp.Orientation+180;
tbl.resp            = tmp.Response+180;
tbl.error           = tmp.error;
tbl.block           = tmp.run;
tbl.cond            = ones(size(tbl.obs));
tbl.rt              = 9*ones(size(tmp.subject_id)); %% missing RT just put 1
tbl.trials          = tmp.trialN;

% add variables 
nobs            = unique(tbl.obs);
n               = numel(nobs);
retbl           = table();
for i = 1:n
    i_tbl       = tbl_subset(tbl,'obs',nobs(i));
    i_tbl.delta = sdp_acute_ang(nbk(i_tbl.theta)-i_tbl.theta,360);
    retbl       = vertcat(retbl,i_tbl);
end
tbl             = retbl;

tbl_name        = 'Chetverikov & Jehee (2023)';
tbl.study       = repmat({tbl_name},numel(tbl.theta),1);
tbl.experiment  = repmat({'Experiment 1'},numel(tbl.theta),1);
tbl.stimulus    = repmat({'Motion'},numel(tbl.theta),1);

% fprintf('Storing the standard table \n\n');

tbl1            = tbl; % to combine later

% cd('..\..\tables\single_experiments')
% save Chetverikov_Jehee_2023_Exp1.mat tbl tbl_name
% writetable(tbl,'Chetverikov_Jehee_2023_Exp1.csv','delimiter',';')


%--------------------------------------------------------------------------
fprintf('Loading experiment 2 (high coherence)...\n\n');

% read experiment 2
tmp                 = readtable('high_coherence.csv');


fprintf('Converting to standard table format...\n\n');

% recode
tbl                 = table;
tbl.obs             = str2double(strrep(strrep(tmp.true_sid, 'S', ''), 'P', ''));
tbl.theta           = tmp.mdir;
tbl.resp            = tmp.resp;
tbl.error           = sdp_acute_ang(tbl.resp-tbl.theta,360);
tbl.block           = tmp.run;
tbl.cond            = 2*ones(size(tbl.obs)); % high coherence
tbl.rt              = 9*ones(size(tmp.resp)); %% missing RT just put 1
tbl.trials          = tmp.trial;

% add variables 
nobs            = unique(tbl.obs);
n               = numel(nobs);
retbl           = table();
for i = 1:n
    i_tbl       = tbl_subset(tbl,'obs',nobs(i));
    i_tbl.delta = sdp_acute_ang(nbk(i_tbl.theta)-i_tbl.theta,360);
    retbl       = vertcat(retbl,i_tbl);
end
tbl             = retbl;

tbl_name        = 'Chetverikov & Jehee (2023)';
tbl.study       = repmat({tbl_name},numel(tbl.theta),1);
tbl.experiment  = repmat({'Experiment 1'},numel(tbl.theta),1);
tbl.stimulus    = repmat({'Motion'},numel(tbl.theta),1);

tbl             = vertcat(tbl1,tbl);

fprintf('Storing the standard table \n\n');

cd('..\..\tables\single_experiments')
save Chetverikov_Jehee_2023.mat tbl tbl_name
writetable(tbl,'Chetverikov_Jehee_2023.csv','delimiter',';')
