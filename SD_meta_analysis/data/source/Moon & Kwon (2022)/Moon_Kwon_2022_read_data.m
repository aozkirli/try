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

cd('221012')
nblocks             = 15;
nsubj               = 32;
tbl                 = table;
for i = 1:nsubj
    for k = 1:nblocks
        load(['E_subj' num2str(i) '_b' num2str(k) '.mat']);
        tmp         = table;        
        tmp.delta   = array(:,1);
        tmp.theta   = array(:,2);
        tmp.resp    = array(:,3);
        tmp.error   = array(:,4);
        tmp.rt      = array(:,5);
        tmp.obs     = i.*ones(size(tmp.delta));
        tmp.block   = k.*ones(size(tmp.delta));
        tmp.range   = range(tmp.delta).*ones(size(tmp.delta));
        tmp.mindelta= min(tmp.delta).*ones(size(tmp.delta));
        tbl         = vertcat(tbl,tmp);
    end
end

tbl                 = tbl_subset(tbl,'mindelta',-180);

tbl_name            = 'Moon & Kwon (2022)';
tbl.study           = repmat({tbl_name},numel(tbl.theta),1);
tbl.experiment      = repmat({'Experiment 1'},numel(tbl.theta),1);
tbl.stimulus        = repmat({'Motion'},numel(tbl.theta),1);

fprintf('Storing the standard table \n\n');

cd('..\..\..\tables\single_experiments')
save Moon_Kwon_2022.mat tbl tbl_name
writetable(tbl,'Moon_Kwon_2022.csv','delimiter',';')

