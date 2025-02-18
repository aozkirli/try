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


fprintf('\nReading Experiment %d, control group\n\n',1);
cd(currentDir);
tbl             = readtable('SD and SZ-dataset-zenodo.csv');
% exclude already excluded subjects and select only the control group
tbl             = tbl_subset(tbl,'group',1,'exclude',0);
fprintf('Converting to standard table format...\n\n');
tmp             = table;
tmp.obs         = tbl.obs;
tmp.theta       = tbl.theta;
tmp.resp        = tbl.resp;
tmp.error       = tbl.error;
tmp.rt          = tbl.rt;
tmp.block       = tbl.session;
tmp.delta       = tbl.delta; % target delta
tmp.trial       = (1:numel(tbl.delta))';
tmp.cond        = ones(size(tmp.theta));
tbl             = tmp;
tbl_name            = 'Pascucci et al. (2024)';
tbl.study           = repmat({tbl_name},numel(tbl.theta),1);
tbl.experiment      = repmat({['Experiment ' num2str(1)]},numel(tbl.theta),1);
tbl.stimulus        = repmat({'Orientation'},numel(tbl.theta),1);
fprintf('Storing the standard table \n\n');

cd('..\..\tables\single_experiments')
fname         = 'Pascucci_et_al_2024';
save(fname,'tbl','tbl_name')
writetable(tbl,[fname '.csv'],'delimiter',';')
