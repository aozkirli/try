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

% loop and read+store
for j = 1:2
    %---------------------------------------------------% read experiment j
    fprintf('\nReading Experiment %d\n\n',j);
    cd(currentDir);
    load(['EC_SD_Unc_Exp' num2str(j) '.mat']);
    
    fprintf('Converting to standard table format...\n\n');
    tmp             = table;
    
    tmp.obs         = tbl.obs;
    tmp.theta       = tbl.theta;
    tmp.resp        = tbl.resp;
    tmp.error       = tbl.error;
    tmp.rt          = tbl.rt;
    tmp.block       = tbl.block;
    tmp.delta       = tbl.delta; % target delta
    tmp.trial       = tbl.trial;
    tmp.cond        = double(tbl.sd==max(tbl.sd));
    
    tbl             = tmp;

    tbl_name            = 'Ozkirli & Pascucci (2023)';
    tbl.study           = repmat({tbl_name},numel(tbl.theta),1);
    tbl.experiment      = repmat({['Experiment ' num2str(j)]},numel(tbl.theta),1);
    tbl.stimulus        = repmat({'Orientation'},numel(tbl.theta),1);

    fprintf('Storing the standard table \n\n');

    cd('..\..\tables\single_experiments')
    fname         = ['Ozkirli_Pascucci_2023_Exp', num2str(j)];
    save(fname,'tbl','tbl_name')
    writetable(tbl,[fname '.csv'],'delimiter',';')
end

%% third experiment with feedback

fprintf('\nReading Experiment %d\n\n',j+1);
cd(currentDir);
load(['EC_SD_Unc_Exp' num2str(j+1) '.mat']);
fprintf('Converting to standard table format...\n\n');
tmp             = table;
tmp.obs         = tbl.obs;
tmp.theta       = tbl.theta;
tmp.resp        = tbl.resp;
tmp.error       = tbl.error;
tmp.rt          = tbl.rt;
tmp.block       = tbl.block;
tmp.delta       = tbl.delta; % target delta
tmp.trial       = tbl.trial;
tmp.cond        = double(tbl.feedback==max(tbl.feedback));
tbl             = tmp;
tbl_name            = 'Ozkirli & Pascucci (2023)';
tbl.study           = repmat({tbl_name},numel(tbl.theta),1);
tbl.experiment      = repmat({['Experiment ' num2str(j+1)]},numel(tbl.theta),1);
tbl.stimulus        = repmat({'Orientation'},numel(tbl.theta),1);
fprintf('Storing the standard table \n\n');

cd('..\..\tables\single_experiments')
fname         = ['Ozkirli_Pascucci_2023_Exp', num2str(j+1)];
save(fname,'tbl','tbl_name')
writetable(tbl,[fname '.csv'],'delimiter',';')
