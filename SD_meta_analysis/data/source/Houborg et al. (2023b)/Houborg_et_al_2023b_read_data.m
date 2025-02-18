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
for i = 1:2
    %---------------------------------------------------% read experiment i
    fprintf('\nReading Experiment %d\n\n',i);
    cd(currentDir);
    cd('Datasets')
    load(['exp' num2str(i) '_data.mat']);
    % remove outlier subjects as in the main paper
    if i == 1
    tbl             = tbl_subset(exp1,'propout','<0.25','ccorr','>.4','macc','>.6');
    else
    tbl             = tbl_subset(exp2,'propout','<0.25','ccorr','>.4','macc','>.6');
    end
    
    fprintf('Converting to standard table format...\n\n');
    % rearrange in table format and clean error
    tmp             = table;
    tmp.obs         = tbl.obs;
    tmp.theta       = tbl.theta;
    tmp.resp        = tbl.resp;
    tmp.error       = tbl.error;
    tmp.rt          = tbl.rt;
    tmp.block       = tbl.block;
    tmp.delta       = tbl.delta;
    tmp.trial       = tbl.trial;

    tbl             = tmp;
    tbl_name        = 'Houborg et al. (2023b)';
    tbl.study       = repmat({tbl_name},numel(tbl.theta),1);
    tbl.experiment  = repmat({['Experiment ' num2str(i)]},numel(tbl.theta),1);
    tbl.stimulus    = repmat({'Orientation'},numel(tbl.theta),1);

    fprintf('Storing the standard table \n\n');

    cd('..\..\..\tables\single_experiments')
    fname         = ['Houborg_et_al_2023b_Exp', num2str(i)];
    save(fname,'tbl','tbl_name')
    writetable(tbl,[fname '.csv'],'delimiter',';')
end


