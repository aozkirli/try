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
addpath(genpath('/Users/ayberkozkirli/Documents/GitHub/BEIM_toolbox'))

exp_folds           = {'MotionCloud', 'RDK', 'RDK_CursorReport'};
exp_id              = {'MC_subj', 'RDK_subj', 'RDK_subj'};

% loop and read+store
for j = 1:3
    %---------------------------------------------------% read experiment j
    fprintf('\nReading Experiment %d\n\n',j);
    path_i              = [currentDir filesep exp_folds{j} filesep];
    cd(path_i);
    
    nblocks             = 15;
    nsubj               = 8;
    tbl                 = table;

    for i = 1:nsubj
        for k = 1:nblocks
            load([exp_id{j} num2str(i) '_b' num2str(k) '.mat']);
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

    tbl_name            = 'Moon et al. (2023)';
    tbl.study           = repmat({tbl_name},numel(tbl.theta),1);
    tbl.experiment      = repmat({['Experiment ' num2str(j)]},numel(tbl.theta),1);
    tbl.stimulus        = repmat({'Motion'},numel(tbl.theta),1);

    fprintf('Storing the standard table \n\n');

    % cd(['..' filesep '..' filesep '..' filesep 'tables' filesep 'single_experiments'])
    cd(['..' filesep '..' filesep '..' filesep 'experiments'])
    fname         = ['Moon_et_al_2023_Exp', num2str(j)];
    save(fname,'tbl','tbl_name')
    writetable(tbl,[fname '.csv'],'delimiter',';')
end


