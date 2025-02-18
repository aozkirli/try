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

list_subj           = {1:23,[1:9,11:25],[1:7,9:25]};
    
% loop and read+store
for i = 1:3
    %---------------------------------------------------% read experiment i
    fprintf('\nReading Experiment %d\n\n',i);
    cd(currentDir);
    path_i          = [currentDir filesep ['Experiment ' num2str(i)] filesep];
    addpath(genpath(path_i));
    % follow the approach in the shared function a_main_analysis_script.m
    subjects            = list_subj{i};
    cfg                 = [];
    cfg.dataPath        = [currentDir filesep ['Experiment ' num2str(i)] filesep 'Data' filesep];
    cfg.subjects        = subjects; % Subjects to be included in the analysis
    [data, log]         = load_data(cfg);

    fprintf('Converting to standard table format...\n\n');
    % rearrange in table format and clean error
    tbl                 = cell2mat(data);
    tbl                 = array2table(tbl);
    tbl.Properties.VariableNames = {'obs','block','trial','theta','resp','rt','cond'};
    if i ==2
        cond = nan(size(tbl.cond));
        cond(tbl.block == nbk(tbl.block) & ~isnan(tbl.resp) & ~isnan(nbk(tbl.resp))) =...
            tbl.cond(tbl.block == nbk(tbl.block) & ~isnan(tbl.resp) & ~isnan(nbk(tbl.resp))) == ...
            nbk(tbl.cond(tbl.block == nbk(tbl.block) & ~isnan(tbl.resp) & ~isnan(nbk(tbl.resp))));
        tbl.cond = cond;
    end
    % recode
    tbl.theta           = round(tbl.theta+90);
    tbl.resp            = tbl.resp+90;
    obs                 = unique(tbl.obs);
    n                   = numel(obs);
    tbl.delta           = nan(size(tbl.theta));
    for j = 1:n
        j_theta         = tbl.theta(tbl.obs==obs(j));
        j_delta         = sdp_acute_ang(nbk(j_theta)-j_theta);
        tbl.delta(tbl.obs==obs(j)) = j_delta;
    end
    tbl.error           = sdp_acute_ang(tbl.resp-tbl.theta);
   

    tbl_name            = 'Fritsche et al. (2020)';
    tbl.study           = repmat({tbl_name},numel(tbl.theta),1);
    tbl.experiment      = repmat({['Experiment ' num2str(i)]},numel(tbl.theta),1);
    tbl.stimulus        = repmat({'Orientation'},numel(tbl.theta),1);

    fprintf('Storing the standard table \n\n');

    cd(['..' filesep '..' filesep 'experiments'])
    fname               = ['Fritsche_et_al_2020_Exp', num2str(i)];
    save(fname,'tbl','tbl_name')
    writetable(tbl,[fname '.csv'],'delimiter',';')
end


