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

cd('data')
folders         = indir();

for k = 1:numel(folders)
    tbl             = table();

    fprintf(['Reading ' folders(k).name ' data...\n\n']);
    cd(currentDir)
    cd('data')
    cd(folders(k).name)
    ss          = indir();
    % there are two sessions per subject, split and give the same subject
    % id
    for j = 1:2
        ss_j    = ss(j:2:end);
        for i = 1:numel(ss_j)
            load(ss_j(i).name)
            obs     = i.*ones(size(RTs'));
            rt      = RTs';
            theta   = mod(TargetOrientations',180);
            resp    = mod(ResponseOrientations',180);
            block   = repmat(1:10,numel(RTs)/10,1); block = block(:);
            % add variables
            delta   = sdp_acute_ang(nbk(theta)-theta);
            error   = sdp_acute_ang(resp-theta);
            cond    = ones(size(delta));
            tbl     = vertcat(tbl,table(obs,theta,resp,delta,error,cond,block,rt));
        end
    end

    fprintf('Converting to standard table format...\n\n');
    
    tbl_name      = 'Kondo et al. (2022)';
    tbl.study     = repmat({tbl_name},numel(tbl.theta),1);
    tbl.experiment= repmat({['Experiment ' num2str(k)]},numel(tbl.theta),1);
    tbl.stimulus  = repmat({'Orientation'},numel(tbl.theta),1);
        
    fprintf('Storing the standard table \n\n');
    
    cd(currentDir)
    cd('..\..\tables\single_experiments')
    fname         = ['Kondo_et_al_2022_Exp', num2str(k)];
    save(fname,'tbl','tbl_name')
    writetable(tbl,[fname '.csv'],'delimiter',';')
end