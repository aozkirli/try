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

ExpNames        = {'ExperimentColorCue_Data.txt',...
                   'ExperimentSerialPositionCue_Data.txt',...
                   'ExperimentSimColCue_Data.txt',...
                   'ExperimentSimPosCue_Data.txt'};

% loop and read+store
for i = 1:4
    cd(currentDir);
    %---------------------------------------------------% read experiment j
    fprintf('\nReading Experiment %d\n\n',i);
    Expdata             = readtable(ExpNames{i});
    obs                 = table2array( Expdata(:,1));
    rt                  = table2array( Expdata(:,17));
    target              = table2array( Expdata(:,12));
    resp                = table2array( Expdata(:,14));
    one                 = table2array( Expdata(:,4));
    two                 = table2array( Expdata(:,5));
    same                = table2array( Expdata(:,end));
    stim                = [one two];
    tar                 = nan(numel(target),1); % the target direction
    for k = 1:numel(target)
        tar(k,:)        = stim(k,target(k));
    %     nontar(k,:)     = stim(k,~ismember(1:2,target(k)));
    end
    % compute errors and delta
    error               = sdp_acute_ang(resp-tar,360);
    delta               = nan(size(error));
    for j = 1:max(obs)
        tmp             = tar(obs==j);
        dd              = sdp_acute_ang(nbk(tmp)-tmp,360);
        delta(obs==j,:) = dd;
    end
    theta               = tar;
    tbl                 = table(obs,delta,resp,error,theta,rt);
       
    fprintf('Converting to standard table format...\n\n');

    tbl_name      = 'Fischer et al. (2020)';
    tbl.study     = repmat({tbl_name},numel(tbl.theta),1);
    tbl.experiment= repmat({['Experiment ' num2str(i)]},numel(tbl.theta),1);
    tbl.stimulus  = repmat({'Motion'},numel(tbl.theta),1);

    fprintf('Storing the standard table \n\n');

    cd('..\..\tables\single_experiments')
    fname         = ['Fischer_et_al_2020_Exp', num2str(i)];
    save(fname,'tbl','tbl_name')
    writetable(tbl,[fname '.csv'],'delimiter',';')
end


