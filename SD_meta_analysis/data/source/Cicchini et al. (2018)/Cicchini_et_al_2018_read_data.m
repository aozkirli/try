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


% subject list from mat files
cd('Experiment2')
ss                  = indir;
for i = 1:numel(ss)
    nm{i}           = ss(i).name(1:10); % extract initials
end
% is kmallrand the same as kmrand?
nm(ismember(nm,'D15OPTkmal')) = {'D15OPTkmra'};
subs                = unique(nm);

% concatenate subjects in a table, check the number of data each
tbl                 = table;
bl                  = 1;
obsidold            = '';
for i = 1:numel(ss)
    load(ss(i).name)
    obsid           = ss(i).name(1:10);
    
    if strcmp(obsid,'D15OPTkmal')
        obsid       = 'D15OPTkmra';
    end
    if strcmp(obsid,obsidold)
        bl          = bl+1;
    else
        bl          = 1;
    end
    obs             = find(ismember(subs,obsid));
    resp            = mod(RESP(:,2),180);
    theta           = mod(RESP(:,3),180);
    error           = sdp_acute_ang(resp-theta);
    delta           = sdp_acute_ang(nbk(theta)-theta);
    obs             = obs.*ones(size(resp));
    block           = bl.*ones(size(resp));
    rt              = RESP(:,6)-RESP(:,5);
    cond            = ones(size(resp));
    tmp             = table(obs,theta,resp,error,delta,block,rt,cond);
    tbl             = vertcat(tbl,tmp);
    obsidold        = obsid;
end

% select the relevant condition
tbl           = tbl_subset(tbl,'cond',1);
tbl_name      = 'Cicchini et al. (2018)';
tbl.study     = repmat({tbl_name},numel(tbl.theta),1);
tbl.experiment= repmat({'Experiment 2'},numel(tbl.theta),1);
tbl.stimulus  = repmat({'Orientation'},numel(tbl.theta),1);

fprintf('Storing the standard table \n\n');

cd('..\..\..\tables\single_experiments')
save Cicchini_et_al_2018.mat tbl tbl_name
writetable(tbl,'Cicchini_et_al_2018.csv','delimiter',';')