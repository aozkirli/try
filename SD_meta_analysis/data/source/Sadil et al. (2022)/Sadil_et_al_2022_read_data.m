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


% read experiment 
cd(currentDir)
cd('events')

fprintf('\nConverting to standard table format...\n\n');

subject_list        = ls('*_data.csv');
% create table
tbl                 = table();
n                   = size(subject_list,1);
for i = 1:n
    tmp             = readtable(subject_list(i,:));
    theta           = mod(tmp.orientation,180);
    resp            = mod(tmp.response,180);
    delta           = sdp_acute_ang(nbk(theta)-theta);
    error           = sdp_acute_ang(resp-theta);
    obs             = tmp.subject;
    trial           = tmp.trial;
    block           = tmp.block;
    rt              = tmp.rt;
    i_tbl           = table(obs,theta,delta,trial,block,rt,error,resp);    
    tbl             = vertcat(tbl,i_tbl);
end

tbl_name        = 'Sadil et al. (2022)';
tbl.study       = repmat({tbl_name},numel(tbl.theta),1);
tbl.experiment  = repmat({'Experiment 1'},numel(tbl.theta),1);
tbl.stimulus    = repmat({'Orientation'},numel(tbl.theta),1);

fprintf('Storing the standard table \n\n');

cd('..\..\..\tables\single_experiments')
save Sadil_et_al_2022.mat tbl tbl_name
writetable(tbl,'Sadil_et_al_2022.csv','delimiter',';')




