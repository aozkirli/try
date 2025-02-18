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
cd('Data')
ss                  = indir('*mat');
tbl                 = table();

fprintf('\nConverting to standard table format...\n\n');

for i = 1:numel(ss)
    load(ss(i).name)
    i_tbl             = table();
    i_tbl.theta       = StimOri;
    i_tbl.delta       = sdp_acute_ang(nbk(i_tbl.theta)-i_tbl.theta);
    i_tbl.obs         = i.*ones(size(i_tbl.theta));
    i_tbl.resp        = RespOri;
    i_tbl.error       = sdp_acute_ang(i_tbl.resp-i_tbl.theta);
    i_tbl.trial       = (1:numel(i_tbl.obs))';
    i_tbl.cond        = double(Conf>median(Conf));
    tbl               = vertcat(tbl,i_tbl);
end

tbl_name        = 'Geurts et al. (2022)';
tbl.study       = repmat({tbl_name},numel(tbl.theta),1);
tbl.experiment  = repmat({'Experiment 1'},numel(tbl.theta),1);
tbl.stimulus    = repmat({'Orientation'},numel(tbl.theta),1);

fprintf('Storing the standard table \n\n');

cd('..\..\..\tables\single_experiments')
save Geurts_et_al_2022.mat tbl tbl_name
writetable(tbl,'Geurts_et_al_2022.csv','delimiter',';')

