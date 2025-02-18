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


subjects = [2,3,4,5,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,25,26,...
    27,28,29,30,31,32,33,35,36,37,38];

cfg = [];
cfg.subjects = subjects;
cfg.cd  = pwd;
[data_ori,~,data_size,~] = load_data(cfg);

fprintf('Converting to standard table format...\n\n');

tbl         = table;
for k = 1:numel(data_ori)
    tmp     = data_ori{k};
    tbl_i   = table;
    tbl_i.obs = tmp(:,1);
    tbl_i.block = tmp(:,2);
    tbl_i.trial = tmp(:,4);
    tbl_i.theta = round(tmp(:,10)+90); % test grating
    tbl_i.delta = sdp_acute_ang(round(tmp(:,5)+90)-tbl_i.theta); % delta wrt to inducer (see readme)
    tbl_i.rt    = tmp(:,12);
    tbl_i.resp  = tmp(:,11)+90;
    tbl_i.error = sdp_acute_ang(tbl_i.resp-tbl_i.theta);
    tbl_i.cond  = double(tmp(:,3)==1); % orientation discrimination task (see readme)
    tbl_ori     = tbl_i;
    
    % size task
    tmp     = data_size{k};
    tbl_i   = table;
    tbl_i.obs = tmp(:,1);
    tbl_i.block = tmp(:,2);
    tbl_i.trial = tmp(:,4);
    tbl_i.theta = round(tmp(:,10)+90); % test grating
    tbl_i.delta = sdp_acute_ang(round(tmp(:,5)+90)-tbl_i.theta); % delta wrt to inducer (see readme)
    tbl_i.rt    = tmp(:,12);
    tbl_i.resp  = tmp(:,11)+90;
    tbl_i.error = sdp_acute_ang(tbl_i.resp-tbl_i.theta);
    tbl_i.cond  = double(tmp(:,3)==1); % orientation discrimination task (see readme)
    tbl_size    = tbl_i;
    
    tbl_both    = vertcat(tbl_ori,tbl_size);
    % combine
    tbl         = vertcat(tbl,tbl_both);
end

tbl_name            = 'Fritsche & deLange (2019)';
tbl.study           = repmat({tbl_name},numel(tbl.theta),1);
tbl.experiment      = repmat({'Experiment 1'},numel(tbl.theta),1);
tbl.stimulus        = repmat({'Orientation'},numel(tbl.theta),1);

fprintf('Storing the standard table \n\n');

cd('..\..\tables\single_experiments')
fname         = ['Fritsche_deLange_2019'];
save(fname,'tbl','tbl_name')
writetable(tbl,[fname '.csv'],'delimiter',';')



function [data_ori, log_ori, data_size, log_size] = load_data(cfg)

subjects = cfg.subjects;

data_ori        = cell(size(subjects,2),1);
log_ori         = cell(size(subjects,2),1);
data_size       = cell(size(subjects,2),1);
log_size        = cell(size(subjects,2),1);


counter = 1;

for iSubject = subjects
    
   
    disp(['Loading data of S' num2str(iSubject) ' ...'])
    
    
    subject_data_ori = [];
    subject_log_ori  = [];
    subject_data_size = [];
    subject_log_size  = [];
    
    
    for iBlock = 1:8
        
        
        tmp = load([cfg.cd '\Data\S' num2str(iSubject) '\data_S' num2str(iSubject) '_Task_1_Block_' num2str(iBlock) '.mat']);
        subject_data_ori = [subject_data_ori; cell2mat(tmp.data)];
        
        tmp2 = load([ cfg.cd '\Data\S' num2str(iSubject) '\Log\log_S' num2str(iSubject) '_Task_1_Block_' num2str(iBlock) '.mat']);
        subject_log_ori{iBlock,1} = tmp2.log;
       
        clear tmp tmp2
        
       
        tmp = load([ cfg.cd '\Data\S' num2str(iSubject) '\data_S' num2str(iSubject) '_Task_2_Block_' num2str(iBlock) '.mat']);
        subject_data_size = [subject_data_size; cell2mat(tmp.data)];
        
        tmp2 = load([ cfg.cd '\Data\S' num2str(iSubject) '\Log\log_S' num2str(iSubject) '_Task_2_Block_' num2str(iBlock) '.mat']);
        subject_log_size{iBlock,1} = tmp2.log;
        
        clear tmp tmp2
        
        
    end
    
    
    
    data_ori{counter} = subject_data_ori;
    log_ori{counter}  = subject_log_ori;
    
    
    data_size{counter} = subject_data_size;
    log_size{counter}  = subject_log_size;
    
    
    counter = counter + 1;
    
end


end




