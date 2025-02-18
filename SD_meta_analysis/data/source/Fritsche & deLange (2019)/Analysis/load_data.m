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
        
        
        tmp = load(['../Data/S' num2str(iSubject) '/data_S' num2str(iSubject) '_Task_1_Block_' num2str(iBlock) '.mat']);
        subject_data_ori = [subject_data_ori; cell2mat(tmp.data)];
        
        tmp2 = load(['../Data/S' num2str(iSubject) '/Log/log_S' num2str(iSubject) '_Task_1_Block_' num2str(iBlock) '.mat']);
        subject_log_ori{iBlock,1} = tmp2.log;
       
        clear tmp tmp2
        
       
        tmp = load(['../Data/S' num2str(iSubject) '/data_S' num2str(iSubject) '_Task_2_Block_' num2str(iBlock) '.mat']);
        subject_data_size = [subject_data_size; cell2mat(tmp.data)];
        
        tmp2 = load(['../Data/S' num2str(iSubject) '/Log/log_S' num2str(iSubject) '_Task_2_Block_' num2str(iBlock) '.mat']);
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