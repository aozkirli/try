function [sd_data, sd_log] = load_data(cfg)

subjects = cfg.subjects;

sd_data = cell(size(subjects,2),1);
sd_log = cell(size(subjects,2),1);

counter = 1;

for iSubject = subjects
    
    
    subject_data = [];
    subject_log  = [];
    
    for iSession = 1:2
        
        for iBlock = 1:6
            
            tmp = load([cfg.dataPath '/Subj_' sprintf( '%03d', iSubject) '/data_S' num2str(iSubject) '_Session_' num2str(iSession) '_Block_' num2str(iBlock) '.mat']);
            subject_data = [subject_data; cell2mat(tmp.data)];
            
            tmp2 = load([cfg.dataPath '/Subj_' sprintf( '%03d', iSubject) '/Log/log_S' num2str(iSubject) '_Session_' num2str(iSession) '_Block_' num2str(iBlock) '.mat']);
            subject_log{iBlock,1} = tmp2.log;
            
        end
        
    end
    
    sd_data{counter} = subject_data;
    sd_log{counter}  = subject_log;
    
    counter = counter + 1;
    
    
end

end