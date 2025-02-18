function [sd_data, sd_log] = load_data(cfg)

subjects = cfg.subjects;

sd_data = cell(size(subjects,2),1);
sd_log = cell(size(subjects,2),1);

for i = subjects
    
    
    
    
    tmp = load([cfg.dataPath 'Data - Subj' num2str(i) filesep 'totaldata_Subj' num2str(i) '.mat']);
    
    sd_subject_data = cell(length(tmp.totaldata.dataPath), 1);
    sd_subject_log = cell(length(tmp.totaldata.dataPath), 1);
    
    for b = 1:length(sd_subject_data)
        
        tmppath1 = tmp.totaldata.dataPath{b};
        tmppath2 = tmp.totaldata.logPath{b};
        if ~any(strfind(pwd,'/'))
            tmppath1(strfind(tmppath1,'/')) = '\';
            tmppath2(strfind(tmppath2,'/')) = '\';
        end
        tmp1 = load([cfg.dataPath, tmppath1]);
        tmp2 = load([cfg.dataPath, tmppath2]);
%         tmp1 = load(['..\Data\', tmppath1]);
%         tmp2 = load(['..\Data\', tmppath2]);
       
        sd_subject_data{b} = cell2mat(tmp1.data);
        sd_subject_log{b} = tmp2.log;
        
    end
    
    sd_subject_data = cell2mat(sd_subject_data);
    
    % necessary to make stimulus angle measure consistent with response
    % angle measure (different definitions of zero angle)
    sd_subject_data(:,4) = mod(-sd_subject_data(:, 4) + 180,180) - 90;
    
    
    
    
    sd_data{i,1} = sd_subject_data;
    sd_log{i,1} = sd_subject_log;
    
end

% delete empty cells in cell arrays
sd_data = sd_data(~cellfun('isempty',sd_data));
sd_log = sd_log(~cellfun('isempty',sd_log));

  
end