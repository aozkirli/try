function sorted_data = sort_data(cfg)

data = cfg.data;

% allocating space for each subject
data_nback          = cell(cfg.nBack,1);
data_nback_short    = cell(cfg.nBack,1);
data_nback_long     = cell(cfg.nBack,1);

for iSubject = 1:size(data,1)
    
    valid_trials = not(isnan(data{iSubject,1}(:,6)));
    
    response_errors = mod(data{iSubject,1}(:,5) - data{iSubject,1}(:,4) + 90, 180) - 90;
    
    if cfg.demean
        response_errors = bsxfun(@minus,response_errors,...
            circ_rad2ang(circ_mean(circ_ang2rad(response_errors(valid_trials)))));
    end
    
    
    %% nback sorting
    for iBack = 1:cfg.nBack
        
        
        %% all trials
        tmp_data = [];
        
        for j = (iBack+1):size(data{iSubject,1},1)
            
            if (data{iSubject,1}(j,2) == data{iSubject,1}(j-iBack,2)) && (~isnan(data{iSubject,1}(j,6)) && ~isnan(data{iSubject,1}(j-iBack,6)))
              
                nback_ori   = data{iSubject,1}(j-iBack,4);
                current_ori = data{iSubject,1}(j,4);
                ori_diff    = mod(nback_ori - current_ori + 90, 180) - 90;
                error       = response_errors(j);
                
                tmp_data = [tmp_data; nback_ori current_ori ori_diff error];
                
            end
            
        end
        
        data_nback{iBack}{iSubject,1} = tmp_data;
        
        
        %% current short trials
        tmp_data = [];
        
        for j = (iBack+1):size(data{iSubject,1},1)
            
            if (data{iSubject,1}(j,2) == data{iSubject,1}(j-iBack,2)) && ...
                    (~isnan(data{iSubject,1}(j,6)) && ~isnan(data{iSubject,1}(j-iBack,6))) && ...
                    (data{iSubject,1}(j,7) == 1)
              
                nback_ori   = data{iSubject,1}(j-iBack,4);
                current_ori = data{iSubject,1}(j,4);
                ori_diff    = mod(nback_ori - current_ori + 90, 180) - 90;
                error       = response_errors(j);
                
                tmp_data = [tmp_data; nback_ori current_ori ori_diff error];
                
            end
            
        end
        
        data_nback_short{iBack}{iSubject,1} = tmp_data;
        
        
        %% current long trials         
        tmp_data = [];
        
        for j = (iBack+1):size(data{iSubject,1},1)
            
            if (data{iSubject,1}(j,2) == data{iSubject,1}(j-iBack,2)) && ...
                    (~isnan(data{iSubject,1}(j,6)) && ~isnan(data{iSubject,1}(j-iBack,6))) && ...
                    (data{iSubject,1}(j,7) == 2)
              
                nback_ori   = data{iSubject,1}(j-iBack,4);
                current_ori = data{iSubject,1}(j,4);
                ori_diff    = mod(nback_ori - current_ori + 90, 180) - 90;
                error       = response_errors(j);
                
                tmp_data = [tmp_data; nback_ori current_ori ori_diff error];
                
            end
            
        end
        
        data_nback_long{iBack}{iSubject,1} = tmp_data;
        
        
        
    end  
    
    sorted_data.data_nback          = data_nback;
    sorted_data.data_nback_short    = data_nback_short;
    sorted_data.data_nback_long     = data_nback_long;
    

end