function sorted_data = sort_data(cfg)

data = cfg.data;

% allocating space for each subject
data_nback           = cell(cfg.nBack,1);

for iSubject = 1:size(data,1)
    
    
    %% nback sorting
    for iBack = 1:cfg.nBack
        
        
        tmp_data = [];
        
        for j = (iBack+1):size(data{iSubject,1},1)
            
            if (data{iSubject,1}(j,2) == data{iSubject,1}(j-iBack,2)) && (~isnan(data{iSubject,1}(j,5)) && ~isnan(data{iSubject,1}(j-iBack,5)))
              
                nback_ori   = data{iSubject,1}(j-iBack,4);
                current_ori = data{iSubject,1}(j,4);
                ori_diff    = mod(nback_ori - current_ori + 90, 180) - 90;
                error       = mod(data{iSubject,1}(j,5) - data{iSubject,1}(j,4) + 90, 180) - 90;
                
                tmp_data = [tmp_data; nback_ori current_ori ori_diff error];
                
            end
            
        end
        
        if cfg.demean
            tmp_data(:,3) = bsxfun(@minus,tmp_data(:,3),circ_rad2ang(circ_mean(circ_ang2rad(tmp_data(:,3)))));
        end
        
        data_nback{iBack}{iSubject,1} = tmp_data;
        
        
    end
    
    
    
end

sorted_data.all = data_nback;

end