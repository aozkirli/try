function moving_averages = compute_moving_averages_variance(cfg)

sd_data     = cfg.data;
bin_width   = cfg.bin_width;

nBack = size(sd_data.all,1);

for iBack = 1:nBack

    for iSubject = 1:size(sd_data.all{iBack},1)
        
        % pooled over locations
        x = sd_data.all{iBack}{iSubject}(:,3);
        y = sd_data.all{iBack}{iSubject}(:,4);
        
        x_padded = [x - 180; x; x + 180];
        y_padded = [y; y; y];
        
        moving_averages_all{iBack,1}{iSubject,1} = zeros(1,181);
        
        for b = 0:180
            moving_averages_all{iBack}{iSubject}(b+1) = ...
                circ_rad2ang(circ_std(circ_ang2rad(y_padded(inrange(x_padded, [1, bin_width] - floor(bin_width/2) - 1 - 90 + b)))));
        end
        
        clear x; clear y; clear x_padded; clear y_padded;
        
        if cfg.demean == true
            
            moving_averages_all{iBack}{iSubject} = ...
                moving_averages_all{iBack}{iSubject} - ...
                mean(moving_averages_all{iBack}{iSubject});
            
        end
        
        
        % same location
        x = sd_data.same_loc{iBack}{iSubject}(:,3);
        y = sd_data.same_loc{iBack}{iSubject}(:,4);
        
        x_padded = [x - 180; x; x + 180];
        y_padded = [y; y; y];
        
        moving_averages_same_loc{iBack,1}{iSubject,1} = zeros(1,181);
        
        for b = 0:180
            moving_averages_same_loc{iBack}{iSubject}(b+1) = ...
                circ_rad2ang(circ_std(circ_ang2rad(y_padded(inrange(x_padded, [1, bin_width] - floor(bin_width/2) - 1 - 90 + b)))));
        end
        
        clear x; clear y; clear x_padded; clear y_padded;
        
        if cfg.demean == true
            
            moving_averages_same_loc{iBack}{iSubject} = ...
                moving_averages_same_loc{iBack}{iSubject} - ...
                mean(moving_averages_same_loc{iBack}{iSubject});
            
        end
        
        % different location
        x = sd_data.diff_loc{iBack}{iSubject}(:,3);
        y = sd_data.diff_loc{iBack}{iSubject}(:,4);
        
        x_padded = [x - 180; x; x + 180];
        y_padded = [y; y; y];
        
        moving_averages_diff_loc{iBack,1}{iSubject,1} = zeros(1,181);
        
        for b = 0:180
            moving_averages_diff_loc{iBack}{iSubject}(b+1) = ...
                circ_rad2ang(circ_std(circ_ang2rad(y_padded(inrange(x_padded, [1, bin_width] - floor(bin_width/2) - 1 - 90 + b)))));
        end
        
        clear x; clear y; clear x_padded; clear y_padded;
        
        if cfg.demean == true
            
            moving_averages_diff_loc{iBack}{iSubject} = ...
                moving_averages_diff_loc{iBack}{iSubject} - ...
                mean(moving_averages_diff_loc{iBack}{iSubject});
            
        end
        
    end

% Grand moving averages


grand_moving_averages.all{iBack,1} = ...
    circ_rad2ang(circ_mean(circ_ang2rad(cell2mat(moving_averages_all{iBack})),[],1));

grand_moving_averages.same_loc{iBack,1} = ...
    circ_rad2ang(circ_mean(circ_ang2rad(cell2mat(moving_averages_same_loc{iBack})),[],1));

grand_moving_averages.diff_loc{iBack,1} = ...
    circ_rad2ang(circ_mean(circ_ang2rad(cell2mat(moving_averages_diff_loc{iBack})),[],1));

end

moving_averages.all         = moving_averages_all;
moving_averages.same_loc    = moving_averages_same_loc;
moving_averages.diff_loc    = moving_averages_diff_loc;

moving_averages.grand_all       = grand_moving_averages.all;
moving_averages.grand_same_loc  = grand_moving_averages.same_loc;
moving_averages.grand_diff_loc  = grand_moving_averages.diff_loc;




end

