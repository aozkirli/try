function moving_averages = compute_moving_averages(cfg)

sd_data     = cfg.data;
bin_width   = cfg.bin_width;

nBack = size(sd_data.all,1);

for iBack = 1:nBack

    for iSubject = 1:size(sd_data.all{iBack},1)
          
        x = sd_data.all{iBack}{iSubject}(:,3);
        y = sd_data.all{iBack}{iSubject}(:,4);
        
        x_padded = [x - 180; x; x + 180];
        y_padded = [y; y; y];
        
        moving_averages_all{iBack,1}{iSubject,1} = zeros(1,181);
        
        for b = 0:180
            moving_averages_all{iBack}{iSubject}(b+1) = ...
                circ_rad2ang(circ_mean(circ_ang2rad(y_padded(inrange(x_padded, [1, bin_width] - floor(bin_width/2) - 1 - 90 + b)))));
        end
        
        clear x; clear y; clear x_padded; clear y_padded;
              
        
    end

% Grand moving averages


grand_moving_averages.all{iBack,1} = ...
    circ_rad2ang(circ_mean(circ_ang2rad(cell2mat(moving_averages_all{iBack})),[],1));


end

moving_averages.all         = moving_averages_all;

moving_averages.grand_all       = grand_moving_averages.all;


end

