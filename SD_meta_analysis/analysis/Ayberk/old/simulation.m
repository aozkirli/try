clc; clear; close all;
load('/Users/ayberkozkirli/Library/CloudStorage/GoogleDrive-ayberk.ozkirli@epfl.ch/.shortcut-targets-by-id/14P9enoQnGTKLDkkJnld8K5Wh-X7GmXcB/EPFL_Ambizione/Experiments/SD_meta_analysis/data/datasets/SD_ma_master_table.mat')
% Define true polynomial coefficients
coefficients = [-0.000679, 0.06389, 9.25]; % Polynomial coefficients
polynomial = @(x) coefficients(1)*x.^2 + coefficients(2)*x + coefficients(3);
tmp = tbl(abs(tbl.delta)<=90,:);

% Define delta values
delta_values = (0:90)';

% Number of participants
num_participants = 1000;

% Initialize arrays for results
super_true_participant = zeros(num_participants, 1);
super_bin_participant = zeros(num_participants, 1);
super_poly_participant = zeros(num_participants, 1);

% Simulate participants
parfor p = 1:num_participants
    % Initialize participant data
    participant_scatter = [];
    delta_list = [];
    weights = [];
    
    for delta = delta_values'
        % Define mean and sigma for the distribution at this delta
        true_mean = polynomial(delta);

        sigma = 1.5 * rand; % Random sigma
        
        % Generate 100 values for this delta
        distribution = true_mean + sigma * randn(10000, 1);
        
        % Randomly sample between 10 and 50 values for this participant
        num_samples = randi([1, 20]);
        sampled_values = datasample(distribution, num_samples, 'Replace', false);
        
        % Store the sampled data and corresponding weights
        participant_scatter = [participant_scatter; sampled_values];
        delta_list = [delta_list; repmat(delta, num_samples, 1)];
        weights = [weights; ones(num_samples, 1) * num_samples];
    end
    plot(delta_list,participant_scatter,'k.');hold on;
    plot(0:90,polynomial(0:90),'b','LineWidth',2);
    
    % True superiority for this participant
    super_true_participant(p) = polynomial(0) - polynomial(90);

    % Bin-wise superiority for this participant
    iso_mean = mean(participant_scatter(delta_list <= 10));
    ortho_mean = mean(participant_scatter(delta_list >= 80));
    super_bin_participant(p) = iso_mean - ortho_mean;
    
    % Polynomial fitting for this participant
    degrees = 2:9;
    participant_scatter = grpstats(participant_scatter,delta_list);
    weights = grpstats(weights,delta_list);
    delta_list = unique(delta_list);
    
    delta_mean = mean(delta_list);
    delta_std = std(delta_list);
    delta_norm = (delta_list - delta_mean) / delta_std;
    
    BIC_values = zeros(size(degrees));
    models = cell(size(degrees)); % Store models for each degree
    
    for d = 1:length(degrees)
        degree = degrees(d);
        fitOptions = fitoptions('Method', 'LinearLeastSquares', 'Weights', weights);
        model = fit(delta_norm, participant_scatter, sprintf('poly%d', degree), fitOptions);
        models{d} = model;
        
        % Calculate residuals and BIC
        y_fit = feval(model, delta_norm);
        residuals = participant_scatter - y_fit;
        SSR = sum((sqrt(weights) .* residuals).^2);
        numParams = degree + 1;
        n = length(delta_list);
        BIC_values(d) = n * log(SSR / n) + numParams * log(n);
    end
    
    % Select the best polynomial degree
    [sortedBIC, sortedIndex] = sort(BIC_values);
    bestDegree = degrees(min(sortedIndex(sortedBIC < sortedBIC(1) + 2)))
    best_model = models{bestDegree-1};
    
    % Compute polynomial superiority for this participant
    delta_test = (0 - delta_mean) / delta_std; % Normalized 0
    delta_90 = (90 - delta_mean) / delta_std; % Normalized 90
    super_poly_participant(p) = feval(best_model, delta_test) - feval(best_model, delta_90);
    plot(0:90,feval(best_model, delta_norm),'r','LineWidth',2);
    hold off;
end

% Combine all participant-level differences
bin_differences = sqrt((super_bin_participant - super_true_participant).^2);
poly_differences = sqrt((super_poly_participant - super_true_participant).^2);

% Calculate mean and standard error for bin and polynomial differences
bin_mean = mean(bin_differences);
poly_mean = mean(poly_differences);
bin_std = 1.96 * std(bin_differences) / (sqrt(num_participants));
poly_std = 1.96 * std(poly_differences) / (sqrt(num_participants));

% Visualization
figure;
errorbar(1:2, [bin_mean, poly_mean], [bin_std, poly_std], 'k', 'LineWidth', 1);
hold on;
plot([0 3], [0 0], 'k--');
xlim([0 3]);
xticks(1:2);
xticklabels({'bin-true', 'poly-true'});
ylabel('Percent Deviation from True Superiority');
xlabel('Method');
title('Comparison of Bin and Polynomial Methods (Participant Level)');
grid on;

%% 
% Combine all participant-level differences
bin_differences = super_bin_participant - super_true_participant;
poly_differences = super_poly_participant - super_true_participant;

% Calculate mean and standard error for true, bin, and polynomial values
true_mean =  mean(super_true_participant) ;
bin_mean =  mean(super_bin_participant) ;
poly_mean = mean(super_poly_participant) ;

true_std =  std(super_true_participant) / sqrt(num_participants);
bin_std =  std(super_bin_participant)  / sqrt(num_participants);
poly_std =  std(super_poly_participant)  / sqrt(num_participants);

% Visualization
figure;
errorbar(1:3, [true_mean, bin_mean, poly_mean], [true_std, bin_std, poly_std], 'k', 'LineWidth', 1);
hold on;
plot([0 4], [0 0], 'k--'); % Plot reference line at 100% (true value baseline)
xlim([0.5 3.5]);
xticks(1:3);
xticklabels({'True', 'Bin', 'Poly'});
ylabel('Percent Deviation from True Superiority');
xlabel('Method');
title('Comparison of True, Bin, and Polynomial Methods');
grid on;

