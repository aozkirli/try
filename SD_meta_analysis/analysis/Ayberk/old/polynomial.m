tbl_scatter = [];
all_curve= [];
pol_deg = [];
aggregated_delta = [];
aggregated_scatt = [];
aggregated_counts = [];
obsnum = [];
oo =1;
error_variable = 'error_ori_deb_sd_deb';
close all;
% figure;
for i = 1:length(unique(tbl.studynum))
    tbl_i             = tbl_subset(tbl,'studynum',i);
    tbl_i             = tbl_i(abs(tbl_i.delta)<=90,:);
    % fix decimals here
    tbl_i.theta       = round(tbl_i.theta);
    tbl_i.delta       = round(tbl_i.delta);
    if any(strcmp(tbl_i.study,'Fischer et al. (2020)'))%length(unique(abs(tbl_i.delta)))<=10
        continue
    end
    % loop over exps and conditions
    nexps           = max(tbl_i.expnum);
    for k = 1:nexps
        tbl_i_k       = tbl_subset(tbl_i,'expnum',k);
        cond          = unique(tbl_i_k.cond);
        ncond         = numel(cond);
        for j = 1:ncond
            tbl_i_k_j     = tbl_subset(tbl_i_k,'cond',cond(j));
            obs           = unique(tbl_i_k_j.obs);
            nobs          = numel(obs);
            for o = 1:nobs
                tbl_i_k_j_o              = tbl_subset(tbl_i_k_j,'obs',obs(o));
                % Polynomial fittingb
                scatt= grpstats(tbl_i_k_j_o.(error_variable),abs(tbl_i_k_j_o.delta),'std');
                [counts,grpp] = groupcounts(abs(tbl_i_k_j_o.delta));
                delta_fit = nanunique(abs(tbl_i_k_j_o.delta));
                scatt_fit = scatt(scatt~=0 & ~isnan(scatt));
                delta_fit = delta_fit(scatt~=0 & ~isnan(scatt));
                counts = counts(scatt~=0 & ~isnan(scatt));
                weights = counts./max(counts);

                % Min and maximum polynomial degree to test
                if length(counts)<=10
                    degrees = [2:min(length(counts)-1,3)];
                else
                    degrees = [2:min(length(counts)-1,9)];
                end

                % Normalize delta_fit
                delta_fit_mean = mean(delta_fit);
                delta_fit_std = std(delta_fit);
                delta_fit_normalized = (delta_fit - delta_fit_mean) / delta_fit_std;

                % Initialize storage for BIC values
                BIC_values = [];

                % Loop over polynomial degrees
                for degree = degrees
                    % Define a weighted polynomial fitting model
                    fitOptions = fitoptions('Method', 'LinearLeastSquares', ...
                        'Weights', weights);  % Apply weights

                    % Perform the weighted fit on normalized data
                    model = fit(delta_fit_normalized, scatt_fit, ['poly' num2str(degree)], fitOptions);

                    % Compute fitted values on normalized data
                    y_fit = feval(model, delta_fit_normalized);

                    % Compute weighted residuals
                    residuals = scatt_fit - y_fit;
                    weighted_residuals = sqrt(counts) .* residuals;

                    % Compute weighted sum of squared residuals
                    SSR = sum(weighted_residuals.^2);

                    % Calculate number of parameters in the model
                    numParams = degree + 1;  % Polynomial degree + 1 for intercept

                    % Compute BIC
                    n = length(delta_fit);  % Number of data points
                    BIC = n * log(SSR / n) + numParams * log(n);
                    BIC_values = [BIC_values;BIC];
                end

                % Find the best model (lowest BIC value)
                [sortedBIC, sortedIndex] = sort(BIC_values);
                bestDegree = degrees(min(sortedIndex(sortedBIC<sortedBIC(1)+2)));

                % Display results
                % fprintf('Best-fitting polynomial degree: %d\n', bestDegree);

                % Perform the weighted fit on normalized data with the best degree
                best_model = fit(delta_fit_normalized, scatt_fit, ['poly' num2str(bestDegree)], fitOptions);

                % Compute fitted values for the normalized range
                best_fit = feval(best_model, ([0:90]'-delta_fit_mean)./delta_fit_std);

                % if bestDegree==6
                %     disp('here')
                % elseif bestDegree>3
                %     disp('here')
                    % Plot results
                    scatter(delta_fit, scatt_fit, 10*counts, 'filled');  % Weighted scatter plot
                    hold on;
                    plot(0:90, best_fit, 'r-', 'LineWidth', 1.5);
                    xlabel('delta\_fit');
                    ylabel('scatt\_fit');
                    title(sprintf('Best Polynomial Fit (Degree %d)', bestDegree));
                    legend('Data', 'Best Fit', 'Location', 'Best');
                    grid on;
                    hold off;
                % end

                if any(best_fit([1 46 91])<=0)
                    disp('negative scatter!')
                    continue
                end
                pol_deg = [pol_deg;bestDegree];
                all_curve = [all_curve best_fit];
                % Collect scatter and delta values
                aggregated_delta = [aggregated_delta; delta_fit];
                aggregated_scatt = [aggregated_scatt; scatt_fit];
                aggregated_counts = [aggregated_counts; counts];
                obsnum            = [obsnum; oo*ones(size(scatt_fit))];
                oo                = oo+1;
                % store
                tmp                      = repmat(tbl_i_k_j_o(1,ismember(tbl_i_k_j_o.Properties.VariableNames,{'obsid','codenum'})),3,1);
                tmp.bin                  = categorical(1:3)';
                tt     = grpstats(tbl_i_k_j_o(:,{'codenum','stimtype','obsid','bin' error_variable}),{'codenum','stimtype','obsid','bin'},'std');
                tmp.poly_scatter          = best_fit([1 46 91]);
                tmp.bin_scatter   = tt.(['std_' error_variable]);
                tbl_scatter = [tbl_scatter;tmp];
            end
        end
    end
end
%%
figure
titles = {'íso' 'mid' 'ortho'}
colors          = linspecer(total_datasets);%
for i = 1:3
    subplot(1,3,i)
    poly = tbl_scatter.poly_scatter(i:3:end);
    binned = tbl_scatter.bin_scatter(i:3:end);
    idx = ~isnan(poly) & ~isnan(binned) 
    [rr,pp] = corr(binned(idx),poly(idx))
    plot(binned(idx),poly(idx),'ko');
    hold on
    l1 = lsline;
    l1.Color = 'k';
    l1.LineWidth = 2;
    % for c = i:3:height(tbl_scatter)
    %     sc = scatter(tbl_scatter.bin_scatter(c),tbl_scatter.poly_scatter(c),'o','Color',colors(tbl_scatter.codenum,:));
    % end
    xlabel('binned')
    ylabel('polynomial fit')
    % xlim([0 80]);ylim([0 80]);
    
    title([titles{i} ': r = ' num2str(round(rr,2)) ', p = ' num2str(pp)])
    set(gca,'fontsize',15)
end
% figure
% ksdensity(tbl_scatter.poly_scatter(3:3:end)-tbl_scatter.poly_scatter(1:3:end))
% hold on;
% ksdensity(tbl_scatter.bin_scatter(3:3:end)-tbl_scatter.bin_scatter(1:3:end))

%% 
figure
x = [0:90]; % X-axis
out = (isoutlier(tbl_scatter.poly_scatter(1:3:end),'median') | isoutlier(tbl_scatter.poly_scatter(2:3:end),'median') | isoutlier(tbl_scatter.poly_scatter(3:3:end),'median'))%zeros(size(all_curve,2),1);%isoutlier(tt.std_error_ori_deb_sd_deb,'median');%isoutlier(mean([tbl_scatter.poly_scatter(1:3:end) tbl_scatter.poly_scatter(2:3:end) tbl_scatter.poly_scatter(3:3:end)],2),'quartiles');
(isoutlier(tbl_scatter.poly_scatter(1:3:end),'median') | isoutlier(tbl_scatter.poly_scatter(2:3:end),'median') | isoutlier(tbl_scatter.poly_scatter(3:3:end),'median'));
weighted_agg_scatt=aggregated_scatt(ismember(obsnum,find(~out)));
weighted_agg_delta=aggregated_delta(ismember(obsnum,find(~out)));
weighted_agg_counts=aggregated_counts(ismember(obsnum,find(~out)));

mean_values = nanmean(all_curve(:,~out),2)'; % Mean values

alpha = 0.05; % Significance level for 95% confidence
% Bootstrap confidence intervals
num_bootstraps = 1000; % Number of bootstrap samples
bootstrap_samples = bootstrp(num_bootstraps, @mean, all_curve(:,~out)'); % Bootstrap resamples

% % Calculate percentiles for the confidence interval
ci_lower = prctile(bootstrap_samples, 100 * alpha / 2);
ci_upper = prctile(bootstrap_samples, 100 * (1 - alpha / 2));

% Calculate the standard deviation for each delta on the aggregate dataset by weighing:
all_scatter = grpstats(weighted_agg_scatt.*weighted_agg_counts,weighted_agg_delta,'sum')./grpstats(weighted_agg_counts,weighted_agg_delta,'sum');
scatter(x,all_scatter(1:91),'filled','MarkerFaceColor','#FFA500','MarkerFaceAlpha',.5)
hold on;
% Plot shaded area and mean line
fill([x, fliplr(x)], [ci_upper, fliplr(ci_lower)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(x, mean_values, 'k', 'LineWidth', 2);
xlabel('Stimulus dissimilarity index'); ylabel('Error Scatter'); title(['Best polynomial fit (' num2str(degrees) ')']);
legend('95% CI', 'mean', 'Location', 'best');
grid on;
set(gca,'FontSize',15)
xlim([-1 91])
xticks([0 45 90])

% 
% fitOptions = fitoptions('Method', 'LinearLeastSquares', ...
%     'Weights', grpstats(aggregated_counts,aggregated_delta,'mean'));  % Apply weights
% % Perform the weighted fit on normalized data
% model = fit(x', all_scatter, 'poly4', fitOptions);
% % Compute fitted values on normalized data
% y_fit = feval(model, x);
% plot(x,y_fit)


lmm             = fitlme(tbl_scatter,'poly_scatter ~ bin + (1+bin|obsid) + (1+bin|codenum) + (1+bin|codenum:obsid)');
disp(lmm)
tbl_scatter.new_bin = repmat(categorical([1 2 1]'),height(tbl_scatter)/3,1);
reduced_lmm = fitlme(tbl_scatter, 'poly_scatter ~ new_bin + (1+new_bin|obsid) + (1+new_bin|codenum) + (1+new_bin|codenum:obsid)');

bic_full = lmm.ModelCriterion.BIC;
bic_reduced = reduced_lmm.ModelCriterion.BIC;
dBIC = bic_reduced-bic_full
BF_reduced_over_full = 1/exp((bic_reduced-bic_full)/2)

lmm             = fitlme(tbl_scatter,'bin_scatter ~ bin + (1+bin|obsid) + (1+bin|codenum) + (1+bin|codenum:obsid)');
disp(lmm)
tbl_scatter.new_bin = repmat(categorical([1 2 1]'),height(tbl_scatter)/3,1);
reduced_lmm = fitlme(tbl_scatter, 'bin_scatter ~ new_bin + (1+new_bin|obsid) + (1+new_bin|codenum) + (1+new_bin|codenum:obsid)');

bic_full = lmm.ModelCriterion.BIC;
bic_reduced = reduced_lmm.ModelCriterion.BIC;
dBIC = bic_reduced-bic_full
BF_reduced_over_full = 1/exp((bic_reduced-bic_full)/2)


