function [cohen_d, p, BF] = weighted_paired_ttest(x, y, weights)
    % Inputs:
    % x - First sample (vector)
    % y - Second sample (vector)
    % weights - Weights for each paired observation
    
    % Check input dimensions
    if numel(x) ~= numel(y) || numel(x) ~= numel(weights)
        error('x, y, and weights must have the same length.');
    end
    
    % % Compute differences
    % d = x - y;
    % 
    % % Weighted mean of differences
    % w_mean = nansum(weights .* d) / nansum(weights);
    % 
    % % Weighted variance of differences
    % w_var = nansum(weights .* (d - w_mean).^2) / nansum(weights);
    % 
    % % Weighted standard deviation
    % w_sd = sqrt(w_var);
    % 
    % % Weighted t-statistic
    % t = w_mean / sqrt(w_var / length(d));
    % 
    % % Degrees of freedom (approximation)
    % df = (nansum(weights)^2) / nansum(weights.^2);
    
    % % Compute p-value using t-distribution
    % p = 2 * (1 - tcdf(abs(t), df));
    % 

    [~,p,~,stats]= ttest(x,y);
    
    % Cohen's d
    cohen_d = round(nanmean(x-y) / nanstd(x-y) ,2);

    % Bayes Factor (BF) calculation (from See Rouder et al, 2009)
    t = stats.tstat;
    n = length(x);
    r = 0.707;
    % Function to be integrated
    F = @(g,t,n,r) (1+n.*g.*r.^2).^(-1./2) .* (1 + t.^2./((1+n.*g.*r.^2).*(n-1))).^(-n./2) .* (2.*pi).^(-1./2) .* g.^(-3./2) .* exp(-1./(2.*g));
    % Bayes factor calculation
    BF = (1 + t^2/(n-1))^(-n/2) / integral(@(g) F(g,t,n,r),0,Inf);
    % Invert Bayes Factor
    BF = 1 / BF;

    if BF>100
        BF = "_{10} > 100";
    elseif BF < 1/100
        BF = "_{01} > 100";
    elseif BF >= 1
        BF = ['_{10} = ' num2str(round(BF,2))]; 
    else
        BF = ['_{01} = ' num2str(round(1/BF,2))]; 
    end
end
