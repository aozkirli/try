clear; clc; close all;

x = linspace(-150, 150, 1000); % orientation axis
sigma = 20; % base sigma for both current and context Gaussians
separations = 0:30:90; % degrees of separation between target and context
noise_level = 0.005; % adjust this for more or less noise

figure('Position', [100, 100, 1600, 600]);

for idx = 1:length(separations)
    sep = separations(idx);

    % Current (target) Gaussian at 0
    target_mu = 0;
    target = normpdf(x, target_mu, sigma);

    % Context Gaussian shifted by +sep degrees
    context_mu = sep;
    w =.5; % as sep increases, context degrades (adjust as needed)
    context_gauss = normpdf(x, context_mu, sigma);
    context_uniform = ones(size(x)) / range(x); % flat prior
    context = w * context_gauss + (1 - w) * context_uniform;
    context = context / trapz(x, context); % normalize

    %% ----- Bayesian Optimal Integration -----
    bayes = target .* (context);
    bayes = bayes / trapz(x, bayes); % normalize to make it a PDF

    % Calculate posterior mean and variance using trapz
    bayes_mu = trapz(x, x .* bayes);
    bayes_var = trapz(x, (x - bayes_mu).^2 .* bayes);
    bayes_sigma = sqrt(bayes_var);

    % Plot Bayesian (top row)
    subplot(2, length(separations), idx)
    hold on;
    fill([x fliplr(x)], [target zeros(size(target))], 'r', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    fill([x fliplr(x)], [context zeros(size(context))], 'b', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    plot(x, bayes, 'k', 'LineWidth', 2, 'Color', [0, 0, 0, 0.5]);
    plot([target_mu target_mu],[0 max(target)], 'k--', 'LineWidth', 2);
    plot([context_mu context_mu],[0 max(context)], 'k--', 'LineWidth', 2);
    yticks([])
    % Sigma text
    text(0.05, 0.85, ['\sigma_{posterior} = ', num2str(bayes_sigma, '%.2f')], 'Units', 'normalized', 'FontSize', 13);
    set(gca,'FontName','Arial','FontSize',15)

    if context_mu == target_mu
        xticks([context_mu])
        xticklabels('\mu_t = \mu_c')
    else
        xticks([target_mu context_mu])
        xticklabels({'\mu_t' '\mu_c'})
    end
    title(['Feature distance = ', num2str(sep), '°']);
    if idx == length(separations)
        legend('Target', 'Context', 'Bayesian');
    end
    ylim([0, 0.05]);

    %% ----- Interference Model (Summation) -----
    context = normpdf(x, context_mu, sigma);
    overlap = trapz(min([target' context']'))/5;
    % tmp =max(target,context);
    % overlap = min(tmp(target_mu+50:context_mu+50))./max(target)
    summed = target + context;
    
    % Add noise to interference summed curve
    summed_noisy = summed + noise_level * randn(size(summed));
    summed_noisy(summed_noisy < 0) = 0; % enforce non-negative PDF
    summed_noisy = summed_noisy./max(summed)*max(target);

    
    % Fit Gaussian to noisy interference (weighted)
    interference_fit = fit(x', summed_noisy', 'gauss1');
    interference_sigma = (1-overlap)*sigma + overlap*interference_fit.c1 / sqrt(2);

    % Plot Interference (bottom row)
    subplot(2, length(separations), idx + length(separations))
    hold on;
    fill([x fliplr(x)], [target zeros(size(target))], 'r', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    fill([x fliplr(x)], [context zeros(size(context))], 'b', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    plot(x, summed_noisy, 'k', 'LineWidth', 2, 'Color', [0, 0, 0, 0.2]);
    plot([target_mu target_mu],[0 max(target)], 'k--', 'LineWidth', 2);
    plot([context_mu context_mu],[0 max(context)], 'k--', 'LineWidth', 2);
    yticks([])
    % Sigma text
    text(0.05, 0.85, ['\sigma_{readout} = ', num2str(interference_sigma, '%.2f')], 'Units', 'normalized', 'FontSize', 13);

    if context_mu == target_mu
        xticks([context_mu])
        xticklabels('\mu_t = \mu_c')
    else
        xticks([target_mu context_mu])
        xticklabels({'\mu_t' '\mu_c'})
    end
    if idx == length(separations)
        legend('Target', 'Context', 'Summed (noisy)');
    end
    ylim([0, 0.05]);
    set(gca,'FontName','Arial','FontSize',15)
end

