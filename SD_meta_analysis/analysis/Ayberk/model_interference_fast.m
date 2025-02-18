clear all;close all;clc;
figure('Units', 'normalized', 'Position', [0 0 0.8 0.3]);

%% constants
curr = 0;
delta = -90:90;
sds         = 5:1:20;
colors = blue2pink(length(sds));

for ss = 1:length(sds)
    %% model parameters
    sigma = sds(ss);
    wp = .2;     % previous weight

    %% functions
    sdk             = @(x) sd2k(deg2rad(x*2));     % std° to Von Mises K
    nrmsum          = @(x) x./sum(x);
    vonmises        = @(x,s) nrmsum(exp(sdk(s).*(cosd(2*x-(2*delta))-1)));
    create_one_dist = @(delta,par) vonmises(par(1),par(2));
    create_two_dist = @(delta,par) vonmises(par(1),par(3)).*wp + vonmises(par(2),par(3)).*(1-wp);
    options         = optimoptions('fmincon','Display','off');

    cc = create_one_dist(delta,[curr sigma]);
    for n = 1
        for prev = -90:90
            dist = create_two_dist(delta,[prev curr sigma]);
            mu = trapz(delta, delta .* dist);
            second_moment = trapz(delta, delta.^2 .* dist); % E(x^2)
            variance = second_moment - mu^2;
            s_fat    = sqrt(variance);
            pp = create_one_dist(delta,[prev sigma]);
            overlap_scaled(prev+91) = trapz(min([pp' cc']'));
            % percent overlap times one distribution  + 1-percent overlap times
            % actual values
            bias(prev+91)           = overlap_scaled(prev+91)*mu     +  (1-overlap_scaled(prev+91))*curr;
            std_dev(prev+91)        = overlap_scaled(prev+91)*s_fat  +  (1-overlap_scaled(prev+91))*sigma;
        end
    end
    subplot(131)
    plot(delta,overlap_scaled,'Color',colors(ss,:))
    xlim([0 90])
    xticks(-90:45:90)
    hold on;
    subplot(132)
    plot(delta,bias,'Color',colors(ss,:))
    xticks(-90:45:90)
    xlim([0 90])
    hold on;
    subplot(133)
    plot(delta,std_dev./sigma,'Color',colors(ss,:))
    xticks(-90:45:90)
    xlim([0 90])
    hold on;
end
