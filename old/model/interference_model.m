function [model,sdp] = interference_model()
theta= [0 20 90];
Se = 10;
Sd = 10;
Pt = .5;
tau = 5;
nback = 1;
bins=21;
meas = 0;

trials = numel(theta);

%% Function Handles
dSd2k = @(x) sd2k(deg2rad(x) * 2); % Convert std° to concentration 'k'
sumnorm = @(x) x ./ sum(x); % Normalize to sum 1
ch = 0:179; % Orientation channels
chresp = @(ori, Se, ch) ...
    sumnorm(exp(dSd2k(Se) .* (cosd(2 * ch - (2 * ori)) - 1))); % Von Mises
uniform = ones(size(ch)) ./ numel(ch); % Uniform distribution
decoder = @(input, ch) mod(rad2deg(angle(sum(exp(1j * deg2rad(2 * ch)) .* input))) / 2, 180); % Decode posterior
decay = @(x, tau) (exp(x .* -1 * (1 ./ tau))); % Exponential decay
nbkn = @(x, n) [NaN(n, 1); x(1:end-n)]; % N-back shifter
fillbf = @(x, y) [y x(:, 1:end-1)]; % Fill memory buffer

%% Transition Model
tmodel = @(x) exp((-1 ./ (2 .* Sd.^2)) .* abs(mod(x + 90, 180) - 90).^2);

%% Initialize
history = uniform;
decision = zeros(trials, 1);
posterior_s = decision;
buffer = repmat(uniform', [1, nback]);
d = tmodel(ch - 90);
d = d ./ sum(d, 2);
tr_prior = Pt .* d + (1 - Pt) .* uniform;

%% Main Loop
for k = 1:trials
    % Encoding
    if numel(Se) > 1
        theta_k = theta(k);
        if meas ~= 0
            theta_k = circ_vmrnd(deg2rad(2 * theta(k)), dSd2k(meas(k)), 1);
            theta_k = mod(rad2deg(theta_k) / 2, 180);
        end
        encoding = chresp(theta_k, Se(k), ch) .* history;
    else
        theta_k = theta(k);
        if meas ~= 0
            theta_k = circ_vmrnd(deg2rad(2 * theta(k)), dSd2k(meas), 1);
            theta_k = mod(rad2deg(theta_k) / 2, 180);
        end
        encoding = chresp(theta_k, Se, ch) .* history;
    end

    % Decoding
    decoding = decoder(encoding', ch');
    [~, s0] = circ_std(deg2rad(2 * ch'), encoding');
    posterior_s(k) = rad2deg(s0) / 2;

    % Determine if one or two distributions should be used
    % Create distributions for current and previous stimuli
    cc = chresp(theta_k, Se, ch);
    prev = nbkn(theta', 1);prev = prev(2:end);
    prev_dist = chresp(prev, Sd, ch);

    % Cost functions for one and two distributions
    costfun_two_dist = @(x) sum((cc - (x(1) * cc + (1 - x(1)) * prev_dist)).^2, 'omitnan');
    costfun_one_dist = @(x) sum((cc - x(1) * cc).^2, 'omitnan');

    % Optimization
    options = optimoptions('fmincon', 'Display', 'off');
    par0_two = [0.5]; % Initial guess for weight between distributions
    par0_one = [1]; % Initial guess for single distribution weight
    lb = 0;
    ub = 1;
    weight_two = fmincon(costfun_two_dist, par0_two, [], [], [], [], lb, ub, [], options);
    weight_one = fmincon(costfun_one_dist, par0_one, [], [], [], [], lb, ub, [], options);

    % Decision based on cost function comparison
    if costfun_two_dist(weight_two) < costfun_one_dist(weight_one)
        % Use two distributions
        prior = weight_two * cc + (1 - weight_two) * prev_dist;
    else
        % Use one distribution
        prior = cc;
    end

    % Update History and Buffer
    buffer = fillbf(buffer, prior');
    wdecay = decay(1:nback, tau);
    history = sumnorm(sum(wdecay .* buffer, 2)');

    % Save Decoded Decision
    decision(k, :) = decoding;
end

%% Output
model.delta = sdp_acute_ang(nbkn(theta, 1) - theta);
model.error = sdp_acute_ang(decision - theta);
model.decoded = decision;
model.theta = theta;
model.sigma = posterior_s;
model.bias = grpstats(model.error, model.delta, 'mean')';
model.scatter = grpstats(model.error, model.delta, 'std')';

if nback > 0
    tmp = table(ones(size(theta)), theta, model.error, ...
        'VariableNames', {'obs', 'theta', 'errorc'});
    deltabk = sdp_nbk_deltas(tmp, nback);
    y_mv = zeros(181, nback);
    for k = 1:nback
        y_mv(:, k) = sdp_mvav(deltabk(:, k), tmp.errorc, bins, [], [], 'mean');
    end
    model.sdpbk = y_mv;
    sdp = y_mv(:);
else
    sdp = model.bias;
end
end
