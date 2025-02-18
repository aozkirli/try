clearvars; close all;
delete(gcp); parpool

%% Load Data
subjs = 1:32;
rawdataIndv = cell(32,1);
dataIndv = cell(32,2);
% Loop through subjects
cd ../..
for s = subjs
    %% Load Data
    rawdata = []; nbackdata = [];
    for b = 1:15
        load(['E_subj' num2str(s) '_b' num2str(b) '.mat'])
        rawdata(length(array)*(b-1)+1:length(array)*b,:,1) = array;
    end
    rawdataIndv{s} = rawdata;
    %% Cardinal Bias Correction
    currdata = rawdata;
    % Cardinal correction using GP regression
%     currdata(:,4) = bound_180( currdata(:,4) - predict( gprMdl{s}, currdata(:,2) ) );
    % Mean centering instead of cardinal correction
    currdata(:,4) = bound_180( currdata(:,4) - c_mean(currdata(:,4)) );
    % Response correction
    currdata(:,3) = bound_360( currdata(:,2) + currdata(:,4) );
    %% Data Preprocessing & Labeling
    % Calculate relative stimulus, response and error on nth previous trial
    nBack = 1;
    for n = 1:nBack
        currdata(1+n:end,5+2*n-1) = bound_180( currdata(1:end-n,2) - currdata(1+n:end,2) );
        currdata(1+n:end,5+2*n) = bound_180( currdata(1:end-n,3) - currdata(1+n:end,2) );
    end
    currdata(2:end,end+1) = currdata(1:end-1,4);
    % Exclude first n trials
    for b = 1:15
        nbackdata( (length(array)-nBack)*(b-1)+1:(length(array)-nBack)*b,:) = currdata( length(array)*(b-1)+1+nBack:length(array)*b, : );
    end
    % Correct some labels
    if s==4 || s==8
        nbackdata( nbackdata(:,1) == 180, 1 ) = -180;
    end
    %% Outlier correction
    % Exclude trials in which error was further than 2.5 s.d. away from the mean
    outliers = abs( nbackdata(:,4) - c_mean(nbackdata(:,4)) ) > 2.5*c_std(nbackdata(:,4));
    data = nbackdata( ~outliers & ~([false(nBack,1); outliers(1:end-nBack)]), : );
    nExcludedTrials(s) = length(nbackdata) - length(data);
    %%%%%%%% IMPORTANT %%%%%%%%%
    % in radian
    data(:,[1:4,6:end]) = deg2rad(data(:,[1:4,6:end]));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Finalize data
    dataIndv{s} = data;
end
cd([pwd,'\MCMC\woCorrection'])


%% Run MCMC
for k = 1:3
    parfor repeat_number = 1:4
        [samples, prop_cov] = DoG_MCMC( dataIndv, repeat_number, k );
    end
end


%% :: FUNCTION ::

%% :: FUNCTION :: DoG_MCMC
function [samples, prop_cov] = DoG_MCMC( dataIndv, repeat_number, k )

switch k
    case 1
        output_file = strcat('DoG_woCorrection_S_r',num2str(repeat_number),'.mat');
    case 2
        output_file = strcat('DoG_woCorrection_R_r',num2str(repeat_number),'.mat');
    case 3
        output_file = strcat('DoG_woCorrection_SR_r',num2str(repeat_number),'.mat');        
end
if ~isempty(dir(output_file))
    error('The file name already exists.')
end

n_samples = 2000; % 2000
thinning = 1000; % 1000
lambda1 = 1; % 0 to 1
lambda2 = 3; % 0 >
prop_factor = 1e-5;
switch k
    case {1,2}
        n_params = 3;
        prop_update = 200; % 200
        update_thinning = 50; % 50
        factor_update = 200; % 200
    case 3
        n_params = 5;
        prop_update = 800; % 200
        update_thinning = 200; % 50
        factor_update = 800; % 200
end

% dataAll{s} = [ theta_n-1 - theta_n, stim_theta, resp_theta, resp_theta - stim_theta, rt ];
n_subjs = length( dataIndv );


switch k
    case {1,2}
        mu_temp = [ 0; 30; 10 ]; % [ aStim; peakLocStim; ySigma ];
        std_temp = [ 1; 5; 2 ]; % population sigma
    case 3
        mu_temp = [ 0; 30; 0; 30; 10 ]; % [ aStim; peakLocStim; aResp; peakLocResp; ySigma ];
        std_temp = [ 1; 5; 1; 5; 2 ]; % population sigma
end
mu_temp = repmat( mu_temp, n_subjs + 1, 1 ); % N subjects + population mu

curr_sample = [ mu_temp; std_temp ];
switch k
    case {1,2,3}
        curr_sample = deg2rad( curr_sample );
end


var_temp = deg2rad(std_temp).^2;
var_temp = repmat( var_temp, n_subjs + 1, 1 );
var_std = 0.25*deg2rad(std_temp).^2;
var_temp = [ var_temp; var_std ];

prop_cov = diag( var_temp );

%
log_prior_curr = cal_log_prior(curr_sample, n_params, n_subjs, k);
log_like_curr = cal_log_lik(curr_sample, dataIndv, n_params, n_subjs, k);
log_post_curr = log_like_curr + log_prior_curr;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
samples = zeros(length(curr_sample), n_samples);
samples_update = zeros(length(curr_sample), prop_update);

update_number = 1;
update_times = 1;
update_accept = 0;
sample_number = 1;
accept = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a = GetSecs;

for i = 1:n_samples*thinning
    
    proposal_sample = sampling_proposal(curr_sample, prop_factor*prop_cov);
    
    log_prior_prop = cal_log_prior(proposal_sample, n_params, n_subjs, k);
    log_like_prop = cal_log_lik(proposal_sample, dataIndv, n_params, n_subjs, k);
    log_post_prop = log_prior_prop + log_like_prop;
    
    acc = log_post_prop - log_post_curr;
    if isnan(log_post_prop)
        acc = -inf;
    end
    
    if log(rand) <= acc
        curr_sample = proposal_sample;
        log_post_curr = log_post_prop;
        accept = accept + 1;
        update_accept = update_accept + 1;
    end
    
    if i>0 && ~rem(i,thinning)
        samples(:,sample_number) = curr_sample;
        sample_number = sample_number + 1;
        disp(round([sample_number-1,100*accept/i]))
    end
        
    gamma1 = (1/(update_times)^(lambda1));
    gamma2 = lambda2*gamma1;
        
    if ~rem(i, update_thinning)
        samples_update(:,update_number) = curr_sample;
        if ~rem(i,update_thinning*prop_update)
            prop_cov = prop_cov + gamma1 * (cov(samples_update') - prop_cov);
            update_number = 1;
            update_times = update_times + 1;
        else
            update_number = update_number + 1;
        end
    end
    
    if ~rem(i,factor_update)
        log_prop_factor = log(prop_factor) + gamma2 * ((update_accept / factor_update) - 0.234);
        prop_factor = exp(log_prop_factor)
        disp(update_accept / factor_update)
        update_accept = 0;
    end
    
end
disp(round([sample_number,100*accept/i]))
b = GetSecs;
duration = b-a

switch k
    case {1,2,3}
        samples = rad2deg( samples );
end

save (output_file, 'samples');


end


%% :: FUNCTION :: cal_log_prior
function log_prior_all = cal_log_prior(sample, n_params, n_subjs, k )

sample_matrix = reshape( sample, n_params, length(sample)/n_params );

% prior
mu_vector = sample_matrix( :, end-1 );
mu_vector = repmat( mu_vector, n_subjs, 1 );

std_prior = sample_matrix( :, end );
std_vector = repmat( std_prior, n_subjs, 1 );

log_prior = sum(log(circ_vmpdf( sample(1:n_params*n_subjs), mu_vector, 1./std_vector.^2 )));

% hyper-prior
switch k
    case {1,2}
        gammaMode = deg2rad( [ .5; 2.5; 1 ] ); % sd(y)/2
    case 3
        gammaMode = deg2rad( [ .5; 2.5; .5; 2.5; 1 ] ); % sd(y)/2
end
gammaSD = 4*gammaMode; % gammaMode % 2*sd(y)

gammaScale = 2*gammaSD.^2 ./ (gammaMode + sqrt( gammaMode.^2 + 4*gammaSD.^2 )); % gammaRate = 1./gammaScale;
gammaShape = 1 + gammaMode./gammaScale;

log_hyper_prior = sum(log(gampdf( std_prior, gammaShape, gammaScale )));

% finalize
switch k
    case {1,2}
        flagImplausible = sum( sample_matrix( 2, 1:end-1 ) <= pi/36 ) > 0 || sum( sample_matrix( 2, 1:end-1 ) >= pi/2 ) > 0 || ...
            sum( sample_matrix( 3, 1:end-1 ) <= eps ) > 0;
    case 3
        flagImplausible = sum(sum( sample_matrix( [2 4], 1:end-1 ) <= pi/36 )) > 0 || sum(sum( sample_matrix( [2 4], 1:end-1 ) >= pi/2 )) > 0 || ...
            sum( sample_matrix( 5, 1:end-1 ) <= eps ) > 0;
end
if flagImplausible
    log_prior_all = -inf;
else
    log_prior_all = log_prior + log_hyper_prior;
end

end



%% :: FUNCTION :: cal_log_lik
function sum_log_likeli = cal_log_lik( curr_sample, dataIndv, n_params, n_subjs, k )

log_likeli_all = zeros( n_subjs, 1 );

params_matrix = reshape( curr_sample, n_params, n_subjs + 2 );

for s = 1:n_subjs
    
    params = params_matrix( :, s );
    data = dataIndv{s};
    
    log_likeli_all(s) = cal_log_lik_subj( params, data, k );
    
end
sum_log_likeli = sum( log_likeli_all );

end


%% :: FUNCTION :: cal_log_likeli_subj
function log_likeli_subj = cal_log_lik_subj( params_subj, data, k )

switch k
    
    case 1
        
        a = params_subj(1);
        peakLoc = params_subj(2);
        ySigma = params_subj(3);
        
        log_likeli_subj = sum(log(circ_vmpdf( data(:,4), DoG( data(:,6), a, peakLoc ), 1./ySigma.^2 )));
        
    case 2
        
        a = params_subj(1);
        peakLoc = params_subj(2);
        ySigma = params_subj(3);
        
        log_likeli_subj = sum(log(circ_vmpdf( data(:,4), DoG( data(:,7), a, peakLoc ), 1./ySigma.^2 )));
        
    case 3
        
        aStim = params_subj(1);
        peakLocStim = params_subj(2);
        aResp = params_subj(3);
        peakLocResp = params_subj(4);
        ySigma = params_subj(5);

        log_likeli_subj = sum(log(circ_vmpdf( data(:,4), DoG( data(:,1), aStim, peakLocStim ) + DoG( data(:,7), aResp, peakLocResp ), 1./ySigma.^2 )));

end
        
end


%% :: FUNCTION :: DoG
function y = DoG( x, a, peakLoc )

w = 1./(sqrt(2)*peakLoc);
y = x.*a.*w.*sqrt(2)./exp(-.5).*exp(-1.*(w.*x).^2);

end

%% :: FUNCTION :: sampling_proposal
function sample_all = sampling_proposal(curr_vector, prop_cov)

sample = mvnrnd( curr_vector', prop_cov );
sample_all = sample';
end


%% :: FUNCTION :: circ_vmpdf
function [p, alpha] = circ_vmpdf(alpha, thetahat, kappa)

if nargin < 1 || isempty(alpha)
    alpha = linspace(0, 2*pi, 101)';
    alpha = alpha(1:end-1);
end
if nargin < 3
    kappa = 1;
end
if nargin < 2
    thetahat = 0;
end

% evaluate pdf
p = exp( kappa.*(cos(alpha-thetahat)-1) ) ./ (2*pi.*besseli(0,kappa,1));

end


%% :: FUNCTION :: c_mean
function mu = c_mean(alpha, w, dim)

if nargin < 3
  dim = find(size(alpha) > 1, 1, 'first');
  if isempty(dim)
    dim = 1;
  end
end

if nargin < 2 || isempty(w)
  % if no specific weighting has been specified
  % assume no binning has taken place
	w = ones(size(alpha));
else
  if size(w,2) ~= size(alpha,2) || size(w,1) ~= size(alpha,1) 
    error('Input dimensions do not match');
  end 
end

% deg2rad
alpha = deg2rad(alpha);

% compute weighted sum of cos and sin of angles
r = sum(w.*exp(1i*alpha),dim,'omitnan');

% obtain mean by
mu = angle(r);

% rad2deg
mu = rad2deg(mu);

end

%% :: FUNCTION :: c_std
function s0 = c_std(alpha, w, d, dim)

if nargin < 4
  dim = find(size(alpha) > 1, 1, 'first');
  if isempty(dim)
    dim = 1;
  end  
end

if nargin < 3 || isempty(d)
  % per default do not apply correct for binned data
  d = 0;
end

if nargin < 2 || isempty(w)
  % if no specific weighting has been specified
  % assume no binning has taken place
	w = ones(size(alpha));
else
  if size(w,2) ~= size(alpha,2) || size(w,1) ~= size(alpha,1) 
    error('Input dimensions do not match');
  end 
end

% deg2rad
alpha = deg2rad(alpha);

% compute mean resultant vector length
r = circ_r(alpha,w,d,dim);

s0 = sqrt(-2*log(r));   % 26.21

% rad2deg
s0 = rad2deg(s0);

end

%% :: FUNCTION :: circ_r
function r = circ_r(alpha, w, d, dim)

if nargin < 4
  dim = find(size(alpha) > 1, 1, 'first');
  if isempty(dim)
    dim = 1;
  end
end

if nargin < 2 || isempty(w) 
  % if no specific weighting has been specified
  % assume no binning has taken place
	w = ones(size(alpha));
else
  if size(w,2) ~= size(alpha,2) || size(w,1) ~= size(alpha,1) 
    error('Input dimensions do not match');
  end 
end

if nargin < 3 || isempty(d)
  % per default do not apply correct for binned data
  d = 0;
end

% compute weighted sum of cos and sin of angles
r = sum(w.*exp(1i*alpha),dim);

% obtain length 
r = abs(r)./sum(w,dim);

% for data with known spacing, apply correction factor to correct for bias
% in the estimation of r (see Zar, p. 601, equ. 26.16)
if d ~= 0
  c = d/2/sin(d/2);
  r = c*r;
end
end

%% :: FUNCTION :: bound_180
function vec = bound_180( vec )

vec( vec > 180 ) = vec( vec > 180 ) - 360;
vec( vec < -180 ) = vec( vec < -180 ) + 360;

end

%% :: FUNCTION :: bound_360
function vec = bound_360( vec )

vec( vec > 360 ) = vec( vec > 360 ) - 360;
vec( vec < 0 ) = vec( vec < 0 ) + 360;

end
