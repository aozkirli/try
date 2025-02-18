clearvars

%%
subjs = 1:32;
M = 4;
N = 1000;
    
filename = 'AR';
n_params = 7;

samples_all = [];
for repeat_number = 1:M
    load(strcat('DoG_',filename,'_r',num2str(repeat_number),'.mat'));
    samples_all(:,end+1:end+N) = samples(:, N+1:end);

    thetahat(:,repeat_number) = mean( samples(:, N+1:end), 2 );
    varhat(:,repeat_number) = var( samples(:, N+1:end), 0, 2 );
end

% subjects
for s = 1:length(subjs)

    doi = samples_all( (s-1)*n_params+1 : s*n_params, : )';

    [n,d] = size(doi);
    bw = std(doi)*(4/((d+2)*n))^(1/(d+4));

    mvkde = @(params) -mvksdensity(doi,params,'Bandwidth',bw);
    params0 = mean(doi);
    options = optimset('MaxIter',1e+5,'MaxFunEvals',1e+5);

    par.modeIndv( s, 1:n_params ) = fminsearch(@(params) mvkde(params), params0, options);
    
end

% population
for i = 1:n_params

    doi = samples_all( end-2*n_params+i, : )';
    
    % mode
    kde = @(params) -ksdensity( doi, params );
    params0 = mean(doi);

    par.modeGroup(i) = fminsearch(@(params) kde(params), params0, options);

    % hdi
    nx = 10000; % even number
    xi = linspace( mean(doi)-5*std(doi), mean(doi)+5*std(doi), nx );
    f = ksdensity( doi, xi );

    pci = 95; params0 = f(3000);
    fval = fminsearchbnd(@(params) hdi( params, xi, f, pci ), params0, eps, max(f), options);

    [~,xl] = min(abs( f( 1 : length(f)/2 ) - fval ));

    [~,minI] = min(abs( f( length(f)/2+1 : end ) - fval ));
    xu = length(f)/2 + minI;

    par.hdi95(i,:) = xi([xl xu]);

end

B = N/(M-1)*sum( (thetahat-mean(thetahat,2)).^2, 2 );
W = 1/M*sum(varhat,2);
Vhat = (N-1)/N*W + (M+1)/(M*N)*B;
par.Rhat = Vhat./W;

save(strcat('par_',filename,'.mat'),'par','samples_all');


%% :: FUNCTION ::
function output = hdi( fval, xi, f, pci )

[~,xl] = min(abs( f( 1 : length(f)/2 ) - fval ));

[~,minI] = min(abs( f( length(f)/2+1 : end ) - fval ));
xu = length(f)/2 + minI;

output = abs( pci/100 - sum(f(xl:xu))*(range(xi)/(length(xi)-1)) ); 

end