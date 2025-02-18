function fits = fit_dog(cfg)

if size(cfg.data,2) == 1 % Moving average
    x = [-90:1:90]';
    y = cfg.data;
elseif size(cfg.data,2) == 2 % Datapoints
    x = cfg.data(:,1);
    y = cfg.data(:,2);
end


if cfg.fixedwidth == false
    
    Model   = @(a,w,x) x.*a.*w.*(sqrt(2)./exp(-0.5)).*exp(-(w.*x).^2);
    
    % determine range for starting positions in parameter space
    upper_a = 10;
    lower_a = 0;

    lower_w = 0.02;
    upper_w = 0.07;

    % determine random starting positions within specified range
    start_a = (upper_a-lower_a).*rand(1,cfg.fittingsteps) + lower_a;
    start_w = (upper_w-lower_w).*rand(1,cfg.fittingsteps) + lower_w;
    
elseif cfg.fixedwidth == true
    
    width   = cfg.width;
    Model   = @(a,x) x.*a.*width.*(sqrt(2)./exp(-0.5)).*exp(-(width.*x).^2);
    
    % determine range for starting positions in parameter space
    upper_a = 10;
    lower_a = 0;
    
    % determine random starting positions within specified range
    start_a = (upper_a-lower_a).*rand(1,cfg.fittingsteps) + lower_a;

end

% Allocating space for the model fits and creating bookkeeping variables
fits.fitobjects = cell(cfg.fittingsteps,1);
fits.gof        = cell(cfg.fittingsteps,1);
sse             = Inf;
curr_min        = Inf;

for i = 1:cfg.fittingsteps
    
    if cfg.fixedwidth == true
        [fitobject, gof] = fit(x,y,Model,'Robust', 'LAR','MaxIter',10000,'MaxFunEvals',10000,'TolFun',10^(-8),'TolX',10^(-8),'StartPoint',start_a(1,i));
    elseif cfg.fixedwidth == false
        [fitobject, gof] = fit(x,y,Model,'Robust', 'LAR','MaxIter',10000,'MaxFunEvals',10000,'TolFun',10^(-8),'TolX',10^(-8),'StartPoint',[start_a(1,i) start_w(1,i)],'Lower',[-Inf lower_w],'Upper',[+Inf upper_w]);
    end
    
    
    fits.fitobjects{i,1}.object = fitobject;
    fits.gof{i,1}               = gof;
    
    if gof.sse < sse
        sse         = gof.sse;
        curr_min    = i;
    end
    
end

fits.bestfit = fits.fitobjects{curr_min,1}.object;
fits.min = curr_min;

if cfg.fixedwidth == true
    fits.coeffs = [coeffvalues(fits.bestfit) cfg.width];
elseif  cfg.fixedwidth == false
    fits.coeffs = coeffvalues(fits.bestfit);
end

end