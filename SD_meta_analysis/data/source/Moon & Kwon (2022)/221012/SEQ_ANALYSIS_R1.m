clearvars;
% close all;
rng default
plot_or_not = true;
DoG = @(x,a,w) x.*a.*w.*sqrt(2)./exp(-.5).*exp(-1.*(w.*x).^2);

%% Load Data
rawdataIndv = cell(32,1);
dataIndv = cell(32,1);
dataCombined = [];
load SEQ_GPRmodels
subjs = 1:32;
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
    currdata(:,4) = bound_180( currdata(:,4) - predict( gprMdl{s}, currdata(:,2) ) );
    currdata(:,3) = bound_360( currdata(:,2) + currdata(:,4) );
    %% Data Preprocessing & Labeling
    % Calculate relative stimulus, response and error on nth previous trial
    nBack = 1;
    for n = 1:nBack
        currdata(1+n:end,5+2*n-1) = bound_180( currdata(1:end-n,2) - currdata(1+n:end,2) );
        currdata(1+n:end,5+2*n) = bound_180( currdata(1:end-n,3) - currdata(1+n:end,2) );
    end
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
    % Some statistics
    SRcorr(s) = circ_corrcc(deg2rad(data(:,2)),deg2rad(data(:,3)));
    % Finalize data
    dataIndv{s} = data;
    dataCombined(end+1:end+length(data),:) = data;
end


%% Analysis & Plot
% For all analyses
binWidth = 10;
xBin = -90:binWidth:90;
% For marginal bias plot
xMoving = -90:90;
MarginalStimBias = nan(length(subjs),length(xBin));
MarginalRespBias = nan(length(subjs),length(xMoving));
errorStimSD = nan(length(subjs),length(xBin));
numTrialsThreshold_MRB = 10;
% For joint bias map
map = nan(length(xBin),length(xBin),length(subjs));
indvmap = nan(length(xBin),length(xBin),length(subjs));
mapHowmany = nan(length(xBin),length(xBin),length(subjs));
isShortIndv = cell(1,length(subjs));
numEmpty = zeros(1,length(subjs));
numShort = zeros(1,length(subjs));
mapAll = nan(length(xBin));
numEmptyAll = 0;
numShortAll = 0;
[X,Y] = meshgrid(1:length(xBin)); X = X(:); Y = Y(:);
numSubjThreshold = 14;
numTrialsThreshold = 5;
% For conditional bias plot
givenResp = -20:binWidth:20;
StimBiasGivenResp = nan(length(givenResp),length(xBin),length(subjs));
RespBiasGivenStim = nan(length(givenResp),length(xBin),length(subjs));
meanSBGR = nan(length(givenResp),length(xBin));
ciSBGR = nan(length(givenResp),length(xBin));
meanRBGS = nan(length(givenResp),length(xBin));
ciRBGS = nan(length(givenResp),length(xBin));
% For DoG prediction
nModel = 3;
DoGmap = zeros(length(xBin),length(xBin),nModel+1);
mdlString = cell(1,nModel+1);
AIC = zeros(length(subjs),3);
AICc = zeros(length(subjs),3);
resi_R = cell(length(subjs),1);
resi_S = cell(length(subjs),1);
ResidualStimBias = nan(length(subjs),length(xBin));
ResidualRespBias = nan(length(subjs),length(xMoving));
ResidualStimBiasCombined = zeros(1,length(xBin));
ResidualRespBiasCombined = zeros(1,length(xMoving));

% Subject Loop
for s = subjs
    % Load either empirical or simulation data
    data = dataIndv{s};
    mapdata = dataCombined;

    %% Marginal Bias
    % Marginal stimulus bias
    for j = 1:length(xBin)
        l = data(:,6) == xBin(j);
        if sum(l) > 0
            MarginalStimBias(s,j) = c_mean(data(l, 4));
            errorStimSD(s,j) = std( data(data(:,1) == xBin(j), 4) - MarginalStimBias(s,j) );
        end
    end
    % Marginal response bias
    for j = 1:length(xMoving)
        l = abs( data(:,7) - xMoving(j) ) <= binWidth/2;
        if sum(l) > numTrialsThreshold_MRB
            MarginalRespBias(s,j) = c_mean(data(l,4));
        end
    end
    
    %% Joint Bias
    % Individual
    for i = 1:length(xBin)
        for j = 1:length(xBin)
            l = ( data(:,6)==xBin(i) ) & ( data(:,7)>=xBin(j)-binWidth/2 & data(:,7)<xBin(j)+binWidth/2 );
            if sum(l) == 0
                numEmpty(s) = numEmpty(s) + 1;
            elseif sum(l) < numTrialsThreshold
                numShort(s) = numShort(s) + 1;
                map(i,j,s) = c_mean(data(l,4));
                mapHowmany(i,j,s) = sum(l);
            else
                indvmap(i,j,s) = c_mean(data(l,4));
                map(i,j,s) = c_mean(data(l,4));
                mapHowmany(i,j,s) = sum(l);
            end
        end
    end
    isShortIndv{s} = rot90(isnan(indvmap(:,:,s)));

    % Mean & Combined
    if s == subjs(end)
        indvmap = rot90(indvmap);
        map = rot90(map);
        mapHowmany = rot90(mapHowmany);
        meanmap = c_mean(map,[],3);
        DoGmap(:,:,end) = meanmap; mdlString{end} = 'Data';
        isShort = sum(~isnan(map),3)<numSubjThreshold;
        for i = 1:length(xBin)
            for j = 1:length(xBin)
                l = ( mapdata(:,6)==xBin(i) ) & ( mapdata(:,7)>=xBin(j)-binWidth/2 & mapdata(:,7)<xBin(j)+binWidth/2 );
                if sum(l) == 0
                    numEmptyAll = numEmptyAll + 1;
                elseif sum(l) < numSubjThreshold
                    numShortAll = numShortAll + 1;
                else
                    mapAll(i,j) = c_mean(mapdata(l,4));
                end
            end
        end
        mapAll = rot90(mapAll);
        isShortAll = isnan(mapAll);
    end

    %% Conditional Bias
    % Stimulus Bias Given Response
    for i = 1:length(givenResp)
        for j = 1:length(xBin)
            l = (data(:,6)==xBin(j)) & (data(:,7)>=givenResp(i)-binWidth/2 & data(:,7)<givenResp(i)+binWidth/2);
            if sum(l) >= numTrialsThreshold
                StimBiasGivenResp(i,j,s) = c_mean(data(l,4));
            end
            if s == subjs(end)
                if sum(~isnan(StimBiasGivenResp(i,j,:)),3) >= numSubjThreshold
                    meanSBGR(i,j) = c_mean(StimBiasGivenResp(i,j,:),3);
                    ciSBGR(i,j) = tinv(0.975,sum(~isnan(StimBiasGivenResp(i,j,:)))-1).*c_std(StimBiasGivenResp(i,j,:),~isnan(StimBiasGivenResp(i,j,:)))./sqrt(sum(~isnan(StimBiasGivenResp(i,j,:)),3));
                end
            end
        end
    end

    % Response Bias Given Stimulus
    for i = 1:length(givenResp)
        for j = 1:length(xBin)
            l = (data(:,6)==givenResp(i)) & (data(:,7)>=xBin(j)-binWidth/2 & data(:,7)<xBin(j)+binWidth/2);
            if sum(l) >= numTrialsThreshold
                RespBiasGivenStim(i,j,s) = c_mean(data(l,4));
            end
            if s == subjs(end)
                if sum(~isnan(RespBiasGivenStim(i,j,:)),3) >= numSubjThreshold
                    meanRBGS(i,j) = c_mean(RespBiasGivenStim(i,j,:));
                    ciRBGS(i,j) = tinv(0.975,sum(~isnan(RespBiasGivenStim(i,j,:)))-1).*c_std(RespBiasGivenStim(i,j,:),~isnan(RespBiasGivenStim(i,j,:)))./sqrt(sum(~isnan(RespBiasGivenStim(i,j,:))));
                end
            end
        end
    end

    %% DoG Curve Behavior Prediction -- static, descriptive DoG
    for mdl = 1:nModel
        switch mdl
            case 1   
                load par_S.mat
                nParams = 3;
                % Model prediction
                predictedBias = DoG( data(:,6), par.modeIndv(s,1), 1/(sqrt(2)*par.modeIndv(s,2)) );
                noiseSigma = par.modeIndv(s,3);
                % Residual of the model
                resi_S{s} = bound_180( data(:,4) - predictedBias );
                for j = 1:length(xMoving)
                    l = abs( data(:,7) - xMoving(j) ) <= binWidth/2;
                    if sum(l) > numTrialsThreshold_MRB
                        ResidualRespBias(s,j) = c_mean(resi_S{s}(l));
                    end
                end
            case 2
                load par_R.mat
                nParams = 3;
                % Model prediction
                predictedBias = DoG( data(:,7), par.modeIndv(s,1), 1/(sqrt(2)*par.modeIndv(s,2)) );
                noiseSigma = par.modeIndv(s,3);
                % Residual of the model
                resi_R{s} = bound_180( data(:,4) - predictedBias );
                for j = 1:length(xBin)
                    l = data(:,6) == xBin(j);
                    if sum(l) > 0
                        ResidualStimBias(s,j) = c_mean(resi_R{s}(l));
                    end
                end
            case 3
                load par_SR.mat
                nParams = 5;
                % Model prediction
                predictedBias = DoG( data(:,6), par.modeIndv(s,1), 1/(sqrt(2)*par.modeIndv(s,2)) ) + DoG( data(:,7), par.modeIndv(s,3), 1/(sqrt(2)*par.modeIndv(s,4)) );
                noiseSigma = par.modeIndv(s,5);
        end
        AIC(s,mdl) = 2*nParams-2*sum(log( circ_vmpdf( deg2rad(data(:,4)), deg2rad(predictedBias), 1/deg2rad(noiseSigma)^2 ) ));
        AICc(s,mdl) = AIC(s,mdl) + ( 2*length(par.modeGroup)^2 + 2*length(par.modeGroup) ) / ( length(subjs) - length(par.modeGroup) - 1 );
        % Population level
        if s == subjs(end)
            switch mdl
                case 1 % Bias to (or away from) stimulus
                    mdlString{mdl} = 'Stimulus';
                    tempmap = DoG( xBin(X), par.modeGroup(1), 1/(sqrt(2)*par.modeGroup(2)) );
                case 2 % Bias to (or away from) response
                    mdlString{mdl} = 'Response';
                    tempmap = DoG( fliplr(xBin(Y)), par.modeGroup(1), 1/(sqrt(2)*par.modeGroup(2)) );
                case 3 % Bias to (or away from) stimulus & Bias to (or away from) response
                    mdlString{mdl} = 'Stimulus & Response';
                    tempmap = DoG( xBin(X), par.modeGroup(1), 1/(sqrt(2)*par.modeGroup(2)) ) + DoG( fliplr(xBin(Y)), par.modeGroup(3), 1/(sqrt(2)*par.modeGroup(4)) );
            end
            DoGmap(:,:,mdl) = reshape(tempmap,length(xBin),length(xBin));
            % Residual of the models
            load par_R.mat
            resi_R_Combined = bound_180( dataCombined(:,4) - DoG( dataCombined(:,7), par.modeGroup(1), 1/(sqrt(2)*par.modeGroup(2)) ) );
            for j = 1:length(xBin)
                l = dataCombined(:,6) == xBin(j);
                if sum(l) > 0
                    ResidualStimBiasCombined(j) = c_mean(resi_R_Combined(l));
                end
            end
            load par_S.mat
            resi_S_Combined = bound_180( dataCombined(:,4) - DoG( dataCombined(:,6), par.modeGroup(1), 1/(sqrt(2)*par.modeGroup(2)) ) );
            for j = 1:length(xMoving)
                l = abs( dataCombined(:,7) - xMoving(j) ) <= binWidth/2;
                if sum(l) > numTrialsThreshold_MRB
                    ResidualRespBiasCombined(j) = c_mean(resi_S_Combined(l));
                end
            end
        end
    end

end


%% :: PLOT ::

%% :: PLOT :: General Settings
set(groot,'DefaultFigureColor', 'w')
set(groot,'DefaultAxesLineWidth', 1)
set(groot,'DefaultAxesXColor', 'k')
set(groot,'DefaultAxesYColor', 'k')
set(groot,'DefaultAxesFontUnits', 'points')
set(groot,'DefaultAxesFontSize', 8/.45) 
set(groot,'DefaultAxesFontName', 'Helvetica')
set(groot,'DefaultLineLineWidth', 2/3)
set(groot,'DefaultTextFontUnits', 'Points')
set(groot,'DefaultTextFontSize', 8/.45)
set(groot,'DefaultTextFontName', 'Helvetica')
set(groot,'DefaultAxesBox', 'off')
set(groot,'DefaultAxesTickLength', [0.025 0.025]);
set(groot,'DefaultAxesTickDir','out');
set(groot,'DefaultAxesTickDirMode','manual');
set(groot,'DefaultAxesLabelFontSizeMultiplier',1);

load bwr.mat

if plot_or_not
    % Axis limits
    XLimCB = [-45 45]; YLimCB = [-10 10]; XTickCB = [-40 40]; YTickCB = [-8 8];
    XLimMB = [-100 100]; YLimMB = YLimCB; XTickMB = [-90 90]; YTickMB = YTickCB; nXTick = 7;
    
    %% :: Figure 1BC :: Marginal Bias
    figure; 
    % Marginal stimulus bias
    subplot(1,2,1); hold on
    plot([XLimMB(1) 0; XTickMB(2) 0],[0 YLimMB(1); 0 YTickMB(2)],'k--','linewidth',2/3);
    meanMSB = c_mean(MarginalStimBias);
    ciMSB = tinv(0.975,sum(~isnan(MarginalStimBias))-1).*c_std(MarginalStimBias,~isnan(MarginalStimBias))./sqrt(sum(~isnan(MarginalStimBias)));
    shadedErrorBar(xBin,meanMSB,ciMSB,'lineprops',{'-','color','k','linewidth',2});
    % Setting
    ax = gca; pbaspect(ax, [1 1.1 1])
    ax.FontSize = 16;
    ax.XLabel.String = ['Relative direction of' newline 'previous stimulus (\circ)'];
    ax.YLabel.String = 'Response error (\circ)';
    ax.XLim = XLimMB; ax.YLim = YLimMB;
    ax.XTick = linspace(XTickMB(1),XTickMB(2),nXTick); ax.YTick = linspace(YTickMB(1),YTickMB(2),5);
    ax.XTickLabel = {'−90','','','0','','','90'}; ax.YTickLabel = {'−8','−4','    0','4','8'};
    offsetAxes(ax,XTickMB,YTickMB)

    % Marginal response bias
    subplot(1,2,2); hold on
    plot([XLimMB(1) 0; XTickMB(2) 0],[0 YLimMB(1); 0 YTickMB(2)],'k--','linewidth',2/3);
    meanMRB = c_mean(MarginalRespBias);
    ciMRB = tinv(0.975,sum(~isnan(MarginalRespBias))-1).*c_std(MarginalRespBias,~isnan(MarginalRespBias))./sqrt(sum(~isnan(MarginalRespBias)));
    shadedErrorBar(xMoving,meanMRB,ciMRB,'lineprops',{'-','color','k','linewidth',2});
    % Setting
    ax = gca; pbaspect(ax, [1 1.1 1])
    ax.FontSize = 16;
    ax.XLabel.String = ['Relative direction of' newline 'previous response (\circ)'];
    ax.XLim = XLimMB; ax.YLim = YLimMB;
    ax.XTick = linspace(XTickMB(1),XTickMB(2),nXTick); ax.YTick = linspace(YTickMB(1),YTickMB(2),5);
    ax.XTickLabel = {'−90','','','0','','','90'}; ax.YTickLabel = {'−8','−4','    0','4','8'};
    offsetAxes(ax,XTickMB,YTickMB)
    
    %% :: Figure S2 :: Marginal Bias -- Representative Individuals
    eg = [13 30 31 16];
    figure; orient('landscape')
    set(gcf,'PaperUnits','centimeters','PaperSize',[34,21]);
    for s = 1:length(eg)
        data = dataIndv{eg(s)};
        
        % Stimulus
        subplot(2,4,s); hold on
        plot([XLimMB(1) 0; XTickMB(2) 0],[0 -25; 0 20],'k--','linewidth',2/3);
        scatter(data(:,6), data(:,4),'markeredgecolor','none','markerfacealpha',.05,'markerfacecolor','k');
        plot(xBin,MarginalStimBias(eg(s),:),'-','color','k','linewidth',2);
        % Setting
        ax = gca; pbaspect(ax, [1 1.1 1])
        ax.FontSize = 16;
        ax.XLabel.String = ['Relative direction of' newline 'previous stimulus (\circ)'];
        if s == 1, ax.YLabel.String = 'Response error (\circ)'; end
        ax.XLim = XLimMB; ax.YLim = [-25 25];
        ax.XTick = linspace(XTickMB(1),XTickMB(2),nXTick); ax.YTick = -20:10:20;
        ax.XTickLabel = {'−90','','','0','','','90'}; ax.YTickLabel = {'    −20','−10','0','10','20'};
        offsetAxes(ax,-90,[-20 20])
        title(['Subject #' num2str(eg(s)) newline '\it{n}\rm = 1,' num2str(size(rawdataIndv{eg(s)},1)-1000) ' trials'],'fontsize',ax.FontSize,'fontweight','normal')
        
        % Response
        subplot(2,4,s+4); hold on
        plot([XLimMB(1) 0; XTickMB(2) 0],[0 -25; 0 20],'k--','linewidth',2/3);
        scatter(data(:,7), data(:,4),'markeredgecolor','none','markerfacealpha',.05,'markerfacecolor','k');
        plot(xMoving,MarginalRespBias(eg(s),:),'-','color','k','linewidth',2);
        % Setting
        ax = gca; pbaspect(ax, [1 1.1 1])
        ax.FontSize = 16;
        ax.XLabel.String = ['Relative direction of' newline 'previous response (\circ)'];
        if s == 1, ax.YLabel.String = 'Response error (\circ)'; end
        ax.XLim = XLimMB; ax.YLim = [-25 25];
        ax.XTick = linspace(XTickMB(1),XTickMB(2),nXTick); ax.YTick = -20:10:20;
        ax.XTickLabel = {'−90','','','0','','','90'}; ax.YTickLabel = {'    −20','−10','0','10','20'};
        offsetAxes(ax,-90,[-20 20])
    end
   
    %% :: Figure S4A :: Joing Bias Map -- All Individuals
    figure; orient('landscape')
    set(gcf,'PaperUnits','centimeters','PaperSize',[45,30]);
    for s = subjs
        subplot(4,8,8*mod(s-1,4)+floor((s-1)/4)+1);
        imagesc(indvmap(:,:,s),[-8 8]); hold on;
        colormap(bwr)
        patch((repmat(X(isShortIndv{s}(:)),1,4)+[-.5 .5 .5 -.5])',(repmat(Y(isShortIndv{s}(:)),1,4)+[-.5 -.5 .5 .5])',[.9 .9 .9],'edgecolor',[.9 .9 .9])
        % Setting
        ax = gca; pbaspect(ax, [1 1 1])
        ax.LineWidth = 2/3;
        ax.FontSize = 20;
        ax.TickDir = 'in'; ax.TickLength = [0 .025];
        ax.XTick = []; ax.YTick = [];
        if sum( s == 4:4:32 ), ax.XTick = 1:9:19; ax.XTickLabel = {'−90','0','90'}; end
        if sum( s == 1:4 ), ax.YTick = 1:9:19; ax.YTickLabel = {'90','0','−90'}; end
    end    
    h = axes(gcf,'visible','off'); h.XLabel.Visible='on'; h.YLabel.Visible='on';
    xlabel(h,'Relative direction of previous stimulus (\circ)','fontsize',20);
    ylabel(h,'Relative direction of previous response (\circ)','fontsize',20);

    %% :: Figure 2A :: Joing Bias Map -- Mean
    figure; 
    clim = [-8 8];
    imagesc(meanmap,clim); hold on;
    colormap(bwr); cbar = colorbar;
    patch((repmat(X(isShort(:)),1,4)+[-.5 .5 .5 -.5])',(repmat(Y(isShort(:)),1,4)+[-.5 -.5 .5 .5])',[.9 .9 .9],'edgecolor',[.9 .9 .9])
    text([4 16],[4 16],'NA','horizontalalignment','center','fontSize',8/.45)
    % Setting
    ax = gca; pbaspect(ax, [1 1 1])
    ax.FontSize = 8/.45;
    ax.XLabel.String = 'Relative direction of previous stimulus (\circ)';
    ax.YLabel.String = 'Relative direction of previous response (\circ)';
    ax.LineWidth = 1; ax.TickDir = 'in'; ax.TickLength = [.01 .025];
    ax.XTick = 1:3:19; ax.XTickLabel = {'−90','−60','−30','0','30','60','90'};
    ax.YTick = 1:3:19; ax.YTickLabel = {'90','60','30','0','−30','−60',' −90'};
    cbar.LineWidth = 1; cbar.Color = 'k'; cbar.FontSize = 8/.45;
    cbar.Ticks = clim(1):range(clim)/4:clim(2); cbar.TickLength = 0;
    cbar.TickLabels = {'−8  ','−4','0','4','8'};
    cbar.Label.String = 'Response error (\circ)';
    cbar.Label.VerticalAlignment = 'middle'; cbar.Label.Rotation = 270;

    %% :: Figure S4B :: Joing Bias Map -- Subject Count
    % The number of subjects
    figure;
    clim = [0 32];
    mapHowmanySubj = sum(~isnan(map),3);
    imagesc(mean(mapHowmanySubj,3),clim);
    colormap(viridis); cbar = colorbar;
    % Setting
    ax = gca; pbaspect(ax, [1 1 1])
    ax.FontSize = 8/.45;
    ax.XLabel.String = 'Relative direction of previous stimulus (\circ)';
    ax.YLabel.String = 'Relative direction of previous response (\circ)';
    ax.LineWidth = 1; ax.TickDir = 'in'; ax.TickLength = [0 .025];
    ax.XTick = 1:3:19; ax.XTickLabel = {'−90','−60','−30','0','30','60','90'};
    ax.YTick = 1:3:19; ax.YTickLabel = {'90','60','30','0','−30','−60','   −90'};
    cbar.LineWidth = 1; cbar.Color = 'k'; cbar.FontSize = 8/.45; 
    cbar.TickLabels = {'0','8','16','24','32     '};
    cbar.Label.String = 'Subject count';
    cbar.Label.VerticalAlignment = 'middle'; cbar.Label.Rotation = 270;
    cbar.Ticks = clim(1):range(clim)/4:clim(2); cbar.TickLength = 0;
    
    %% :: Figure S4CD :: Joint Bias Map -- Trial Count
    % Proportion of trials
    figure;
    subplot(21,20,1:320);
    mapHowmanyTrials = mean( mapHowmany ./ reshape(repmat([300 165 90 45],1,8),1,1,32).*ones(1,1,length(subjs)), 3, 'omitnan');
    clim = [0 .5];
    imagesc(mapHowmanyTrials,clim);
    % Setting
    ax = gca; pbaspect(ax, [1 1 1])
    ax.FontSize = 8/.45;
    ax.YLabel.String = 'Relative direction of previous response (\circ)';
    ax.LineWidth = 1; ax.TickDir = 'in'; ax.TickLength = [0 .025];
    ax.XTick = [];
    ax.YTick = 1:3:19; ax.YTickLabel = {'90','60','30','0','−30','−60','   −90'};
    colormap(ax,cividis); cbar = colorbar;
    cbar.LineWidth = 1; cbar.Color = 'k'; cbar.FontSize = 8/.45; 
    cbar.TickLabels = {'0','0.1','0.2','0.3','0.4','0.5     '};
    cbar.Label.String = 'Proportion of trials';
    cbar.Label.VerticalAlignment = 'middle'; cbar.Label.Rotation = 270;
    cbar.Ticks = clim(1):range(clim)/5:clim(2); cbar.TickLength = 0;
    % text
    mapHowmanyTrials_linear = mapHowmanyTrials(:);
    threshold = .3;
    plotThis = find(mapHowmanyTrials_linear >= threshold);
    t = text(X(plotThis),Y(plotThis),num2str(mapHowmanyTrials_linear(plotThis),'%.2f'),'fontsize',5/.45,'horizontalalignment','center','color','k');
    for i = 1:length(plotThis)
        t(i).String = t(i).String(2:end);
    end
    plotThis = find(mapHowmanyTrials_linear < threshold);
    t = text(X(plotThis),Y(plotThis),num2str(mapHowmanyTrials_linear(plotThis),'%.2f'),'fontsize',5/.45,'horizontalalignment','center','color','w');
    for i = 1:length(plotThis)
        t(i).String = t(i).String(2:end);
    end
    plotThis = find(isnan(mapHowmanyTrials_linear));
    t = text(X(plotThis),Y(plotThis),'0','fontsize',5/.45,'horizontalalignment','center','color','w');
    % Trial count
    nFill = [5 9 17 19]; nRep = [300 165 90 45];
    for i = 1:4
        subplot(21,20,342+(i-1)*20:357+(i-1)*20);
        imagesc([zeros(1,(19-nFill(i))/2) nRep(i)*ones(1,nFill(i)) zeros(1,(19-nFill(i))/2) ],[0 300]);
        ax = gca; colormap(gca,repmat(linspace(1,0,300)',1,3))
        ax.FontSize = 8/.45;
        ax.LineWidth = 1; ax.TickDir = 'in'; ax.TickLength = [0 .025];
        ax.YTick = []; ax.XTick = 1:3:19;
        ax.XTickLabel = {};
    end
    % Setting
    ax.XLabel.String = 'Relative direction of previous stimulus (\circ)';
    ax.XTickLabel = {'−90','−60','−30','0','30','60','90'};
    % Colorbar
    figure; subplot(21,20,[342:359 362:379 382:399 402:419])
    colormap(repmat(linspace(1,0,300)',1,3))
    cbar = colorbar;
    cbar.LineWidth = 1; cbar.Color = 'k'; cbar.FontSize = 8/.45;
    cbar.TickLabels = {'0','300    '};
    cbar.Label.String = 'Trial count';
    cbar.Label.VerticalAlignment = 'middle'; cbar.Label.Rotation = 270;
    cbar.Ticks = 0:1;cbar.TickLength = 0;

    %% :: Figure 2BC :: Conditional Bias
    figure;
    % Stimulus given response
    subplot(1,2,1); hold on
    plot([XLimCB(1) 0; XTickCB(2) 0],[0 YLimCB(1); 0 YTickCB(2)],'k--','linewidth',1/3); 
    blues = [156,198,99;86,177,218;42,111,246;60,79,156;42,45,101]/255;
    for i = 1:length(givenResp)
        shadedErrorBar(xBin,meanSBGR(i,:),ciSBGR(i,:), ...
            'lineprops',{'-','color',blues(i,:),'linewidth',2.5}, 'patchsaturation',.2);
    end
    % Setting
    ax = gca; pbaspect(ax, [1 1.1 1])
    ax.FontSize = 16;
    ax.XLabel.String = ['Relative direction of' newline 'previous stimulus (\circ)'];
    ax.YLabel.String = 'Response error (\circ)';
    ax.XLim = XLimCB; ax.YLim = YLimCB;
    ax.XTick = linspace(XTickCB(1),XTickCB(2),5); ax.YTick = linspace(YTickCB(1),YTickCB(2),5);
    ax.XTickLabel = {'−40','−20','0','20','40'}; ax.YTickLabel = {'−8','−4','    0','4','8'};
    offsetAxes(ax,XTickCB,YTickCB)

    % Response given stimulus
    subplot(1,2,2); hold on
    plot([XLimCB(1) 0; XTickCB(2) 0],[0 YLimCB(1); 0 YTickCB(2)],'k--','linewidth',1/3);
    reds = [242,186,69;232,144,62;224,101,55;218,79,52;178,50,47]/255;
    for i = 1:length(givenResp)
        shadedErrorBar(xBin,meanRBGS(i,:),ciRBGS(i,:), ...
            'lineprops',{'-','color',reds(i,:),'linewidth',2.5}, 'patchsaturation',.2);
    end
    % Setting
    ax = gca; pbaspect(ax, [1 1.1 1])
    ax.FontSize = 16;
    ax.XLabel.String = ['Relative direction of' newline 'previous response (\circ)'];
    ax.XLim = XLimCB; ax.YLim = YLimCB;
    ax.XTick = linspace(XTickCB(1),XTickCB(2),5); ax.YTick = linspace(YTickCB(1),YTickCB(2),5);
    ax.XTickLabel = {'−40','−20','0','20','40'}; ax.YTickLabel = {'−8','−4','    0','4','8'};
    offsetAxes(ax,XTickCB,YTickCB)

    %% :: Figure S5 :: Residual of Models
    figure; orient('landscape')
    set(gcf,'PaperUnits','centimeters','PaperSize',[40,21]);
    
    % Residual of the Response model
    subplot(1,21,1:4); hold on
    plot([XLimMB(1) 0; XTickMB(2) 0],[0 -25; 0 25],'k--','linewidth',2/3);
    scatter(dataCombined(:,6), resi_R_Combined,'markeredgecolor','none','markerfacealpha',.01,'markerfacecolor','k');
    plot(xBin,ResidualStimBiasCombined,'c','linewidth',2)
    % Setting
    ax = gca; pbaspect(ax, [1 1.1 1])
    ax.FontSize = 16;
    ax.XLabel.String = ['Relative direction of' newline 'previous stimulus (\circ)'];
    ax.YLabel.String = ['Residual of' newline 'the Response model (\circ)'];
    ax.XLim = XLimMB; ax.YLim = [-25 25];
    ax.XTick = linspace(XTickMB(1),XTickMB(2),nXTick); ax.YTick = -20:10:20;
    ax.XTickLabel = {'−90','','','0','','','90'}; ax.YTickLabel = {'  −20','−10','0','10','20'};
    offsetAxes(ax,-90,[-20 20])
    
    % Parameters
    load par_residual_S.mat
    subplot(1,21,7:9); hold on
    bar(1.5,par.modeGroup(1),'facecolor',[.8 .8 .8],'edgecolor','none','barwidth',.5,'showbaseline','off');
    plot([1 2],[0 0],'k-','linewidth',1)
    plot(normrnd(1.5,.1,32,1),par.modeIndv(:,1),'o','markeredgecolor',.5*ones(1,3),'linewidth',1,'markersize',8)
    errorbar(1.5,par.modeGroup(1),par.hdi95(1,1)-par.modeGroup(1),par.hdi95(1,2)-par.modeGroup(1),'k','linewidth',3,'capsize',0)
    [~,p] = ttest(par.modeIndv(:,1));
    mysigstar(gca, 1.5, -2.25, p);
    % Setting
    ax = gca; pbaspect([1,2,1]) 
    ax.FontSize = 16;
    ax.YLabel.String = 'Bias (\circ)';
    ax.TickLength = [.02 .025];
    ax.XLim = [.75 2.25]; ax.YLim = [-2 1];
    ax.XTick = []; ax.YTick = -2:0;
    ax.YTickLabel = {'−2','1','0'};
    offsetAxes(ax,[1 1],[-2 .5])
    
    % Residual from the Stimulus model
    subplot(1,21,13:16); hold on
    plot([XLimMB(1) 0; XTickMB(2) 0],[0 -25; 0 25],'k--','linewidth',2/3);
    scatter(dataCombined(:,7), resi_S_Combined,'markeredgecolor','none','markerfacealpha',.01,'markerfacecolor','k');
    plot(xMoving,ResidualRespBiasCombined,'y','linewidth',2)
    % Setting
    ax = gca; pbaspect(ax, [1 1.1 1])
    ax.FontSize = 16;
    ax.XLabel.String = ['Relative direction of' newline 'previous response (\circ)'];
    ax.YLabel.String = ['Residual of' newline 'the Stimulus model (\circ)'];
    ax.XLim = XLimMB; ax.YLim = [-25 25];
    ax.XTick = linspace(XTickMB(1),XTickMB(2),nXTick); ax.YTick = -20:10:20;
    ax.XTickLabel = {'−90','','','0','','','90'}; ax.YTickLabel = {'  −20','−10','0','10','20'};
    offsetAxes(ax,-90,[-20 20])
    
    % Parameters
    load par_residual_R.mat
    subplot(1,21,19:21); hold on
    bar(1.5,par.modeGroup(1),'facecolor',[.8 .8 .8],'edgecolor','none','barwidth',.5,'showbaseline','off');
    plot([1 2],[0 0],'k-','linewidth',1)
    plot(normrnd(1.5,.1,32,1),par.modeIndv(:,1),'o','markeredgecolor',.5*ones(1,3),'linewidth',1,'markersize',8)
    errorbar(1.5,par.modeGroup(1),par.hdi95(1,1)-par.modeGroup(1),par.hdi95(1,2)-par.modeGroup(1),'k','linewidth',3,'capsize',0)
    [~,p] = ttest(par.modeIndv(:,1));
    mysigstar(gca, 1.5, 8, p);
    % Setting
    ax = gca; pbaspect([1,2,1]) 
    ax.FontSize = 16;
    ax.YLabel.String = 'Bias (\circ)';
    ax.TickLength = [.02 .025];
    ax.XLim = [.75 2.25]; ax.YLim = [0 10];
    ax.XTick = []; ax.YTick = 0:2:8;
    ax.YTickLabel = {'   0','2','4','6','8'};
    offsetAxes(ax,[1 1],[0 8])
    
    %% :: Figure 3 :: DoG Prediction
    figure; orient(gcf,'landscape');
    set(gcf,'PaperUnits','centimeters','PaperSize',[34,21]);
    for mdl = 1:nModel
        subplot(2,nModel,mdl)
        clim = [-8 8];
        imagesc(DoGmap(:,:,mdl),clim); hold on;
        colormap(bwr);
        patch((repmat(X(isShort(:)),1,4)+[-.5 .5 .5 -.5])',(repmat(Y(isShort(:)),1,4)+[-.5 -.5 .5 .5])',[.9 .9 .9],'edgecolor',[.9 .9 .9])
        text([4 16],[4 16],'NA','horizontalalignment','center','fontsize',15)
        % Setting
        ax = gca; pbaspect(ax, [1 1 1]); ax.FontSize = 8/.45;
        ax.LineWidth = 1; ax.TickDir = 'in'; ax.TickLength = [.01 .025];
        ax.XTick = 1:3:19; ax.YTick = 1:3:19;
        ax.XTickLabel = {'−90','','','0','','','90'}; ax.YTickLabel = {};
        if mdl == 1
            ax.YTickLabel = {'90','','','0','','','   −90'};
            ax.YLabel.String = ['Relative direction of' newline 'previous response (\circ)'];
        elseif mdl == 2
            ax.XLabel.String = 'Relative direction of previous stimulus (\circ)';
        end
        title(mdlString{mdl},'fontweight','normal','fontsize',8/.45)
    end

    % DoG comparison
    subplot(2,10,12:20); hold on 
    plot([0 0],[.5 size(AIC,2)+.5],'k','linewidth',1)
    colors = [192 203 67; 114 179 72; 5 164 154; 10 196 216; 255 128 12]/255;
    [~,bestIdx] = min(mean(AIC));
    dAIC = AIC-AIC(:,bestIdx);
    dAICci = zeros(2,size(AIC,2));
    for i = 1:size(AIC,2)
        plot(dAIC(:,i),i*ones(size(dAIC,1),1),'o','markersize',8,'markeredgecolor','w','markerfacecolor',colors(i,:));
        dAICci(:,i) = bootci(10000,@mean,dAIC(:,i));
        errorbar(mean(dAIC(:,i)),i,dAICci(1,i)-mean(dAIC(:,i)),dAICci(2,i)-mean(dAIC(:,i)),'horizontal','k','linewidth',2,'capsize',13.5); % 13.5
        plot(mean(dAIC(:,i)),i,'ko','markersize',15,'linewidth',2,'markerfacecolor',colors(i,:));
        text(-52.5,i,mdlString{i},'horizontalalignment','right','fontsize',8/.45)
    end
    % Break in axis
    daxis = 400;
    fake = dAIC(1,1)-daxis;
    plot(fake,1,'o','markersize',8,'markeredgecolor','w','markerfacecolor',colors(1,:));
    % Setting
    ax = gca; pbaspect(ax,[3.5 1 1]);
    ax.FontSize = 8/.45;
    ax.XLim = [-50 700]; ax.YLim = [.5 size(AIC,2)+.5];
    ax.TickLength = [.01 .025];
    ax.XTick = 0:100:700; ax.YTick = 1:size(dAIC,2);
    ax.XTickLabel = {'0','100','200','300','400','500',num2str(600+daxis),num2str(700+daxis)};
    ax.YTickLabel = mdlString(1:end-1);
    ax.YColor = 'none';
    ax.XLabel.String = 'AIC difference';
    % Break in axis
    axes('Position',[.76 .1075 .015 .025]);
    py=[1 5];
    px1=[1 2];
    height=1;
    px2=px1+height;
    plot(px1,py,'k','LineWidth',2);hold all;
    plot(px2,py,'k','LineWidth',2);hold all;
    fill([px1 flip(px2)],[py flip(py)],'w','EdgeColor','none');
    box off;
    axis off;
    
    %% :: Figure 4 :: DoG Parameters -- Scatter connected by lines
    figure; load par_SR.mat
    red = [200 58 52]/255; blue = [42 45 116]/255;

    % Amplitude
    subplot(1,2,1); hold on;
    b = bar([par.modeGroup(1),par.modeGroup(3)]);
    b.FaceColor = [.8 .8 .8];
    b.EdgeColor = 'none';
    b.BarWidth = .5;
    b.ShowBaseLine = 'off';
    plot([.75 2.25],[0 0],'k-','linewidth',1)
    plot(repmat([1.125 1.875]',1,length(subjs)),par.modeIndv(:,[1 3])','color',[.5 .5 .5],'linewidth',2/3)
    plot(1.125*ones(length(subjs),1),par.modeIndv(:,1),'o','markersize',8,'markerfacecolor',blue,'markeredgecolor','w'); % size 80
    plot(1.875*ones(length(subjs),1),par.modeIndv(:,3),'o','markersize',8,'markerfacecolor',red,'markeredgecolor','w')
    errorbar(1:2,par.modeGroup([1,3]),par.modeGroup([1,3])-par.hdi95([1,3],1)',par.modeGroup([1,3])-par.hdi95([1,3],2)','k','linewidth',2,'linestyle','none','capsize',0)
    [~,p1] = ttest( par.modeIndv(:,1) );
    [~,p2] = ttest( par.modeIndv(:,3) );
    [~,p3] = ttest( par.modeIndv(:,6) );
    mysigstar(gca, 1, 4, p1);
    mysigstar(gca, 2, -4, p2);
    mysigstar(gca, [1 2], 20, p3, 'none', .15);
    % Setting
    ax = gca; pbaspect([3,5,1])
    ax.YLabel.String = 'Bias (\circ)';
    ax.XTickLabel = {'Stimulus','Response'};
    ax.TickLength = [.03 .025];
    ax.XLim = [.5 2.5]; ax.YLim = [-12 20];
    ax.XTick = 1:2; ax.YTick = -8:8:16;
    ax.YTickLabel = {'−8','0','8','16'};
    offsetAxes(ax,1,[ax.YTick(1) ax.YTick(end)])

    % Peak location
    subplot(1,2,2); hold on;
    b = bar([par.modeGroup(2),par.modeGroup(4)]);
    b.FaceColor = [.8 .8 .8];
    b.EdgeColor = 'none';
    b.BarWidth = .5;
    b.ShowBaseLine = 'off';
    plot([.75 2.25],[0 0],'k-','linewidth',1)
    [~,p] = ttest( par.modeIndv(:,7) );
    mysigstar(gca, [1 2], 80, p, 'none', .15);
    plot(repmat([1.125 1.875]',1,length(subjs)),par.modeIndv(:,[2 4])','color',[.5 .5 .5],'linewidth',2/3)
    plot(1.125*ones(length(subjs),1),par.modeIndv(:,2),'o','markersize',8,'markerfacecolor',blue,'markeredgecolor','w'); % size 80
    plot(1.875*ones(length(subjs),1),par.modeIndv(:,4),'o','markersize',8,'markerfacecolor',red,'markeredgecolor','w')    
    errorbar(1:2,par.modeGroup([2,4]),par.modeGroup([2,4])-par.hdi95([2,4],1)',par.modeGroup([2,4])-par.hdi95([2,4],2)','k','linewidth',2,'linestyle','none','capsize',0)
    % Setting
    ax = gca; pbaspect([3,5,1])
    ax.YLabel.String = 'Peak location (\circ)';
    ax.XTickLabel = {'Stimulus','Response'};
    ax.TickLength = [.03 .025];
    ax.XLim = [.5 2.5]; ax.YLim = [0 80];
    ax.XTick = 1:2; ax.YTick = 0:20:80;
    ax.YTickLabel = {'0','20','  40','60','80'};
    offsetAxes(ax,.75,[ax.YTick(1) ax.YTick(end)])
    
    %% :: Figure 5 :: Beta Coefficients
    figure;orient(gcf,'landscape');
    set(gcf,'PaperUnits','centimeters','PaperSize',[34,11]);
    load par_AR.mat
    red = [219 161 55]/255; blue = [110 177 227]/255; green = [0.4660 0.6740 0.1880];

    % Beta coefficient
    subplot(131); hold on;
    plot(repmat([.75 2.25],2,1)',[0 0; 1 1]','k-','linewidth',1)
    plot(2,par.modeGroup(1),'o','markersize',16,'markerfacecolor',blue+(1-blue)*(1-1/3),'markeredgecolor',blue,'linewidth',2)
    errorbar(2,par.modeGroup(1),par.hdi95(1,1)-par.modeGroup(1),par.hdi95(1,2)-par.modeGroup(1),...
        'color','k','linewidth',2,'capsize',0)
    plot(1,par.modeGroup(4),'o','markersize',16,'markerfacecolor',red+(1-red)*(1-1/3),'markeredgecolor',red,'linewidth',2)
    errorbar(1,par.modeGroup(4),par.hdi95(4,1)-par.modeGroup(4),par.hdi95(4,2)-par.modeGroup(4),...
        'color','k','linewidth',2,'capsize',0)
    % Setting
    ax = gca; pbaspect([1 1 1])
    ax.FontSize = 16;
    ax.YLabel.String = 'Beta coeffcients';
    ax.XLim = [.6 2.4]; ax.YLim = [-1.75 1.75];
    ax.XTick = 1:2; ax.YTick = -1.5:.5:1.5;
    ax.XTickLabel = {'Attraction','Repulsion'};
    ax.YTickLabel = {'','   −1','','0','','1',''};
    offsetAxes(ax,[.75 2.25],[-1.5 1.5])
    
    %
    subplot(132); hold on
    plot(par.modeIndv(:,5),par.modeIndv(:,4),'o','markeredgecolor','k','markerfacecolor',red,'markersize',10,'linewidth',1)
    plot(par.modeIndv(1,5)-3,par.modeIndv(1,4),'o','markeredgecolor','k','markerfacecolor',red,'markersize',10,'linewidth',1)
    % Setting
    ax = gca; pbaspect([1 1 1])
    ax.FontSize = 16;
    ax.XLabel.String = 'Attraction bias (\circ)';        
    ax.YLabel.String = ['Beta coefficients' newline 'for attraction curve'];
    ax.XLim = [1 13]; ax.YLim = [.65 1.35];
    ax.XTick = 2:2:12; ax.YTick = 0.7:.1:1.3;
    ax.XTickLabel = {'2','','6','','10','15'};
    ax.YTickLabel = {'','   0.8','','1','','1.2',''};
    offsetAxes(ax,2,[.7 1.3])
    
    %
    subplot(133); hold on
    plot(par.modeIndv(:,2),par.modeIndv(:,1),'o','markeredgecolor','k','markerfacecolor',blue,'markersize',10,'linewidth',1)
    mdl = fitlm(par.modeIndv(:,2),par.modeIndv(:,1));
    xx = [-6 0];
    plot(xx,mdl.Coefficients.Estimate(1)+mdl.Coefficients.Estimate(2)*xx,'k','linewidth',2)
    % Setting
    ax = gca; pbaspect([1 1 1])
    ax.FontSize = 16;
    ax.XLabel.String = 'Repulsion bias (\circ)';        
    ax.YLabel.String = ['Beta coefficients' newline 'for repulsion curve'];
    ax.XLim = [-6.5 0.5]; ax.YLim = [-1.25 -.55];
    ax.XTick = -6:2:0; ax.YTick = -1.2:.2:-.6;
    ax.XTickLabel = {'−6','−4','−2','0'};
    ax.YTickLabel = {'  −1.2','−1','−0.8','−0.6'};
    offsetAxes(ax,-6,[-1.2 -.6])
    
    % Break in axis
    axes('Position',[.585 .23 .015 .025]);
    py=[1 5];
    px1=[1 2];
    height=1;
    px2=px1+height;
    plot(px1,py,'k','LineWidth',2);hold all;
    plot(px2,py,'k','LineWidth',2);hold all;
    fill([px1 flip(px2)],[py flip(py)],'w','EdgeColor','none');
    box off;
    axis off;
    
    %% :: Figure S6 :: Group Difference (boxplot)
    figure; orient(gcf,'landscape');
    set(gcf,'PaperUnits','centimeters','PaperSize',[36,24]);

    % Stimulus model -- Bias magnitude
    load par_S.mat
    subplot(2,3,1); hold on;
    parMat = reshape(par.modeIndv(:,1),4,8)';
    plot(median(parMat),'k-','linewidth',1)
    h = boxplot(parMat,'boxstyle','filled','colors','k','medianstyle','target','plotstyle','compact','width',.25,'notch','off');
    h = findobj(h);
    box off
    for i = 1:4
        h(1+(i-1)*5).LineWidth = 1;
        h(3+(i-1)*5).MarkerSize = 15;
        h(3+(i-1)*5).LineWidth = 1;
        h(4+(i-1)*5).MarkerSize = 30;
        h(5+(i-1)*5).MarkerFaceColor = 'k';
        h(5+(i-1)*5).XData = i;
    end
    % Setting
    ax = gca; pbaspect([1.5,1,1])
    ax.FontSize = 16;
    ax.XTickLabel = {'20','40','80','180'};
    ax.YTickLabel={'','0','','4','','     8',''};
    ax.XLim = [.5 4.5];
    ax.YLim = [-2.5 10.5];
    ax.XTick = 1:4; ax.YTick = -2:2:10;
    offsetAxes(ax,1,-2)

    % Stimulus model -- Peak location
    subplot(2,3,4); hold on; 
    parMat = reshape(par.modeIndv(:,2),4,8)';
    plot(median(parMat),'k-','linewidth',1)
    h = boxplot(parMat,'boxstyle','filled','colors','k','medianstyle','target','plotstyle','compact','width',.25,'notch','off','symbol','');
    h = findobj(h);
    box off
    for i = 1:4
        h(1+(i-1)*5).LineWidth = 1;
        h(3+(i-1)*5).MarkerSize = 15;
        h(3+(i-1)*5).LineWidth = 1;
        h(4+(i-1)*5).MarkerSize = 30;
        h(5+(i-1)*5).MarkerFaceColor = 'k';
        h(5+(i-1)*5).XData = i;
    end
    % Setting
    ax = gca; pbaspect([1.5,1,1])
    ax.FontSize = 16;
    ax.XTickLabel = {'20','40','80','180'};
    ax.YTickLabel={'16','20','24','    28'};
    ax.XLim = [.5 4.5];
    ax.YLim = [15 29];
    ax.XTick = 1:4; ax.YTick = 16:4:28;
    offsetAxes(ax,1,16)

    % SR model Repulsion -- Bias magnitude
    load par_SR.mat
    subplot(2,3,2); hold on; 
    red = [200 58 52]/255; blue = [0, 0.4470, 0.7410];
    parMat = reshape(par.modeIndv(:,1),4,8)';
    plot(median(parMat),'color',blue,'linewidth',1)
    h = boxplot(parMat,'boxstyle','filled','colors',blue,'medianstyle','target','plotstyle','compact','width',.25,'notch','off','symbol','');
    h = findobj(h);
    box off
    for i = 1:4
        h(1+(i-1)*5).LineWidth = 1;
        h(3+(i-1)*5).MarkerSize = 15;
        h(3+(i-1)*5).LineWidth = 1;
        h(4+(i-1)*5).MarkerSize = 30;
        h(4+(i-1)*5).MarkerEdgeColor = blue;
        h(5+(i-1)*5).MarkerFaceColor = blue;
        h(5+(i-1)*5).XData = i;
    end
    % Setting
    ax = gca; pbaspect([1.5,1,1])
    ax.FontSize = 16;
    ax.XTickLabel = {'20','40','80','180'};
    ax.YTickLabel = {'−10','−8','−6','−4','−2'};
    ax.XLim = [.5 4.5]; 
    ax.YLim = [-10.5 -1.5];
    ax.XTick = 1:4; ax.YTick = -10:2:-2;
    offsetAxes(ax,1,-10);

    % SR model Repulsion -- Peak location
    subplot(2,3,5); hold on; 
    parMat = reshape(par.modeIndv(:,2),4,8)';
    plot(median(parMat),'color',blue,'linewidth',1)
    h = boxplot(parMat,'boxstyle','filled','colors',blue,'medianstyle','target','plotstyle','compact','width',.25,'notch','off','symbol','');
    h = findobj(h);
    box off
    for i = 1:4
        h(1+(i-1)*5).LineWidth = 1;
        h(3+(i-1)*5).MarkerSize = 15;
        h(3+(i-1)*5).LineWidth = 1;
        h(4+(i-1)*5).MarkerSize = 30;
        h(4+(i-1)*5).MarkerEdgeColor = blue;
        h(5+(i-1)*5).MarkerFaceColor = blue;
        h(5+(i-1)*5).XData = i;
    end
    % Setting
    ax = gca; pbaspect([1.5,1,1])
    ax.FontSize = 16;
%     ax.XLabel.String = 'Propagation half width (\circ)';
    ax.XTickLabel = {'20','40','80','180'};
    ax.YTickLabel = {'25','30','35','40','   45'};
    ax.XLim = [.5 4.5];
    ax.YLim = [23 47];
    ax.XTick = 1:4; ax.YTick = 25:5:45;
    offsetAxes(ax,1,25);

    % SR model Attraction -- Bias magnitude
    subplot(2,3,3); hold on;
    parMat = reshape(par.modeIndv(:,3),4,8)';
    plot(median(parMat),'color',red,'linewidth',1)
    h = boxplot(parMat,'boxstyle','filled','colors',red,'medianstyle','target','plotstyle','compact','width',.25,'notch','off','symbol','');
    h = findobj(h);
    box off
    for i = 1:4
        h(1+(i-1)*5).LineWidth = 1;
        h(3+(i-1)*5).MarkerSize = 15;
        h(3+(i-1)*5).LineWidth = 1;
        h(4+(i-1)*5).MarkerSize = 30;
        h(4+(i-1)*5).MarkerEdgeColor = red;
        h(5+(i-1)*5).MarkerFaceColor = red;
        h(5+(i-1)*5).XData = i;
    end
    % Setting
    ax = gca; pbaspect([1.5,1,1])
    ax.FontSize = 16;
    ax.XTickLabel = {'20','40','80','180'};
    ax.XLim = [.5 4.5];
    ax.YLim = [3 13];
    ax.XTick = 1:4; ax.YTick = 4:2:12;
    offsetAxes(ax,1,4);

    % SR model Attraction -- Peak location
    subplot(2,3,6); hold on;
    parMat = reshape(par.modeIndv(:,4),4,8)';
    plot(median(parMat),'color',red,'linewidth',1)
    h = boxplot(parMat,'boxstyle','filled','colors',red,'medianstyle','target','plotstyle','compact','width',.25,'notch','off','symbol','');
    h = findobj(h);
    box off
    for i = 1:4
        h(1+(i-1)*5).LineWidth = 1;
        h(3+(i-1)*5).MarkerSize = 15;
        h(3+(i-1)*5).LineWidth = 1;
        h(4+(i-1)*5).MarkerSize = 30;
        h(4+(i-1)*5).MarkerEdgeColor = red;
        h(5+(i-1)*5).MarkerFaceColor = red;
        h(5+(i-1)*5).XData = i;
    end
    % Setting
    ax = gca; pbaspect([1.5,1,1])
    ax.FontSize = 16;
    ax.XTickLabel = {'20','40','80','180'};
    ax.XLim = [.5 4.5];
    ax.YLim = [18 52];
    ax.XTick = 1:4; ax.YTick = 20:10:50;
    offsetAxes(ax,1,20);
    
end

%% :: FUNCTION ::

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

% omitnan
% alpha = alpha(~isnan(alpha));
  
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
r = sum(w.*exp(1i*alpha),dim,'omitnan');

% obtain length 
r = abs(r)./sum(w,dim);

% for data with known spacing, apply correction factor to correct for bias
% in the estimation of r (see Zar, p. 601, equ. 26.16)
if d ~= 0
  c = d/2/sin(d/2);
  r = c*r;
end
end

%% :: FUNCTION :: circ_corrcc
function [rho, pval] = circ_corrcc(alpha1, alpha2)

if size(alpha1,2) > size(alpha1,1)
	alpha1 = alpha1';
end

if size(alpha2,2) > size(alpha2,1)
	alpha2 = alpha2';
end

if length(alpha1)~=length(alpha2)
  error('Input dimensions do not match.')
end

% compute mean directions
n = length(alpha1);
alpha1_bar = circ_mean(alpha1);
alpha2_bar = circ_mean(alpha2);

% compute correlation coeffcient from p. 176
num = sum(sin(alpha1 - alpha1_bar) .* sin(alpha2 - alpha2_bar));
den = sqrt(sum(sin(alpha1 - alpha1_bar).^2) .* sum(sin(alpha2 - alpha2_bar).^2));
rho = num / den;	

% compute pvalue
l20 = mean(sin(alpha1 - alpha1_bar).^2);
l02 = mean(sin(alpha2 - alpha2_bar).^2);
l22 = mean((sin(alpha1 - alpha1_bar).^2) .* (sin(alpha2 - alpha2_bar).^2));

ts = sqrt((n * l20 * l02)/l22) * rho;
pval = 2 * (1 - normcdf(abs(ts)));
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

%% :: FUNCTION :: offsetAxes
function offsetAxes( ax, x, y )
if ~exist('ax', 'var'), ax = gca; end
if ~exist('x', 'var')
    x = [ax.XLim(1)+(ax.XTick(2)-ax.XTick(1))/4 ax.XLim(2)-(ax.XTick(2)-ax.XTick(1))/4];
elseif numel(x)<2
    x = [x ax.XLim(2)-(x-ax.XLim(1))];
end
if ~exist('y', 'var')
    y = [ax.YLim(1)+(ax.YTick(2)-ax.YTick(1))/4 ax.YLim(2)-(ax.YTick(2)-ax.YTicm(1))/4];
elseif numel(y)<2
    y = [y ax.YLim(2)-(y-ax.YLim(1))];
end
% remove the x and y ticks if they are outside the visible X and YLim
ax.XTick = ax.XTick( ax.XTick >= x(1) & ax.XTick <= x(2) );
ax.YTick = ax.YTick( ax.YTick >= y(1) & ax.YTick <= y(2) );
% position the X and Y labels on the center of the visible axes
% ax.XLabel.Position(1) = mean(x);
% ax.YLabel.Position(2) = mean(y);
% this will keep the changes constant even when resizing axes
addlistener (ax, 'MarkedClean', @(obj,event)resetVertex(ax,x,y));
end
function resetVertex( ax, x, y )
% extract the x axis vertext data
% X, Y and Z row of the start and end of the individual axle.
ax.XRuler.Axle.VertexData(1,1) = x(1);
ax.XRuler.Axle.VertexData(1,2) = x(2);
% repeat for Y (set 2nd row)
ax.YRuler.Axle.VertexData(2,1) = y(1);
ax.YRuler.Axle.VertexData(2,2) = y(2);
end

%% :: FUNCTION :: mysigstar
function h = mysigstar(ax, xpos, ypos, pval, whichWay, howMuch)
% replaces sigstar, which doesnt work anymore in matlab 2014b

if ~exist('ax', 'var'); ax = gca; end
if ~exist('whichWay', 'var'), whichWay = 'none'; end
if ~exist('howMuch', 'var'), howMuch = .05*range(get(gca,'ylim')); end

if numel(ypos) > 1
    assert(ypos(1) == ypos(2), 'line wont be straight!');
    ypos = ypos(1);
end

% draw line
hold on;
if numel(xpos) > 1
    % plot the horizontal line
    p = plot(ax, [xpos(1), xpos(2)], ...
        [ypos ypos], '-', 'LineWidth', 1, 'color', 'k');
    
    % also add small downward ticks
    switch whichWay
        case 'down'
            plot([xpos(1) xpos(1)], [ypos ypos-howMuch], '-k', 'LineWidth', 1.5);
            plot([xpos(2) xpos(2)], [ypos ypos-howMuch], '-k', 'LineWidth', 1.5);
        case 'up'
            plot([xpos(1) xpos(1)], [ypos ypos+howMuch], '-k', 'LineWidth', 1.5);
            plot([xpos(2) xpos(2)], [ypos ypos+howMuch], '-k', 'LineWidth', 1.5);
        case 'none'
    end
    txtBg = 'w';
else
    txtBg = 'none';
end

fz = 30; fontweight = 'normal';
ypos = ypos-0.025*range(get(gca,'ylim'));
if pval < 1e-3
    txt = '***';
elseif pval < 1e-2
    txt = '**';
elseif pval < 0.05
    txt = '*';
elseif ~isnan(pval)
    % this should be smaller
    txt = 'n.s.';
    %txt = '';
%     fz = fz * 2/3;
    fz = 15;
    ypos = ypos+0.025*range(get(gca,'ylim'));
    fontweight = 'normal';
else
    return
end

% draw the stars in the bar
h = text(mean(xpos), mean(ypos), txt, ...
    'horizontalalignment', 'center', 'backgroundcolor', ...
    txtBg, 'margin', 6, 'fontsize', fz, 'fontweight', fontweight, 'Parent', ax);
end

