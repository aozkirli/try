clearvars;
% close all;
rng default
plot_or_not = true;
DoG = @(x,a,w) x.*a.*w.*sqrt(2)./exp(-.5).*exp(-1.*(w.*x).^2);

%% Load Data
cd ..
rawdataIndv = cell(32,2);
dataIndv = cell(32,2);
dataCombined = [];
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
    currdata(:,4) = bound_180( currdata(:,4) - c_mean(currdata(:,4)) );
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
cd([pwd,'/woCorrection'])


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
%% Subject Loop
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
                load par_woCorrection_S.mat
                predictedBias = DoG( data(:,6), par.modeIndv(s,1), 1/(sqrt(2)*par.modeIndv(s,2)) );
                noiseSigma = par.modeIndv(s,3);
            case 2
                load par_woCorrection_R.mat
                predictedBias = DoG( data(:,7), par.modeIndv(s,1), 1/(sqrt(2)*par.modeIndv(s,2)) );
                noiseSigma = par.modeIndv(s,3);
            case 3
                load par_woCorrection_SR.mat
                predictedBias = DoG( data(:,6), par.modeIndv(s,1), 1/(sqrt(2)*par.modeIndv(s,2)) ) + DoG( data(:,7), par.modeIndv(s,3), 1/(sqrt(2)*par.modeIndv(s,4)) );
                noiseSigma = par.modeIndv(s,5);
        end
        AIC(s,mdl) = 2*length(par.modeGroup)-2*sum(log( circ_vmpdf( deg2rad(data(:,4)), deg2rad(predictedBias), 1/deg2rad(noiseSigma)^2 ) ));
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
                    mdlString{mdl} = ['Stimulus &' newline 'Response'];
                    tempmap = DoG( xBin(X), par.modeGroup(1), 1/(sqrt(2)*par.modeGroup(2)) ) + DoG( fliplr(xBin(Y)), par.modeGroup(3), 1/(sqrt(2)*par.modeGroup(4)) );
            end
            DoGmap(:,:,mdl) = reshape(tempmap,length(xBin),length(xBin));
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
    XLimCB = [-45 45]; YLimCB = [-12.5 12.5]; XTickCB = [-40 40]; YTickCB = [-10 10];
    XLimMB = [-100 100]; YLimMB = YLimCB; XTickMB = [-90 90]; YTickMB = YTickCB; nXTick = 7;
    
    %% :: Figure S6AB :: Marginal Bias
    figure; tiledlayout(1,2,'tilespacing','compact');
    % Marginal stimulus bias
    nexttile; hold on
    plot([XLimMB(1) 0; XTickMB(2) 0],[0 YLimMB(1); 0 YTickMB(2)],'k--','linewidth',2/3);
    meanMSB = c_mean(MarginalStimBias);
    ciSigma = tinv(0.975,sum(~isnan(MarginalStimBias))-1).*c_std(MarginalStimBias,~isnan(MarginalStimBias))./sqrt(sum(~isnan(MarginalStimBias)));
    shadedErrorBar(xBin,meanMSB,ciSigma,'lineprops',{'-','color','k','linewidth',2});
    % Setting
    ax = gca; pbaspect(ax, [1 1.1 1])
    ax.FontSize = 16;
    ax.XLabel.String = ['Relative direction of' newline 'previous stimulus (\circ)'];
    ax.YLabel.String = 'Response error (\circ)';
    ax.XLim = XLimMB; ax.YLim = YLimMB;
    ax.XTick = linspace(XTickMB(1),XTickMB(2),nXTick); ax.YTick = linspace(YTickMB(1),YTickMB(2),5);
    ax.XTickLabel = {'−90','','','0','','','90'}; ax.YTickLabel = {'−10','−5','    0','5','10'};
    offsetAxes(ax,XTickMB,YTickMB)

    % Marginal response bias
    nexttile; hold on
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
    ax.XTickLabel = {'−90','','','0','','','90'}; ax.YTickLabel = {'−10','−5','    0','5','10'};
    offsetAxes(ax,XTickMB,YTickMB)

    %% :: Figure S6CD :: Conditional Bias
    figure; tiledlayout(1,2,'tilespacing','compact');
    % Stimulus given response
    nexttile; hold on
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
    ax.XLim = XLimCB; ax.YLim = YLimCB;
    ax.XTick = linspace(XTickCB(1),XTickCB(2),5); ax.YTick = linspace(YTickCB(1),YTickCB(2),5);
    ax.XTickLabel = {'−40','−20','0','20','40'}; ax.YTickLabel = {'−10','−5','    0','5','10'};
    offsetAxes(ax,XTickCB,YTickCB)

    % Response given stimulus
    nexttile; hold on
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
    ax.XTickLabel = {'−40','−20','0','20','40'}; ax.YTickLabel = {'−10','−5','    0','5','10'};
    offsetAxes(ax,XTickCB,YTickCB)

    %% :: Figure S6EFGHI :: DoG Prediction
    figure; orient(gcf,'landscape');
    set(gcf,'PaperUnits','centimeters','PaperSize',[40,24]);
    for mdl = 1:nModel+1
        subplot(2,nModel+1,mod(mdl,4)+1)
        clim = [-10 10];
        imagesc(DoGmap(:,:,mdl),clim); hold on;
        colormap(bwr);
        patch((repmat(X(isShort(:)),1,4)+[-.5 .5 .5 -.5])',(repmat(Y(isShort(:)),1,4)+[-.5 -.5 .5 .5])',[.9 .9 .9],'edgecolor',[.9 .9 .9])
        text([4 16],[4 16],'NA','horizontalalignment','center','fontsize',15)
        % Setting
        ax = gca; pbaspect(ax, [1 1 1]); ax.FontSize = 16;
        ax.LineWidth = 1; ax.TickDir = 'in'; ax.TickLength = [.01 .025];
        ax.XTick = 1:3:19; ax.YTick = 1:3:19;
        ax.XTickLabel = {'−90','','','0','','','90'}; ax.YTickLabel = {};
        if mdl == 4
            ax.YTickLabel = {'90','','','0','','','   −90'};
            ax.YLabel.String = ['Relative direction of' newline 'previous response (\circ)'];
        elseif mdl == 2
        text(-20,23.5,'Relative direction of previous stimulus (\circ)','fontsize',16);
        end
        switch mdl
            case {1,2,4}
                title(mdlString{mdl},'fontweight','normal','fontsize',16)
            case 3
                title('Stimulus & Response','fontweight','normal','fontsize',16)
        end
    end

    % DoG comparison
    subplot(2,10,11:14); hold on 
    plot([0 0],[0 size(AIC,2)+.5],'k','linewidth',1)
    colors = [192 203 67; 114 179 72; 5 164 154; 10 196 216; 255 128 12]/255;
    [~,bestIdx] = min(mean(AIC));
    dAIC = AIC-AIC(:,bestIdx);
    dAICci = zeros(2,size(AIC,2));
    for i = 1:size(AIC,2)
        plot(dAIC(:,i),i*ones(size(dAIC,1),1),'o','markersize',8,'markeredgecolor','w','markerfacecolor',colors(i,:));
        dAICci(:,i) = bootci(10000,@mean,dAIC(:,i));
        errorbar(mean(dAIC(:,i)),i,dAICci(1,i)-mean(dAIC(:,i)),dAICci(2,i)-mean(dAIC(:,i)),'horizontal','k','linewidth',2,'capsize',13.5); % 13.5
        plot(mean(dAIC(:,i)),i,'ko','markersize',15,'linewidth',2,'markerfacecolor',colors(i,:));
        text(-52.5,i,mdlString{i},'horizontalalignment','right','fontsize',16)
    end
    % Break in axis
    daxis = 600;
    fake = dAIC(1,1)-daxis;
    plot(fake,1,'o','markersize',8,'markeredgecolor','w','markerfacecolor',colors(1,:));

    % Setting
    ax = gca; pbaspect(ax,[1.25 1 1]);
    ax.FontSize = 16;
    ax.XLim = [-50 800]; ax.YLim = [.5 size(AIC,2)+.5];
    ax.TickLength = [.01 .025];
    ax.XTick = 0:100:800; ax.YTick = 1:size(dAIC,2);
    ax.XTickLabel = {'0','100','200','300','400','500','600',num2str(600+daxis),num2str(700+daxis)};
    ax.YTickLabel = mdlString(1:end-1);
    ax.YColor = 'none';
    ax.XLabel.String = 'AIC difference';
    % Break in axis
    axes('Position',[.37 .1075 .015 .025]);
    py=[1 5];
    px1=[1 2];
    height=1;
    px2=px1+height;
    plot(px1,py,'k','LineWidth',2);hold all;
    plot(px2,py,'k','LineWidth',2);hold all;
    fill([px1 flip(px2)],[py flip(py)],'w','EdgeColor','none');
    box off;
    axis off;

    %% :: Figure S6JK :: DoG Parameters -- Scatter connected by lines
    figure; 
    load par_woCorrection_SR.mat
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
    mysigstar(gca, 1, 5, p1);
    mysigstar(gca, 2, -5, p2);
    mysigstar(gca, [1 2], 25, p3, 'none', .15);
    % Setting
    ax = gca; pbaspect([3,4.25,1])
    ax.YLabel.String = 'Bias (\circ)';
    ax.XTickLabel = {'Stimulus','Response'};
    ax.TickLength = [.03 .025];
    ax.XLim = [.5 2.5]; ax.YLim = [-15 25];
    ax.XTick = 1:2; ax.YTick = -10:10:20;
    ax.YTickLabel = {'−10','0','10','20'};
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
    ax = gca; pbaspect([3,4.25,1])
    ax.YLabel.String = 'Peak location (\circ)';
    ax.XTickLabel = {'Stimulus','Response'};
    ax.TickLength = [.03 .025];
    ax.XLim = [.5 2.5]; ax.YLim = [0 80];
    ax.XTick = 1:2; ax.YTick = 0:20:80;
    ax.YTickLabel = {'0','20','  40','60','80'};
    offsetAxes(ax,.75,[ax.YTick(1) ax.YTick(end)])
    
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

