addpath(genpath('C:\Users\chare\Google Drive\Work\09_Code\BEIM_toolbox\matlabv\'))
addpath('C:\Users\chare\Google Drive\Work\11_WorkInProgress\EPFL_Ambizione\Scripts\modeling\SD_models')
addpath('C:\Users\chare\Google Drive\Work\11_WorkInProgress\EPFL_Ambizione\Experiments\Functions\sdp_tool')
addpath('/Users/ayberkozkirli/Library/CloudStorage/GoogleDrive-ayberk.ozkirli@epfl.ch/.shortcut-targets-by-id/14P9enoQnGTKLDkkJnld8K5Wh-X7GmXcB/EPFL_Ambizione/Scripts/modeling/SD_models')
addpath('/Users/ayberkozkirli/Library/CloudStorage/GoogleDrive-ayberk.ozkirli@epfl.ch/.shortcut-targets-by-id/14P9enoQnGTKLDkkJnld8K5Wh-X7GmXcB/EPFL_Ambizione/Experiments/Functions/sdp_tool')
addpath(genpath('/Users/ayberkozkirli/Library/CloudStorage/GoogleDrive-ayberk.ozkirli@epfl.ch/My Drive/EPFL-PhD/Projects/SerialDependence/BEIM_toolbox'))
free
figure
del = -90:90;
sds = 8:3:20;
colors      = cbrewer('seq','PuBu',numel(sds)+10);
colors      = colors(11:end,:);
for k = 1:numel(sds)
    theta       = datasample(0:179,1000)';
    Se          = sds(k)/2;
    Sd          = 10;
    Pt          = 0.5;
    tau         = 1.2;
    [sdp,model] = SD_model_BayesianDecoding(Se,Sd,Pt,tau,theta)

    [delta,idx] = sort(model.delta);

    g = grpstats(model.error(idx),delta,'mean');
    g = sdp_mvav(nanunique(delta),g,21,[],[]);

    subplot(223)
    plot(g.x,g.m,'-','linewidth',2,'color',colors(k,:));hold on
    vline(0,'k');hline(0,'k');
    ylim([-10 10])

    g = grpstats(model.sigma(idx),delta,'mean');
    g = sdp_mvav(nanunique(delta),g,21,[],[]);

    subplot(224)
    plot(g.x,g.m,'-','linewidth',2,'color',colors(k,:));hold on
    h=hline(Se*2,'k--');
    h.Color = [0 0 0 .4];
    ylim([5 23])
    drawnow

    sigma = sds(k);
    wp = 1./(2+(del/sigma).^2);
    bias = wp.*del;
    sigma_dec = sqrt(wp.^2.*sigma.^2+(1-wp).^2.*sigma^2);

    subplot(221)
    plot(del,bias,'-','linewidth',2,'color',colors(k,:));hold on
    vline(0,'k');hline(0,'k');
    ylim([-10 10])

    subplot(222)
    plot(del,sigma_dec,'-','linewidth',2,'color',colors(k,:));hold on
    h=hline(Se*2,'k--');
    h.Color = [0 0 0 .4];
    ylim([5 23])
    drawnow
end
subplot(221)
xlabel('\Delta (°)');
ylabel('Bias (°)');
set(gca,'FontSize',15,'FontName','Times New Roman')

subplot(222)
xlabel('\Delta (°)');
ylabel('Error Scatter (°)');
set(gca,'FontSize',15,'FontName','Times New Roman')
annotation('textbox', [0.1, 0.95, 0.8, 0.05], 'String', 'SIMPLE CUE INTEGRATION', ...
    'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold', ...
    'LineStyle', 'none');

subplot(223)
xlabel('\Delta (°)');
ylabel('Bias (°)');
set(gca,'FontSize',15,'FontName','Times New Roman')

subplot(224)
xlabel('\Delta (°)');
ylabel('Error Scatter (°)');
set(gca,'FontSize',15,'FontName','Times New Roman')
annotation('textbox', [0.1, 0.48, 0.8, 0.05], 'String', 'BAYESIAN INTEGRATION', ...
    'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold', ...
    'LineStyle', 'none');
% 
% subplot(131)
% plot(delta,bias,'LineWidth',2,'Color',colors(ss,:));
% hold on
% subplot(132)
% h(ss)=plot(delta,(sigma_dec),'LineWidth',2,'Color',colors(ss,:));
% hold on
% m = sdp_mvav(model.delta,model.error,[],2,180)
% figure
% plot(m.x,m.std)