codeobsid = unique(tmp.codeobsid);
strong_bias_subj_iqr = nan(181,length(codeobsid));
strong_scatt_subj_iqr = nan(181,length(codeobsid));
for co = 1:length(codeobsid)
    
    tt = strong(strcmp(strong.codeobsid,codeobsid(co)),:);  
    
    biass = nan(1,181);
    scatt = nan(1,181);
    biass(ismember([-90:90],unique(tt.delta))) = grpstats(tt.error_iqr,tt.delta,'mean');
    scatt(ismember([-90:90],unique(tt.delta))) = grpstats(tt.error_iqr,tt.delta,'std');scatt(scatt==0)=NaN;
    strong_bias_subj_iqr(:,co)   = movmean(biass,30,'omitnan')-mean(movmean(biass,30,'omitnan'));
    strong_scatt_subj_iqr(:,co)  = movmean(scatt,30,'omitnan')-mean(movmean(scatt,30,'omitnan'));
end
codeobsid = unique(weak.codeobsid);
weak_bias_subj_iqr = nan(181,length(codeobsid));
weak_scatt_subj_iqr = nan(181,length(codeobsid));
for co = 1:length(codeobsid)
    tt = weak(strcmp(weak.codeobsid,codeobsid(co)),:);
    
    biass = nan(1,181);
    scatt = nan(1,181);
    biass(ismember([-90:90],unique(tt.delta))) = grpstats(tt.error_iqr,tt.delta,'mean');
    scatt(ismember([-90:90],unique(tt.delta))) = grpstats(tt.error_iqr,tt.delta,'std');scatt(scatt==0)=NaN;
    weak_bias_subj_iqr(:,co)   = movmean(biass,30,'omitnan')-mean(movmean(biass,30,'omitnan'));
    weak_scatt_subj_iqr(:,co)  = movmean(scatt,30,'omitnan')-mean(movmean(scatt,30,'omitnan'));
end
subplot(142)
[p2,f2] = plot_mean_ci(repmat(-90:90,1,length(codeobsid)),strong_bias_subj_iqr(:),1,1,'#FF7F0E','mean',0);
[p3,f3] = plot_mean_ci(repmat(-90:90,1,length(codeobsid)),weak_bias_subj_iqr(:),1,1,'#17BECF','mean',0);
% set(f2,'Vertices',[f2.Vertices(:,1), f2.Vertices(:,2)- p2.YData(1)])
% set(p2, 'YData', p2.YData-p2.YData(1));
% set(f3,'Vertices',[f3.Vertices(:,1), f3.Vertices(:,2)- p3.YData(1)])
% set(p3, 'YData', p3.YData-p3.YData(1));
subplot(143)
[p2,f2] = plot_mean_ci(repmat(-90:90,1,length(codeobsid)),strong_scatt_subj_iqr(:),1,1,'#FF7F0E','mean',0);
[p3,f3] = plot_mean_ci(repmat(-90:90,1,length(codeobsid)),weak_scatt_subj_iqr(:),1,1,'#17BECF','mean',0);
% set(f2,'Vertices',[f2.Vertices(:,1), f2.Vertices(:,2)- p2.YData(1)])
% set(p2, 'YData', p2.YData-p2.YData(1));
% set(f3,'Vertices',[f3.Vertices(:,1), f3.Vertices(:,2)- p3.YData(1)])
% set(p3, 'YData', p3.YData-p3.YData(1));