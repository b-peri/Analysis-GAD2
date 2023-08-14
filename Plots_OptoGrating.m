%% Figure 1: Plots OptoGratings

%% Colors

cols.off = [0 0 0];
cols.on_GAD2 = [0 0.4470 0.7410];
cols.on_Cont = [0.9290 0.6940 0.1250];
cols.xline = [143 143 143]/255;
cols.refline = [0.35 0.35 0.35];
cols.err_bar = [26 35 126]/255;
cols.err_bar_Cont = [199, 143, 0]/255;
cols.opto = [0.3010 0.7450 0.9330];
cols.Act = [0.8500 0.3250 0.0980];
cols.Inh = [0.9290 0.6940 0.1250];

%% Leonie's Figures pt. 2: Summary Stats

figure;

% --- GAD2 Mice ---

% Normalized Response Summary All Cells - Normalized to Peak
subplot(2,3,1); hold on;
shadedErrorBar(DataOut_OG_GAD2.Overall.PSTHBinCenters, DataOut_OG_GAD2.Overall.PSTHMean_Off_norm, DataOut_OG_GAD2.Overall.PSTHSEM_Off_norm, 'lineProps', {'Color', cols.off});
shadedErrorBar(DataOut_OG_GAD2.Overall.PSTHBinCenters, DataOut_OG_GAD2.Overall.PSTHMean_On_norm, DataOut_OG_GAD2.Overall.PSTHSEM_On_norm, 'lineProps', {'Color', cols.on_GAD2});
xline(0,'--','Color',cols.xline); xline(1,'--','Color',cols.xline);
% xline(-0.03,'g-'); xline(0.98, 'g-');
ylabel('Normalized Response');
xlabel('Time from stim. onset (s)');
xlim([min(DataOut_OG_GAD2.Overall.PSTHtime) max(DataOut_OG_GAD2.Overall.PSTHtime)]);
% ylim([-0.1 1]);
yline(0, 'Color', cols.xline);
% text(0.6, 2, sprintf('%g units', DataOut_OG_GAD2.Overall.NCells), 'FontSize',15);
legend('Laser Off', 'Laser On');
annotation('rectangle',[0.155 0.926 0.154 0.008], 'FaceColor', cols.opto, 'LineStyle', 'none');
% title('All Cells', 'Units', 'normalized', 'Position', [0.5 1.03]);
fixfig;
drawnow;
hold off;

% Scatter Plot
subplot(2,3,2); hold on;

ER_OptoOn_Clean = DataOut_OG_GAD2.ClusterData.ER_OptoOn;
ER_OptoOn_Clean(DataOut_OG_GAD2.ClusterData.ER_OptoOn > 50) = 30;

plot([-20:1:100],[-20:1:100],'--', 'Color', cols.refline); % Reference Line
scatter(DataOut_OG_GAD2.ClusterData.ER_OptoOff, ER_OptoOn_Clean, 'MarkerEdgeColor', cols.off);
ylim([-5 30]);
xlim([-5 30]);
yticks([0 10 20 30]);
yticklabels(["0" "10" "20" ">30"]);
ylabel('Evoked Rate Opto On (spks/s)');
xlabel('Evoked Rate Opto Off (spks/s)');
grid on
% axis square
fixfig;
hold off

% Barplot
outliers_GAD2 = abs(zscore(DataOut_OG_GAD2.ClusterData.ER_OptoOn)) > 3;

subplot(2,3,3); hold on;
bar_pl = bar([DataOut_OG_GAD2.Overall.ER_OptoOff DataOut_OG_GAD2.Overall.ER_OptoOn],'FaceColor', [cols.off]);
bar_pl.FaceColor = 'flat';
bar_pl.CData(2,:) = cols.on_GAD2;
bar_err = errorbar([1,2],[DataOut_OG_GAD2.Overall.ER_OptoOff DataOut_OG_GAD2.Overall.ER_OptoOn], [DataOut_OG_GAD2.Overall.SE_OptoOff DataOut_OG_GAD2.Overall.SE_OptoOn]);
bar_err.Color = cols.err_bar;
bar_err.LineStyle = 'none';
bar_err.LineWidth = 1.25;
sigstar({[1 2]}, DataOut_OG_GAD2.Overall.p_val);
ylabel('Evoked Rate (spks/s)');
ylim([0 4.5]);
xticks([1 2]);
xticklabels({"Off", "On"});
% xtickangle(15);
fixfig;
hold off

%% --- Comparing Response Shape

% --- Example Cell Downregulated ---
cell_sel = [(DataOut_OG_GAD2.ClusterDataDown.Subject + "_" + DataOut_OG_GAD2.ClusterDataDown.RecDate + "_AP.mat") ...
    DataOut_OG_GAD2.ClusterDataDown.ClusterN];
idx = find(DataOut_OG_GAD2.ClusterDataDown.p_val == min(DataOut_OG_GAD2.ClusterDataDown.p_val));

% A. Rasterplot
% Load Stim Data
sAP1 = load(cell_sel(idx,1)); sAP1 = sAP1.sAP;
st = sAP1.sCluster(str2num(cell_sel(idx,2))).SpikeTimes;
laser_trials = sAP1.cellBlock{1}.vecOptoOn;
stim = sAP1.cellBlock{1}.vecStimOnTime - 0.2;
maxDur = 1.4;

%Draw Plot
subplot(4,4,9); hold on;
plotRasterSplit(st,stim,maxDur,~laser_trials,[cols.Inh; cols.off]); 
xline(0.2, '--', 'Color', cols.xline); xline(1.2, '--', 'Color', cols.xline);
% xline(0.170, '-g'); xline(1.18,'-g');
xlim([0 1.4]); xlabel("Time from stim. onset (s)");
xticks([0.2 0.7 1.2]);
xticklabels(["0" "0.5" "1"]);
annotation('rectangle',[0.148 0.488 0.113 0.008], 'FaceColor', cols.opto, 'LineStyle', 'none');
fixfig; hold off;

% B. PSTH
subplot(4,4,13); hold on;
plot(DataOut_OG_GAD2.ClusterDataDown.PSTHBinCenters_Off(idx,:), DataOut_OG_GAD2.ClusterDataDown.PSTHMean_Off(idx,:),'Color', cols.off);
plot(DataOut_OG_GAD2.ClusterDataDown.PSTHBinCenters_Off(idx,:), DataOut_OG_GAD2.ClusterDataDown.PSTHMean_On(idx,:),'Color',cols.Inh);
xlim([-0.2 1.2]); xlabel("Time from stim. onset (s)"); ylabel("Spiking Rate (spks/s)");
xline(0, '--', 'Color', cols.xline); xline(1, '--', 'Color', cols.xline);
% xline(-0.03, '-g'); xline(0.98, 'g');
annotation('rectangle',[0.148 0.269 0.113 0.008], 'FaceColor', cols.opto, 'LineStyle', 'none');
fixfig; hold off;

% Mean Response Downregulated Cells
PSTHMeanDown_Off_BL = DataOut_OG_GAD2.ClusterDataDown.PSTHMean_Off - DataOut_OG_GAD2.ClusterDataDown.SpontRate;
PSTHMeanDown_On_BL = DataOut_OG_GAD2.ClusterDataDown.PSTHMean_On - DataOut_OG_GAD2.ClusterDataDown.SpontRate;
norm = max(PSTHMeanDown_Off_BL,[],2);
PSTHMeanDown_Off_norm = mean(PSTHMeanDown_Off_BL./norm, 1);
PSTHSEMDown_Off_norm = std(PSTHMeanDown_Off_BL./norm, 1)/sqrt(DataOut_OG_GAD2.Overall.NCells_Reduced);
PSTHMeanDown_On_norm = mean(PSTHMeanDown_On_BL./norm, 1);
PSTHSEMDown_On_norm = std(PSTHMeanDown_On_BL./norm, 1)/sqrt(DataOut_OG_GAD2.Overall.NCells_Reduced);

subplot(2,4,6); hold on;
shadedErrorBar(DataOut_OG_GAD2.Overall.PSTHBinCenters, PSTHMeanDown_Off_norm, ...
    PSTHSEMDown_Off_norm, 'lineProps', {'Color', cols.off});
shadedErrorBar(DataOut_OG_GAD2.Overall.PSTHBinCenters, PSTHMeanDown_On_norm, ...
    PSTHSEMDown_On_norm, 'lineProps', {'Color', cols.Inh});
ylabel('Normalized Response');
xlabel('Time from stim. onset (s)');
xlim([min(DataOut_OG_GAD2.Overall.PSTHtime) max(DataOut_OG_GAD2.Overall.PSTHtime)]);
% % ylim([-0.1 1]);
xline(0, '--', 'Color', cols.xline); xline(1, '--', 'Color', cols.xline);
% xline(-0.03,'g-'); xline(0.980,'g');
yline(0, 'Color', cols.xline);
text(0.82, 0.56, sprintf('%g units', DataOut_OG_GAD2.Overall.NCells_Reduced), 'FontSize',15);
legend('Laser On', 'Laser Off');
% title("Downregulated Cells",'Units', 'normalized', "Position",[0.5 1.04])
annotation('rectangle',[0.354 0.453 0.113 0.008], 'FaceColor', cols.opto, 'LineStyle', 'none');
fixfig;
drawnow;
hold off;

% --- Example Cell Upregulated ---
cell_sel = [(DataOut_OG_GAD2.ClusterDataUp.Subject + "_" + DataOut_OG_GAD2.ClusterDataUp.RecDate + "_AP.mat") ...
    DataOut_OG_GAD2.ClusterDataUp.ClusterN];
idx = find(DataOut_OG_GAD2.ClusterDataUp.p_val == min(DataOut_OG_GAD2.ClusterDataUp.p_val));

% A. Rasterplot
% Load Stim Data
sAP1 = load(cell_sel(idx,1)); sAP1 = sAP1.sAP;
st = sAP1.sCluster(str2num(cell_sel(idx,2))).SpikeTimes;
laser_trials = sAP1.cellBlock{1}.vecOptoOn;
stim = sAP1.cellBlock{1}.vecStimOnTime - 0.2;
maxDur = 1.4;

%Draw Plot
subplot(4,4,11); hold on;
plotRasterSplit(st,stim,maxDur,~laser_trials,[cols.Act; cols.off]); 
xline(0.2, '--', 'Color', cols.xline); xline(1.2, '--', 'Color', cols.xline);
% xline(0.170, '-g'); xline(1.18,'-g');
xlim([0 1.4]); xlabel("Time from stim. onset (s)");
xticks([0.2 0.7 1.2]);
xticklabels(["0" "0.5" "1"]);
annotation('rectangle',[0.560 0.488 0.113 0.008], 'FaceColor', cols.opto, 'LineStyle', 'none');
fixfig; hold off;

% B. PSTH
subplot(4,4,15); hold on;
plot(DataOut_OG_GAD2.ClusterDataUp.PSTHBinCenters_Off(idx,:), DataOut_OG_GAD2.ClusterDataUp.PSTHMean_Off(idx,:),'Color', cols.off);
plot(DataOut_OG_GAD2.ClusterDataUp.PSTHBinCenters_Off(idx,:), DataOut_OG_GAD2.ClusterDataUp.PSTHMean_On(idx,:),'Color',cols.Act);
xlim([-0.2 1.2]); xlabel("Time from stim. onset (s)"); ylabel("Spiking Rate (spks/s)");
xline(0, '--', 'Color', cols.xline); xline(1, '--', 'Color', cols.xline);
% xline(-0.03, '-g'); xline(0.98, 'g');
annotation('rectangle',[0.560 0.269 0.113 0.008], 'FaceColor', cols.opto, 'LineStyle', 'none');
fixfig; hold off;


% --- Mean Response Upregulated Cells ---
PSTHMeanUp_Off_BL = DataOut_OG_GAD2.ClusterDataUp.PSTHMean_Off - DataOut_OG_GAD2.ClusterDataUp.SpontRate;
PSTHMeanUp_On_BL = DataOut_OG_GAD2.ClusterDataUp.PSTHMean_On - DataOut_OG_GAD2.ClusterDataUp.SpontRate;
norm = max(PSTHMeanUp_Off_BL,[],2);
PSTHMeanUp_Off_norm = mean(PSTHMeanUp_Off_BL./norm, 1);
PSTHSEMUp_Off_norm = std(PSTHMeanUp_Off_BL./norm, 1)/sqrt(DataOut_OG_GAD2.Overall.NCells_Increased);
PSTHMeanUp_On_norm = mean(PSTHMeanUp_On_BL./norm, 1);
PSTHSEMUp_On_norm = std(PSTHMeanUp_On_BL./norm, 1)/sqrt(DataOut_OG_GAD2.Overall.NCells_Increased);

subplot(2,4,8); hold on;
shadedErrorBar(DataOut_OG_GAD2.Overall.PSTHBinCenters, PSTHMeanUp_Off_norm, ...
    PSTHSEMUp_Off_norm, 'lineProps', {'Color', cols.off});
shadedErrorBar(DataOut_OG_GAD2.Overall.PSTHBinCenters, PSTHMeanUp_On_norm, ...
    PSTHSEMUp_On_norm, 'lineProps', {'Color', cols.Act});
ylabel('Normalized Response');
xlabel('Time from stim. onset (s)');
xlim([min(DataOut_OG_GAD2.Overall.PSTHtime) max(DataOut_OG_GAD2.Overall.PSTHtime)]);
% % ylim([-0.1 1]);
xline(0, '--', 'Color', cols.xline); xline(1, '--', 'Color', cols.xline);
% xline(-0.03,'g-'); xline(0.980,'g');
yline(0, 'Color', cols.xline);
text(0.82, 7, sprintf('%g units', DataOut_OG_GAD2.Overall.NCells_Increased), 'FontSize',15);
legend('Laser On', 'Laser Off');
% title("Upregulated Cells", 'Units', 'normalized', "Position", [0.5 1.04])
annotation('rectangle',[0.766 0.453 0.113 0.008], 'FaceColor', cols.opto, 'LineStyle', 'none');
fixfig;
drawnow;
hold off;

%% Save Figure Output

% saveas(gcf, 'C:\Software and Code\Analysis-GAD2\Plots\Figure1.png');
% savefig(gcf, 'C:\Software and Code\Analysis-GAD2\Plots\Figure1.fig');
saveas(gcf, 'D:\NIN\Analysis-GAD2\Plots\Figure2.png');
saveas(gcf, 'D:\NIN\Analysis-GAD2\Plots\Figure2.fig');

%% Figure S1. Control Mice

figure;

% --- Example Cell Controls ---
cell_sel = [(DataOut_OG_Cont.ClusterDataNonSig.Subject + "_" + DataOut_OG_Cont.ClusterDataNonSig.RecDate + "_AP.mat") ...
    DataOut_OG_Cont.ClusterDataNonSig.ClusterN];
idx = 17;

% A. Rasterplot
% Load Stim Data
sAP2 = load(cell_sel(idx,1)); sAP2 = sAP2.sAP;
st = sAP2.sCluster(str2num(cell_sel(idx,2))).SpikeTimes;
laser_trials = sAP2.cellBlock{1}.vecOptoOn;
stim = sAP2.cellBlock{1}.vecStimOnTime - 0.2;
maxDur = 1.4;

%Draw Plot
subplot(2,4,1); hold on;
plotRasterSplit(st,stim,maxDur,~laser_trials,[cols.on_Cont - 0.1; cols.off]); 
xline(0.2, '--', 'Color', cols.xline); xline(1.2, '--', 'Color', cols.xline);
xlim([0 1.4]); xlabel("Time from stim. onset (s)");
xticks([0.2 0.7 1.2]);
xticklabels(["0" "0.5" "1"]);
annotation('rectangle',[0.148 0.926 0.114 0.008], 'FaceColor', cols.opto, 'LineStyle', 'none');
fixfig; hold off;

% B. PSTH
subplot(2,4,5); hold on;
plot(DataOut_OG_Cont.ClusterDataNonSig.PSTHBinCenters_Off(idx,:), DataOut_OG_Cont.ClusterDataNonSig.PSTHMean_Off(idx,:),'Color', cols.off);
plot(DataOut_OG_Cont.ClusterDataNonSig.PSTHBinCenters_Off(idx,:), DataOut_OG_Cont.ClusterDataNonSig.PSTHMean_On(idx,:),'Color',cols.on_Cont);
xlim([-0.2 1.2]); xlabel("Time from stim. onset (s)"); ylabel('Spiking Rate (spks/s)');
xline(0, '--', 'Color', cols.xline); xline(1, '--', 'Color', cols.xline);
annotation('rectangle',[0.148 0.453 0.114 0.008], 'FaceColor', cols.opto, 'LineStyle', 'none');
fixfig; hold off;

% --- Summary Stats Control Mice ---
% Mean PSTH - Normalized to Peak
subplot(1,4,2); hold on;
shadedErrorBar(DataOut_OG_Cont.Overall.PSTHBinCenters, DataOut_OG_Cont.Overall.PSTHMean_Off_norm, DataOut_OG_Cont.Overall.PSTHSEM_Off_norm, 'lineProps', {'Color', cols.off});
shadedErrorBar(DataOut_OG_Cont.Overall.PSTHBinCenters, DataOut_OG_Cont.Overall.PSTHMean_On_norm, DataOut_OG_Cont.Overall.PSTHSEM_On_norm, 'lineProps', {'Color', cols.on_Cont});
xline(0,'--','Color',cols.xline); xline(1,'--','Color',cols.xline);
% xline(-0.03, 'g'); xline(0.98, 'g');
ylabel('Normalized Response');
xlabel('Time from stim. onset (s)');
xlim([min(DataOut_OG_Cont.Overall.PSTHtime) max(DataOut_OG_Cont.Overall.PSTHtime)]);
% ylim([-0.1 1]);
yline(0, 'Color', cols.xline);
text(0.6, 2, sprintf('%g units', DataOut_OG_Cont.Overall.NCells), 'FontSize',15);
annotation('rectangle',[0.354 0.926 0.114 0.008], 'FaceColor', cols.opto, 'LineStyle', 'none');
legend('Laser Off', 'Laser On');
fixfig;
drawnow;
hold off;

% Scatter Plot
subplot(1,4,3); hold on;
plot([-20:1:100],[-20:1:100],'--', 'Color', cols.refline); % Reference Line
scatter(DataOut_OG_Cont.ClusterData.ER_OptoOff, DataOut_OG_Cont.ClusterData.ER_OptoOn, 'MarkerEdgeColor', cols.off);
ylim([-5 30]);
xlim([-5 30]);
ylabel('Evoked Rate Opto On (spks/s)');
xlabel('Evoked Rate Opto Off (spks/s)');
grid on
% axis square
fixfig;
hold off

% Barplot
subplot(1,4,4); hold on;
bar_pl = bar([DataOut_OG_Cont.Overall.ER_OptoOff DataOut_OG_Cont.Overall.ER_OptoOn],'FaceColor', [cols.off]);
bar_pl.FaceColor = 'flat';
bar_pl.CData(2,:) = cols.on_Cont;

outliers_cont = abs(zscore(DataOut_OG_Cont.ClusterData.ER_OptoOn)) > 3;
bar_err = errorbar([1,2],[DataOut_OG_Cont.Overall.ER_OptoOff DataOut_OG_Cont.Overall.ER_OptoOn], [DataOut_OG_Cont.Overall.SE_OptoOff DataOut_OG_Cont.Overall.SE_OptoOn]);
bar_err.Color = cols.err_bar_Cont;
bar_err.LineStyle = 'none';
bar_err.LineWidth = 1.25;
ylabel('Evoked Rate (spks/s)');
xticks([1 2]);
xticklabels({"Off", "On"});
% xtickangle(15);
ylim([0 3.5]);
fixfig;
hold off

%% Unpaired T-tests

[~, pval_ER_OptoOff] = ttest2(DataOut_OG_GAD2.ClusterData.ER_OptoOff(~outliers_GAD2),DataOut_OG_Cont.ClusterData.ER_OptoOff(~outliers_cont));
[~, pval_ER_OptoOn] = ttest2(DataOut_OG_GAD2.ClusterData.ER_OptoOn(~outliers_GAD2),DataOut_OG_Cont.ClusterData.ER_OptoOn(~outliers_cont));

%% Save
saveas(gcf, 'D:\NIN\Analysis-GAD2\Plots\S_Fig1.png');
saveas(gcf, 'D:\NIN\Analysis-GAD2\Plots\S_Fig1.fig');

