%% Figure 1: Plots Opto Gratings

% TO-DO:
% [] - Remove Outliers from Scatterplot GAD2

%% Colors

cols.off = [0 0 0];
cols.on = [0 0.4470 0.7410];
cols.xline = [143 143 143]/255;
cols.err_bar = [26 35 126]/255;

figure;
%% A-D: Invidiual Example Neurons
 
% --- Example Cell GAD2 ---
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
subplot(4,4,1); hold on;
plotRasterSplit(st,stim,maxDur,~laser_trials,[cols.on; cols.off]); xline(0.2, '--', 'Color', cols.xline);
xlim([0 1.4]); xlabel("Time from stimulus onset (s)"); ylabel("Spiking Rate (spks/s)");
% xline(0.110, 'g--'); xline(0.130, 'g--');
% title(DataOut_OG_GAD2.ClusterDataDown.Subject(idx) + "; Cl: " + cell_sel(idx,2) + " Opto Off");
fixfig; hold off;

% B. PSTH
subplot(4,4,5); hold on;
plot(DataOut_OG_GAD2.ClusterDataDown.PSTHBinCenters_Off(idx,:), DataOut_OG_GAD2.ClusterDataDown.PSTHMean_Off(idx,:),'Color', cols.off);
plot(DataOut_OG_GAD2.ClusterDataDown.PSTHBinCenters_Off(idx,:), DataOut_OG_GAD2.ClusterDataDown.PSTHMean_On(idx,:),'Color',cols.on);
xlim([-0.2 1.2]); xlabel("Time from stimulus onset (s)")
xline(0, '--', 'Color', cols.xline);
fixfig; hold off;

% --- Example Cell Controls ---
cell_sel = [(DataOut_OG_Cont.ClusterDataNonSig.Subject + "_" + DataOut_OG_Cont.ClusterDataNonSig.RecDate + "_AP.mat") ...
    DataOut_OG_Cont.ClusterDataNonSig.ClusterN];
idx = 10;

% A. Rasterplot
% Load Stim Data
sAP2 = load(cell_sel(idx,1)); sAP2 = sAP2.sAP;
st = sAP2.sCluster(str2num(cell_sel(idx,2))).SpikeTimes;
laser_trials = sAP2.cellBlock{1}.vecOptoOn;
stim = sAP2.cellBlock{1}.vecStimOnTime - 0.2;
maxDur = 1.4;

%Draw Plot
subplot(4,4,9); hold on;
plotRasterSplit(st,stim,maxDur,~laser_trials,[cols.on; cols.off]); xline(0.2, '--', 'Color', cols.xline);
xlim([0 1.4]); xlabel("Time from stimulus onset (s)");
% xline(0.110, 'g--'); xline(0.130, 'g--');
% title(DataOut_OG_GAD2.ClusterDataDown.Subject(idx) + "; Cl: " + cell_sel(idx,2) + " Opto Off");
fixfig; hold off;

% B. PSTH
subplot(4,4,13); hold on;
plot(DataOut_OG_Cont.ClusterDataNonSig.PSTHBinCenters_Off(idx,:), DataOut_OG_Cont.ClusterDataNonSig.PSTHMean_Off(idx,:),'Color', cols.off);
plot(DataOut_OG_Cont.ClusterDataNonSig.PSTHBinCenters_Off(idx,:), DataOut_OG_Cont.ClusterDataNonSig.PSTHMean_On(idx,:),'Color',cols.on);
xlim([-0.2 1.2]); xlabel("Time from stimulus onset (s)"); ylabel('Spiking Rate (Spks/S)');
xline(0, '--', 'Color', cols.xline);
fixfig; hold off;

%% Leonie's Figures pt. 2: Summary Stats

% --- GAD2 Mice ---
% TabNoUp = [DataOut_OG_GAD2.ClusterDataDown; DataOut_OG_GAD2.ClusterDataNonSig];

% Mean Response Summary GAD2; Z-Score
subplot(2,4,2); hold on;
shadedErrorBar(DataOut_OG_GAD2.Overall.PSTHBinCenters, DataOut_OG_GAD2.Overall.PSTHMean_Off_z, DataOut_OG_GAD2.Overall.PSTHSEM_Off_z, 'lineProps', {'Color', cols.off});
shadedErrorBar(DataOut_OG_GAD2.Overall.PSTHBinCenters, DataOut_OG_GAD2.Overall.PSTHMean_On_z, DataOut_OG_GAD2.Overall.PSTHSEM_On_z, 'lineProps', {'Color', cols.on});
xline(0,'--','Color',cols.xline);
ylabel('Normalized Response (Z-Score)');
xlabel('Time from stimulus onset (s)');
xlim([min(DataOut_OG_GAD2.Overall.PSTHtime) max(DataOut_OG_GAD2.Overall.PSTHtime)]);
% ylim([-0.1 1]);
yline(0);
text(0.6, 2, sprintf('%g units', DataOut_OG_GAD2.Overall.NCells), 'FontSize',15);
legend('No Opto', 'Opto');
% title('Mean Response Overall');
fixfig;
drawnow;
hold off;

% Scatter Plot
subplot(2,4,3); hold on;
scatter(DataOut_OG_GAD2.ClusterData.ER_OptoOff, DataOut_OG_GAD2.ClusterData.ER_OptoOn, 'MarkerEdgeColor', cols.off);
plot([-20:1:100],[-20:1:100],'--', 'Color', cols.xline); % Reference Line
% ylim([-15 100]);
% xlim([-15 100]);
ylabel('Evoked Rate Opto On (spks/s)');
xlabel('Evoked Rate Opto Off (spks/s)');
grid on
% axis square
fixfig;
hold off

% Barplot
% figure; subplot(2,4,4); hold on;
% bar_pl = bar([DataOut_OG_GAD2.Overall.ER_OptoOff DataOut_OG_GAD2.Overall.ER_OptoOn],'FaceColor', [cols.off]);
% bar_pl.FaceColor = 'flat';
% bar_pl.CData(2,:) = cols.on;
% % dot_pl = plot([1 2],[DataOut_OG_GAD2.ClusterData.ER_OptoOff DataOut_OG_GAD2.ClusterData.ER_OptoOn],'-','Color',[0, 0, 0, 0.3]);
% % scatter([1 2],[DataOut_OG_GAD2.ClusterData.ER_OptoOff DataOut_OG_GAD2.ClusterData.ER_OptoOn],'MarkerEdgeColor', [0 0 0], 'MarkerEdgeAlpha',0.5);
% bar_err = errorbar([1,2],[DataOut_OG_GAD2.Overall.ER_OptoOff DataOut_OG_GAD2.Overall.ER_OptoOn], [DataOut_OG_GAD2.Overall.SE_OptoOff DataOut_OG_GAD2.Overall.SE_OptoOn]);
% bar_err.Color = [26 35 126]/255;
% bar_err.LineStyle = 'none';
% bar_err.LineWidth = 1.25;
% ylabel('Evoked Rate (spks/s)');
% xticks([1 2]);
% xticklabels({"Laser OFF", "Laser ON"});
% fixfig;
% hold off

% Barplot w/o Upregulated Cells
ER_OptoOn = mean([DataOut_OG_GAD2.ClusterDataDown.ER_OptoOn; DataOut_OG_GAD2.ClusterDataNonSig.ER_OptoOn]);
ER_OptoOff = mean([DataOut_OG_GAD2.ClusterDataDown.ER_OptoOff; DataOut_OG_GAD2.ClusterDataNonSig.ER_OptoOff]);
SE_OptoOn = std([DataOut_OG_GAD2.ClusterDataDown.ER_OptoOn; DataOut_OG_GAD2.ClusterDataNonSig.ER_OptoOn], [], 1)/sqrt(height(DataOut_OG_GAD2.ClusterDataDown));
SE_OptoOff = std([DataOut_OG_GAD2.ClusterDataDown.ER_OptoOff; DataOut_OG_GAD2.ClusterDataNonSig.ER_OptoOff], [], 1)/sqrt(height(DataOut_OG_GAD2.ClusterDataDown));

subplot(2,4,4); hold on;
bar_pl = bar([ER_OptoOff ER_OptoOn],'FaceColor', [cols.off]);
bar_pl.FaceColor = 'flat';
bar_pl.CData(2,:) = cols.on;
% dot_pl = plot([1 2],[DataOut_OG_GAD2.ClusterData.ER_OptoOff DataOut_OG_GAD2.ClusterData.ER_OptoOn],'-','Color',[0, 0, 0, 0.3]);
% scatter([1 2],[DataOut_OG_GAD2.ClusterData.ER_OptoOff DataOut_OG_GAD2.ClusterData.ER_OptoOn],'MarkerEdgeColor', [0 0 0], 'MarkerEdgeAlpha',0.5);
bar_err = errorbar([1,2],[ER_OptoOff ER_OptoOn], [SE_OptoOff SE_OptoOn]);
bar_err.Color = cols.err_bar;
bar_err.LineStyle = 'none';
bar_err.LineWidth = 1.25;
ylabel('Evoked Rate (spks/s)');
xticks([1 2]);
xticklabels({"Laser OFF", "Laser ON"});
fixfig;
hold off

% --- Control Mice ---
% Mean PSTH - Z-Scored
subplot(2,4,6); hold on;
shadedErrorBar(DataOut_OG_Cont.Overall.PSTHBinCenters, DataOut_OG_Cont.Overall.PSTHMean_Off_z, DataOut_OG_Cont.Overall.PSTHSEM_Off_z, 'lineProps', {'Color', cols.off});
shadedErrorBar(DataOut_OG_Cont.Overall.PSTHBinCenters, DataOut_OG_Cont.Overall.PSTHMean_On_z, DataOut_OG_Cont.Overall.PSTHSEM_On_z, 'lineProps', {'Color', cols.on});
xline(0,'--','Color',cols.xline);
ylabel('Normalized Response (Z-Score)');
xlabel('Time from stimulus onset (s)');
xlim([min(DataOut_OG_Cont.Overall.PSTHtime) max(DataOut_OG_Cont.Overall.PSTHtime)]);
% ylim([-0.1 1]);
yline(0);
text(0.6, 2, sprintf('%g units', DataOut_OG_Cont.Overall.NCells), 'FontSize',15);
legend('No Opto', 'Opto');
% title('Mean Response Overall');
fixfig;
drawnow;
hold off;

% Scatter Plot
subplot(2,4,7); hold on;
scatter(DataOut_OG_Cont.ClusterData.ER_OptoOff, DataOut_OG_Cont.ClusterData.ER_OptoOn, 'MarkerEdgeColor', cols.off);
plot([-20:1:100],[-20:1:100],'--', 'Color', cols.xline); % Reference Line
ylim([-5 30]);
xlim([-5 30]);
ylabel('Evoked Rate Opto On (spks/s)');
xlabel('Evoked Rate Opto Off (spks/s)');
grid on
% axis square
fixfig;
hold off

% Barplot
subplot(2,4,8); hold on;
bar_pl = bar([DataOut_OG_Cont.Overall.ER_OptoOff DataOut_OG_Cont.Overall.ER_OptoOn],'FaceColor', [cols.off]);
bar_pl.FaceColor = 'flat';
bar_pl.CData(2,:) = cols.on;
bar_err = errorbar([1,2],[DataOut_OG_Cont.Overall.ER_OptoOff DataOut_OG_Cont.Overall.ER_OptoOn], [DataOut_OG_Cont.Overall.SE_OptoOff DataOut_OG_Cont.Overall.SE_OptoOn]);
bar_err.Color = cols.err_bar;
bar_err.LineStyle = 'none';
bar_err.LineWidth = 1.25;
ylabel('Evoked Rate (spks/s)');
xticks([1 2]);
xticklabels({"Laser OFF", "Laser ON"});
fixfig;
hold off

%% Supplementary Plots for Leonie

% % Mean Response Downregulated Cells
% PSTHoff_Down_z = (DataOut_OG_GAD2.ClusterDataDown.PSTHMean_Off - mean(DataOut_OG_GAD2.ClusterDataDown.SpontRate))./std(DataOut_OG_GAD2.ClusterDataDown.SpontRate);
% PSTHon_Down_z = (DataOut_OG_GAD2.ClusterDataDown.PSTHMean_On - mean(DataOut_OG_GAD2.ClusterDataDown.SpontRate))./std(DataOut_OG_GAD2.ClusterDataDown.SpontRate);
% 
% PSTHMean_Off_z = mean(PSTHoff_Down_z, 1);
% PSTHSEM_Off_z = std(PSTHoff_Down_z, 1)/sqrt(height(DataOut_OG_GAD2.ClusterDataUp));
% PSTHMean_On_z = mean(PSTHon_Down_z, 1);
% PSTHSEM_On_z = std(PSTHon_Down_z, 1)/sqrt(height(DataOut_OG_GAD2.ClusterDataUp));
% 
% figure; hold on;
% shadedErrorBar(DataOut_OG_GAD2.Overall.PSTHBinCenters, mean(PSTHoff_Down_z, 1), ...
%     std(PSTHoff_Down_z, 1)/sqrt(height(DataOut_OG_GAD2.ClusterDataUp)), 'lineProps', 'k');
% shadedErrorBar(DataOut_OG_GAD2.Overall.PSTHBinCenters, mean(PSTHon_Down_z, 1), ...
%     std(PSTHon_Down_z, 1)/sqrt(height(DataOut_OG_GAD2.ClusterDataUp)), 'lineProps', 'b');
% xline(0,'r--');
% ylabel('Normalized Response (Z-Score)');
% xlabel('Time from stimulus onset (s)');
% xlim([min(DataOut_OG_GAD2.Overall.PSTHtime) max(DataOut_OG_GAD2.Overall.PSTHtime)]);
% % ylim([-0.1 1]);
% yline(0);
% text(0.94, 2.2, sprintf('%g units', DataOut_OG_GAD2.Overall.NCells_Reduced), 'FontSize',15);
% legend('No Opto', 'Opto');
% title('Downregulated Cells');
% fixfig;
% drawnow;
% hold off;
