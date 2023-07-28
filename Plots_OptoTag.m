% Plots for Optotagging 
DOT = DataOut_OT_50ms.Overall;

GadTab = DataOut_OT_50ms.ClusterData(DataOut_OT_50ms.ClusterData.CellClass == "GAD2+",:);
ActTab = DataOut_OT_50ms.ClusterData(DataOut_OT_50ms.ClusterData.CellClass == "Activated",:);
InhTab = DataOut_OT_50ms.ClusterData(DataOut_OT_50ms.ClusterData.CellClass == "Inhibited",:);

%% Colors

cols.off = [0 0 0];
cols.GAD2 = [0 0.4470 0.7410];
cols.Act = [0.8500 0.3250 0.0980];
cols.Inh = [0.9290 0.6940 0.1250];
% cols.Oth
cols.xline = [143 143 143]/255;
cols.xreg = [200 200 200]/255;
cols.xreg2 = [184 202 214]/255;
cols.err_bar = [26 35 126]/255;

figure;
%% Fig 3A-C: Example RasterPlots

pulseSelect = 0.05

% --- GAD2+ Cell ---
sAP = load("79155_20230512_AP.mat"); sAP = sAP.sAP;
GAD2_Ex = 353;

st = sAP.sCluster(GAD2_Ex).SpikeTimes;
stim = sAP.cellBlock{1, 3}.vecLaserOnTime(sAP.cellBlock{1, 3}.PulseDurTrial == pulseSelect) - 0.05;
maxDur = 0.15;

subplot(3,3,1); hold on;
plotRaster(st,stim,maxDur, [], cols.GAD2);
xregion(0.05, pulseSelect+0.05, 'FaceColor', cols.xreg) % xline(0.060, 'g--'); xline(0.080, 'g--');
xticks([0 0.05 0.1 0.15]);
xticklabels(["-0.05" "0" "0.05" "0.1"]);
xlabel('Time from laser onset (s)');
% title(string(sAP.sJson.subject) + " Cl: " + string(GAD2_Ex) + " Pulse: " + pulseSelect);
title("GAD2+");
fixfig; hold off;

% --- Inh Cell ---
sAP = load("79154_20230511_AP.mat"); sAP = sAP.sAP;
Inh_Ex = 446;

st = sAP.sCluster(Inh_Ex).SpikeTimes;
stim = sAP.cellBlock{1, 3}.vecLaserOnTime(sAP.cellBlock{1, 3}.PulseDurTrial == pulseSelect) - 0.05;
maxDur = 0.15;

subplot(3,3,4); hold on;
plotRaster(st,stim,maxDur, [], cols.Inh); 
xregion(0.05, pulseSelect+0.05, 'FaceColor', cols.xreg);
xticks([0 0.05 0.1 0.15]);
xticklabels(["-0.05" "0" "0.05" "0.1"]);
xlabel('Time from laser onset (s)');
% title(string(sAP.sJson.subject) + " Cl: " + string(Inh_Ex) + " Pulse: " + pulseSelect);
title("Inhibited");
fixfig; hold off;

% --- Act Cell ---
sAP = load("79159_20230426_AP.mat"); sAP = sAP.sAP;
Act_Ex = 377;

st = sAP.sCluster(Act_Ex).SpikeTimes;
stim = sAP.cellBlock{1, 3}.vecLaserOnTime(sAP.cellBlock{1, 3}.PulseDurTrial == pulseSelect) - 0.05;
maxDur = 0.15;

subplot(3,3,7); hold on;
plotRaster(st,stim,maxDur, [], cols.Act);
xregion(0.05, pulseSelect+0.05, 'FaceColor', cols.xreg);
xticks([0 0.05 0.1 0.15]);
xticklabels(["-0.05" "0" "0.05" "0.1"]);
xlabel('Time from laser onset (s)');
% title(string(sAP.sJson.subject) + " Cl: " + string(Act_Ex) + " Pulse: " + pulseSelect);
title("Activated");
fixfig; hold off;

%% Fig 3 D-E: Overall Plots

BinCenters = DOT.PSTHBinCenters(1,:);
excl_idx = (BinCenters > -0.001 & BinCenters < 0.002) | (BinCenters > 0.049 & BinCenters < 0.052);

% Heatmap (-0.05s - 0.1s); Organized by 1. GAD2, 2. Inhibited, 3. Activated
ReorderedTab = [GadTab; InhTab; ActTab];
PSTH_BL = ReorderedTab.PSTHMean - ReorderedTab.SpontRate;
PSTH_norm = PSTH_BL./max(PSTH_BL,[],2);
PSTH_norm(:,excl_idx) = NaN;
subplot(2,3,[2 3]); imagesc(PSTH_norm);
colormap('bone');
colorbar(gca);
y_label = ylabel('Neuron #');
y_label.Position(1) = 23;
set(gca, 'YTick', []);
xticks([1 26 51 76 100]);
xticklabels(["-0.1" "-0.05" "0" "0.05" "0.1"]);
xlim([26 100]);
xlabel('Time from laser onset (s)');
annotation('rectangle',[0.4 0.88 0.0075 0.045], 'FaceColor', cols.GAD2, 'LineStyle', 'none')
annotation('rectangle',[0.4 0.675 0.0075 0.202], 'FaceColor', cols.Inh, 'LineStyle', 'none')
annotation('rectangle',[0.4 0.585 0.0075 0.088], 'FaceColor', cols.Act, 'LineStyle', 'none')

% PSTH (Baseline-Subtracted; Norm. to Peak)
PSTHMean_GAD2_norm = DOT.PSTHMean_GAD2_norm;
PSTHMean_GAD2_norm(excl_idx) = NaN;
PSTHSEM_GAD2_norm = DOT.PSTHSEM_GAD2_norm;
PSTHSEM_GAD2_norm(excl_idx) = NaN;

PSTHMean_Inh_norm = DOT.PSTHMean_Inh_norm;
PSTHMean_Inh_norm(excl_idx) = NaN;
PSTHSEM_Inh_norm = DOT.PSTHSEM_Inh_norm;
PSTHSEM_Inh_norm(excl_idx) = NaN;

PSTHMean_Act_norm = DOT.PSTHMean_Act_norm;
PSTHMean_Act_norm(excl_idx) = NaN;
PSTHSEM_Act_norm = DOT.PSTHSEM_Act_norm;
PSTHSEM_Act_norm(excl_idx) = NaN;

subplot(2,3,[5 6]); hold on;
shadedErrorBar(DOT.PSTHBinCenters, PSTHMean_GAD2_norm, PSTHSEM_GAD2_norm, 'lineProps', {'Color',cols.GAD2});
shadedErrorBar(DOT.PSTHBinCenters, PSTHMean_Inh_norm, PSTHSEM_Inh_norm, 'lineProps', {'Color',cols.Inh});
shadedErrorBar(DOT.PSTHBinCenters, PSTHMean_Act_norm, PSTHSEM_Act_norm, 'lineProps', {'Color',cols.Act});
xregion(0,pulseSelect,"FaceColor",cols.xreg,"LineStyle","none");
ylim([-0.2 1.2]);
xticks([-0.1 -0.05 0 0.05 0.1]);
xlim([-0.05 0.1]);
legend(["GAD2+ (N = "+height(GadTab)+")", "Inhibited (N = "+height(InhTab)+")", "Activated (N = "+height(ActTab)+")"]);
ylabel('Mean Response (Norm. to Peak)');
xlabel('Time from laser onset (ms)');
% title('BL-Subtracted and Normalized to Peak');
fixfig; hold off;

%% Save Figure

saveas(gcf, 'C:\Software and Code\Analysis-GAD2\Plots\Figure2.png');
savefig(gcf, 'C:\Software and Code\Analysis-GAD2\Plots\Figure2.fig');


%% Suppl. Figure 1

% Maybe some plots comparing different pulse durations?