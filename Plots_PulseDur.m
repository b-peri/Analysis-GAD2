% Plots Supplementary Figures: Pulse Duration

%% Colors
cols.off = [0 0 0];
cols.GAD2 = [0 0.4470 0.7410];
cols.GAD2_20ms = [0 0.2470 0.5410];
cols.GAD2_10ms = [0 0.0470 0.3410];
cols.Act = [0.8500 0.3250 0.0980];
cols.Inh = [0.9290 0.6940 0.1250];
% cols.Oth
cols.xline = [143 143 143]/255;
cols.xreg = [200 200 200]/255;
cols.xreg2 = [184 202 214]/255;
cols.err_bar = [26 35 126]/255;

%% Calculate Reliability

sAP = load("79155_20230512_AP.mat"); sAP = sAP.sAP;
GAD2_Ex1 = 353;
GAD2_Ex2 = 366;
% GAD2_cells = [353 366]
st1 = sAP.sCluster(GAD2_Ex1).SpikeTimes;
st2 = sAP.sCluster(GAD2_Ex2).SpikeTimes;
BinEdge = [0.001 0.005];
n_trials = (sAP.cellBlock{1,3}.intTrialNum)/numel(sAP.cellBlock{1,3}.vecPulseDuration);

stim_50ms = sAP.cellBlock{1, 3}.vecLaserOnTime(sAP.cellBlock{1, 3}.PulseDurTrial == 0.05);
sCounts_50ms = zeros(numel(stim_50ms),2);
for intTrial=1:n_trials
    vecTheseEdges = BinEdge + stim_50ms(intTrial); % Add stim onset time for this trial to rel. bin edges to get absolute bin edges
    [vecCounts1,edges1] = histcounts(st1,vecTheseEdges);
    [vecCounts2,edges2] = histcounts(st2,vecTheseEdges);
    sCounts_50ms(intTrial,1) = vecCounts1; % Counts Spontaneous Rate
    sCounts_50ms(intTrial,2) = vecCounts2; % Count for Visual Response
end

stim_20ms = sAP.cellBlock{1, 3}.vecLaserOnTime(sAP.cellBlock{1, 3}.PulseDurTrial == 0.02);
sCounts_20ms = zeros(numel(stim_20ms),2);
for intTrial=1:n_trials
    vecTheseEdges = BinEdge + stim_20ms(intTrial); % Add stim onset time for this trial to rel. bin edges to get absolute bin edges
    [vecCounts1,edges1] = histcounts(st1,vecTheseEdges);
    [vecCounts2,edges2] = histcounts(st2,vecTheseEdges);
    sCounts_20ms(intTrial,1) = vecCounts1; % Counts Spontaneous Rate
    sCounts_20ms(intTrial,2) = vecCounts2; % Count for Visual Response
end

stim_10ms = sAP.cellBlock{1, 3}.vecLaserOnTime(sAP.cellBlock{1, 3}.PulseDurTrial == 0.01);
sCounts_10ms = zeros(numel(stim_10ms),2);
for intTrial=1:n_trials
    vecTheseEdges = BinEdge + stim_10ms(intTrial); % Add stim onset time for this trial to rel. bin edges to get absolute bin edges
    [vecCounts1,edges1] = histcounts(st1,vecTheseEdges);
    [vecCounts2,edges2] = histcounts(st2,vecTheseEdges);
    sCounts_10ms(intTrial,1) = vecCounts1; % Counts Spontaneous Rate
    sCounts_10ms(intTrial,2) = vecCounts2; % Count for Visual Response
end

Reliability.f50ms = [sum(sCounts_50ms(:,1) ~= 0)/100 sum(sCounts_50ms(:,2) ~= 0)/100];
Reliability.f20ms = [sum(sCounts_20ms(:,1) ~= 0)/100 sum(sCounts_20ms(:,2) ~= 0)/100];
Reliability.f10ms = [sum(sCounts_10ms(:,1) ~= 0)/100 sum(sCounts_10ms(:,2) ~= 0)/100];

%% --- Rasters Example Neuron ---

figure;

% 50 ms
pulseSelect = 0.05;
stim = sAP.cellBlock{1, 3}.vecLaserOnTime(sAP.cellBlock{1, 3}.PulseDurTrial == pulseSelect) - 0.05;
maxDur = 0.15;

subplot(2,3,1);
plotRaster2(st1,stim,maxDur, [], cols.GAD2);
xregion(0.05, pulseSelect+0.05, 'FaceColor', cols.xreg) % xline(0.060, 'g--'); xline(0.080, 'g--');
xticks([0 0.05 0.1 0.15]);
xticklabels(["-0.05" "0" "0.05" "0.1"]);
xlabel('Time from laser onset (s)');
title("50 ms Pulse");
fixfig; hold off;

% subplot(2,3,2);
% plotRaster2(st2,stim,maxDur, [], cols.GAD2);
% xregion(0.05, pulseSelect+0.05, 'FaceColor', cols.xreg) % xline(0.060, 'g--'); xline(0.080, 'g--');
% xticks([0 0.05 0.1 0.15]);
% xticklabels(["-0.05" "0" "0.05" "0.1"]);
% xlabel('Time from laser onset (s)');
% fixfig; hold off;

% 20 ms
pulseSelect = 0.02
stim = sAP.cellBlock{1, 3}.vecLaserOnTime(sAP.cellBlock{1, 3}.PulseDurTrial == pulseSelect) - 0.05;
maxDur = 0.15;

subplot(2,3,2);
plotRaster2(st1,stim,maxDur, [], cols.GAD2);
xregion(0.05, pulseSelect+0.05, 'FaceColor', cols.xreg) % xline(0.060, 'g--'); xline(0.080, 'g--');
xticks([0 0.05 0.1 0.15]);
xticklabels(["-0.05" "0" "0.05" "0.1"]);
xlabel('Time from laser onset (s)');
title("20 ms Pulse");
fixfig; hold off;

% subplot(4,5,7);
% plotRaster2(st2,stim,maxDur, [], cols.GAD2);
% xregion(0.05, pulseSelect+0.05, 'FaceColor', cols.xreg) % xline(0.060, 'g--'); xline(0.080, 'g--');
% xticks([0 0.05 0.1 0.15]);
% xticklabels(["-0.05" "0" "0.05" "0.1"]);
% xlabel('Time from laser onset (s)');
% fixfig; hold off;

% 10 ms
pulseSelect = 0.01;
stim = sAP.cellBlock{1, 3}.vecLaserOnTime(sAP.cellBlock{1, 3}.PulseDurTrial == pulseSelect) - 0.05;
maxDur = 0.15;

subplot(2,3,3);
plotRaster2(st1,stim,maxDur, [], cols.GAD2);
xregion(0.05, pulseSelect+0.05, 'FaceColor', cols.xreg) % xline(0.060, 'g--'); xline(0.080, 'g--');
xticks([0 0.05 0.1 0.15]);
xticklabels(["-0.05" "0" "0.05" "0.1"]);
xlabel('Time from laser onset (s)');
title("10 ms Pulse");
fixfig; hold off;

% subplot(4,5,8);
% plotRaster2(st2,stim,maxDur, [], cols.GAD2);
% xregion(0.05, pulseSelect+0.05, 'FaceColor', cols.xreg) % xline(0.060, 'g--'); xline(0.080, 'g--');
% xticks([0 0.05 0.1 0.15]);
% xticklabels(["-0.05" "0" "0.05" "0.1"]);
% xlabel('Time from laser onset (s)');
% % title(string(sAP.sJson.subject) + " Cl: " + string(GAD2_Ex) + " Pulse: " + pulseSelect);
% fixfig; hold off;

%% --- Response Dynamics ---

BinCenters = DataOut_OT_50ms.Overall.PSTHBinCenters(1,:);

% 50 ms
pulseSelect = 0.05;
excl_idx = (BinCenters >= -0.001 & BinCenters <= 0.002) | (BinCenters >= (pulseSelect - 0.001) & (BinCenters <= pulseSelect + 0.002));
% PSTH (Baseline-Subtracted; Norm. to Peak)
PSTHMean_GAD2_norm_50ms = DataOut_OT_50ms.Overall.PSTHMean_GAD2_norm;
PSTHMean_GAD2_norm_50ms(excl_idx) = NaN;
PSTHSEM_GAD2_norm_50ms = DataOut_OT_50ms.Overall.PSTHSEM_GAD2_norm;
PSTHSEM_GAD2_norm_50ms(excl_idx) = NaN;

PSTHMean_Inh_norm_50ms = DataOut_OT_50ms.Overall.PSTHMean_Inh_norm;
PSTHMean_Inh_norm_50ms(excl_idx) = NaN;
PSTHSEM_Inh_norm_50ms = DataOut_OT_50ms.Overall.PSTHSEM_Inh_norm;
PSTHSEM_Inh_norm_50ms(excl_idx) = NaN;

PSTHMean_Act_norm_50ms = DataOut_OT_50ms.Overall.PSTHMean_Act_norm;
PSTHMean_Act_norm_50ms(excl_idx) = NaN;
PSTHSEM_Act_norm_50ms = DataOut_OT_50ms.Overall.PSTHSEM_Act_norm;
PSTHSEM_Act_norm_50ms(excl_idx) = NaN;

subplot(2,3,4); hold on;
shadedErrorBar(DataOut_OT_50ms.Overall.PSTHBinCenters, PSTHMean_GAD2_norm_50ms, PSTHSEM_GAD2_norm_50ms, 'lineProps', {'Color',cols.GAD2});
shadedErrorBar(DataOut_OT_50ms.Overall.PSTHBinCenters, PSTHMean_Inh_norm_50ms, PSTHSEM_Inh_norm_50ms, 'lineProps', {'Color',cols.Inh});
shadedErrorBar(DataOut_OT_50ms.Overall.PSTHBinCenters, PSTHMean_Act_norm_50ms, PSTHSEM_Act_norm_50ms, 'lineProps', {'Color',cols.Act});
xregion(0,pulseSelect,"FaceColor",cols.xreg,"LineStyle","none");
ylim([-0.2 1.2]);
xticks([-0.1 -0.05 0 0.05 0.1]);
xlim([-0.05 0.1]);
% legend(["GAD2+ (N = "+DataOut_OT_50ms.Overall.NCells_GAD2+")", "Inhibited (N = "+DataOut_OT_50ms.Overall.NCells_Inh+")", "Activated (N = "+DataOut_OT_50ms.Overall.NCells_Act+")"]);
legend(["GAD2+", "Inhibited", "Activated"]);
ylabel('Normalized Response');
xlabel('Time from laser onset (ms)');
% title('BL-Subtracted and Normalized to Peak');
fixfig; hold off;

% 20 ms

pulseSelect = 0.02;
excl_idx = (BinCenters >= -0.001 & BinCenters <= 0.002) | (BinCenters >= (pulseSelect - 0.001) & BinCenters <= (pulseSelect + 0.002));
% PSTH (Baseline-Subtracted; Norm. to Peak)
PSTHMean_GAD2_norm_20ms = DataOut_OT_20ms.Overall.PSTHMean_GAD2_norm;
PSTHMean_GAD2_norm_20ms(excl_idx) = NaN;
PSTHSEM_GAD2_norm_20ms = DataOut_OT_20ms.Overall.PSTHSEM_GAD2_norm;
PSTHSEM_GAD2_norm_20ms(excl_idx) = NaN;

PSTHMean_Inh_norm_20ms = DataOut_OT_20ms.Overall.PSTHMean_Inh_norm;
PSTHMean_Inh_norm_20ms(excl_idx) = NaN;
PSTHSEM_Inh_norm_20ms = DataOut_OT_20ms.Overall.PSTHSEM_Inh_norm;
PSTHSEM_Inh_norm_20ms(excl_idx) = NaN;

PSTHMean_Act_norm_20ms = DataOut_OT_20ms.Overall.PSTHMean_Act_norm;
PSTHMean_Act_norm_20ms(excl_idx) = NaN;
PSTHSEM_Act_norm_20ms = DataOut_OT_20ms.Overall.PSTHSEM_Act_norm;
PSTHSEM_Act_norm_20ms(excl_idx) = NaN;

subplot(2,3,5); hold on;
shadedErrorBar(DataOut_OT_20ms.Overall.PSTHBinCenters, PSTHMean_GAD2_norm_20ms, PSTHSEM_GAD2_norm_20ms, 'lineProps', {'Color',cols.GAD2});
shadedErrorBar(DataOut_OT_20ms.Overall.PSTHBinCenters, PSTHMean_Inh_norm_20ms, PSTHSEM_Inh_norm_20ms, 'lineProps', {'Color',cols.Inh});
shadedErrorBar(DataOut_OT_20ms.Overall.PSTHBinCenters, PSTHMean_Act_norm_20ms, PSTHSEM_Act_norm_20ms, 'lineProps', {'Color',cols.Act});
xregion(0,pulseSelect,"FaceColor",cols.xreg,"LineStyle","none");
ylim([-0.2 1.2]);
xticks([-0.1 -0.05 0 0.05 0.1]);
xlim([-0.05 0.1]);
% legend(["GAD2+ (N = "+DataOut_OT_20ms.Overall.NCells_GAD2+")", "Inhibited (N = "+DataOut_OT_20ms.Overall.NCells_Inh+")", "Activated (N = "+DataOut_OT_20ms.Overall.NCells_Act+")"]);
legend(["GAD2+", "Inhibited", "Activated"]);
ylabel('Normalized Response');
xlabel('Time from laser onset (ms)');
% title('BL-Subtracted and Normalized to Peak');
fixfig; hold off;

% 10 ms
pulseSelect = 0.01;
excl_idx = (BinCenters >= -0.001 & BinCenters <= 0.002) | (BinCenters >= (pulseSelect - 0.001) & BinCenters <= (pulseSelect + 0.002));
% PSTH (Baseline-Subtracted; Norm. to Peak)
PSTHMean_GAD2_norm_10ms = DataOut_OT_10ms.Overall.PSTHMean_GAD2_norm;
PSTHMean_GAD2_norm_10ms(excl_idx) = NaN;
PSTHSEM_GAD2_norm_10ms = DataOut_OT_10ms.Overall.PSTHSEM_GAD2_norm;
PSTHSEM_GAD2_norm_10ms(excl_idx) = NaN;

PSTHMean_Inh_norm_10ms = DataOut_OT_10ms.Overall.PSTHMean_Inh_norm;
PSTHMean_Inh_norm_10ms(excl_idx) = NaN;
PSTHSEM_Inh_norm_10ms = DataOut_OT_10ms.Overall.PSTHSEM_Inh_norm;
PSTHSEM_Inh_norm_10ms(excl_idx) = NaN;

PSTHMean_Act_norm_10ms = DataOut_OT_10ms.Overall.PSTHMean_Act_norm;
PSTHMean_Act_norm_10ms(excl_idx) = NaN;
PSTHSEM_Act_norm_10ms = DataOut_OT_10ms.Overall.PSTHSEM_Act_norm;
PSTHSEM_Act_norm_10ms(excl_idx) = NaN;

subplot(2,3,6); hold on;
shadedErrorBar(DataOut_OT_10ms.Overall.PSTHBinCenters, PSTHMean_GAD2_norm_10ms, PSTHSEM_GAD2_norm_10ms, 'lineProps', {'Color',cols.GAD2});
shadedErrorBar(DataOut_OT_10ms.Overall.PSTHBinCenters, PSTHMean_Inh_norm_10ms, PSTHSEM_Inh_norm_10ms, 'lineProps', {'Color',cols.Inh});
% shadedErrorBar(DataOut_OT_10ms.Overall.PSTHBinCenters, PSTHMean_Act_norm_10ms, PSTHSEM_Act_norm_10ms, 'lineProps', {'Color',cols.Act});
plot(DataOut_OT_10ms.Overall.PSTHBinCenters, PSTHMean_Act_norm_10ms, 'Color', cols.Act);
xregion(0,pulseSelect,"FaceColor",cols.xreg,"LineStyle","none");
ylim([-0.2 1.2]);
xticks([-0.1 -0.05 0 0.05 0.1]);
xlim([-0.05 0.1]);
% legend(["GAD2+ (N = "+DataOut_OT_10ms.Overall.NCells_GAD2+")", "Inhibited (N = "+DataOut_OT_10ms.Overall.NCells_Inh+")", "Activated (N = "+DataOut_OT_10ms.Overall.NCells_Act+")"]);
legend(["GAD2+", "Inhibited", "Activated"]);
ylabel('Normalized Response');
xlabel('Time from laser onset (ms)');
% title('BL-Subtracted and Normalized to Peak');
fixfig; hold off;

%% Compare Response
% subplot(3,3,[1 2 3]); hold on;
% shadedErrorBar(DataOut_OT_50ms.Overall.PSTHBinCenters, PSTHMean_GAD2_norm_50ms, PSTHSEM_GAD2_norm_50ms, 'lineProps', {'-', 'Color',cols.GAD2});
% shadedErrorBar(DataOut_OT_20ms.Overall.PSTHBinCenters, PSTHMean_GAD2_norm_20ms, PSTHSEM_GAD2_norm_20ms, 'lineProps', {'Color',cols.GAD2_20ms});
% shadedErrorBar(DataOut_OT_10ms.Overall.PSTHBinCenters, PSTHMean_GAD2_norm_10ms, PSTHSEM_GAD2_norm_10ms, 'lineProps', {'Color',cols.GAD2_10ms});
% xticks([-0.1 -0.05 0 0.05 0.1]);
% xlim([-0.05 0.1]);
% legend(["50 ms", "20 ms", "10 ms"]);
% ylabel('Normalized Response');
% xlabel('Time from laser onset (ms)');
% % title('BL-Subtracted and Normalized to Peak');
% fixfig; hold off;

%% Save Figure

% saveas(gcf, 'C:\Software and Code\Analysis-GAD2\Plots\S_Fig1.png');
% savefig(gcf, 'C:\Software and Code\Analysis-GAD2\Plots\S_Fig1.fig');
saveas(gcf, 'D:\NIN\Analysis-GAD2\Plots\S_Fig2.png');
savefig(gcf, 'D:\NIN\Analysis-GAD2\Plots\S_Fig2.fig');
