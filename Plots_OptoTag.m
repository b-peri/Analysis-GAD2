% Plots for Optotagging 
DOT = DataOut_OT.Overall;

%% Fig 3 A-C: Individual Plots - Raster & Instantaneous FR Plot

% -- A. GAD2 Cell Raster & Instantaneous FR (PSTHish) --
% RasterPlot (79159_20230426: 446, 473)
GadTab = DataOut_OT.ClusterData(DataOut_OT.ClusterData.CellClass == "GAD2+",:);

sAP1 = load("79155_20230512_AP.mat"); sAP1 = sAP1.sAP;
GAD2_Ex = 366;
% GAD2_Ex_Tab = DataOut_OT.ClusterData(DataOut_OT.ClusterData.RecDate == "20230512" & DataOut_OT.ClusterData.ClusterN == 366, :);

st = sAP1.sCluster(GAD2_Ex).SpikeTimes;
stim = sAP1.cellBlock{1, 3}.vecLaserOnTime(sAP1.cellBlock{1, 3}.PulseDurTrial == 0.02) - 0.1;
maxDur = 0.2;

figure; hold on;
plotRaster(st,stim,maxDur); xline(0.1, 'r--'); xline(0.101, 'b--');
xline(0.110, 'g--'); xline(0.130, 'g--');
fixfig; hold off;
% IFR
figure;
plot(GAD2_Ex_Tab.IFR_Time{:}, GAD2_Ex_Tab.IFR_Rate{:});

% B. Inhibited Cell

% C. Activated/Disinhibited Cell

%% Fig 3 D-E: Overall Plots

% Heatmap (-0.05s - 0.1s); Organized by 1. GAD2, 2. Inhibited, 3.
% Activated/Disinhibited
PSTHMean_20ms_BL = DOT.ClusterData.PSTHMean_20ms - DOT.SpontRate;
norm = max(PSTHMean_20ms_BL,[],2);

ReorderedPSTH = []
imagesc

% PSTH Overlay Gad2, Inhibited, Activated response dynamics
figure; hold on;
plot(DOT.Overall.PSTHBinCenters, DOT.Overall.PSTHMean_20ms_GAD2);
plot(DOT.Overall.PSTHBinCenters, DOT.Overall.PSTHMean_20ms_Act);
plot(DOT.Overall.PSTHBinCenters, DOT.Overall.PSTHMean_20ms_Inh);
legend(["GAD2+","Activated","Inhibited"]);
ylabel('Mean Response (Z-Score)');
xlabel('Time from laser onset (ms)');
title('z-scored');
hold off;

figure; hold on;
plot(DOT.PSTHBinCenters, DOT.PSTHMean_20ms_GAD2_norm);
plot(DOT.PSTHBinCenters, DOT.PSTHMean_20ms_Act_norm);
plot(DOT.PSTHBinCenters, DOT.PSTHMean_20ms_Inh_norm);
legend(["GAD2+","Activated","Inhibited"]);
ylabel('Mean Response (Normalized to Peak)');
xlabel('Time from laser onset (ms)');
title('BL-Subtracted and Normalized to Peak');
hold off;

%% Suppl. Figure 1

% Maybe some plots comparing different pulse durations?