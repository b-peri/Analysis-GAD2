% Plots for Optotagging 
DOT = DataOut_OT.AllMice;

%% Fig 3 A-C: Individual Plots - Raster & Instantaneous FR Plot

% A. GAD2 Cell Raster & Instantaneous FR (PSTHish)

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
plot(DOT.Overall.PSTHBinCenters, DOT.Overall.PSTHMean_20ms_GAD2_norm);
plot(DOT.Overall.PSTHBinCenters, DOT.Overall.PSTHMean_20ms_Act_norm);
plot(DOT.Overall.PSTHBinCenters, DOT.Overall.PSTHMean_20ms_Inh_norm);
legend(["GAD2+","Activated","Inhibited"]);
ylabel('Mean Response (Normalized to Peak)');
xlabel('Time from laser onset (ms)');
title('BL-Subtracted and Normalized to Peak');
hold off;

%% Suppl. Figure 1

% Maybe some plots comparing different pulse durations?