% Plots for Optotagging 
DOT = DataOut_OT.Overall;

GadTab = DataOut_OT.ClusterData(DataOut_OT.ClusterData.CellClass == "GAD2+",:);
ActTab = DataOut_OT.ClusterData(DataOut_OT.ClusterData.CellClass == "Activated",:);
InhTab = DataOut_OT.ClusterData(DataOut_OT.ClusterData.CellClass == "Inhibited",:);

%% Fig 3 A-C: Individual Plots - Raster & Instantaneous FR Plot

pulseSelect = 0.05

% -- A. GAD2 Cell Raster & Instantaneous FR (PSTHish) --
% RasterPlot (79159_20230426: 446, 473)
cell_sel_gad = [(GadTab.Subject + "_" + GadTab.RecDate + "_AP.mat") GadTab.ClusterN]

figure; 
for i = 1:length(cell_sel_gad)
    sAP = load(cell_sel_gad(i,1)); sAP = sAP.sAP;
    GAD2_Ex = str2double(cell_sel_gad(i,2));
    
    st = sAP.sCluster(GAD2_Ex).SpikeTimes;
    stim = sAP.cellBlock{1, 3}.vecLaserOnTime(sAP.cellBlock{1, 3}.PulseDurTrial == pulseSelect) - 0.1;
    maxDur = 0.2;
    
    % RasterPlot
    subplot(3,length(cell_sel_gad),i); hold on;
    plotRaster(st,stim,maxDur); xline(0.1, 'r--'); xline(0.101, 'b--');
    xline(0.110, 'g--'); xline(0.130, 'g--');
    title(GadTab(i,:).Subject + "; Cl: " + cell_sel_gad(i,2) + " Pulse: " + pulseSelect);
    fixfig; hold off;
    % IFR
    subplot(3,length(cell_sel_gad),i+length(cell_sel_gad)); hold on;
    plot(GadTab(i,:).IFR_Time{:}, GadTab(i,:).IFR_Rate{:});
    fixfig; hold off;
    % Average Waveform
    % subplot(3,length(cell_sel_gad),i+6); hold on;
    % plot(sAP.sCluster(GAD2_Ex).Waveform);
    % fixfig; hold off;
    subplot(3,length(cell_sel_gad),i+(length(cell_sel_gad)*2)); hold on;
    histogram((sAP.sCluster(GAD2_Ex).SpikeTimes),100)
    xline(sAP.cellBlock{2}.vecStimOnTime(1),'r-')
    xline(sAP.cellBlock{2}.vecStimOffTime(end),'r-')
end

% B. Activated/Disinhibited Cell
cell_sel_Act = [(ActTab.Subject + "_" + ActTab.RecDate + "_AP.mat") ActTab.ClusterN]

for i = 1:length(cell_sel_Act)
    sAP = load(cell_sel_Act(i,1)); sAP = sAP.sAP;
    Act_Ex = str2double(cell_sel_Act(i,2));
    
    st = sAP.sCluster(Act_Ex).SpikeTimes;
    stim = sAP.cellBlock{1, 3}.vecLaserOnTime(sAP.cellBlock{1, 3}.PulseDurTrial == pulseSelect) - 0.1;
    maxDur = 0.2;
    
    figure; subplot(3,1,1); hold on;
    plotRaster(st,stim,maxDur); xline(0.1, 'r--'); xline(0.101, 'b--');
    xline(0.110, 'g--'); xline(0.130, 'g--');
    title(ActTab(i,:).Subject + "; Cl: " + cell_sel_Act(i,2) + " Pulse: " + pulseSelect);
    fixfig; hold off;
    % IFR
    subplot(3,1,2); hold on;
    plot(ActTab(i,:).IFR_Time{:}, ActTab(i,:).IFR_Rate{:});
    fixfig; hold off;
    % Average Waveform
    subplot(3,1,3); hold on;
    histogram((sAP.sCluster(Act_Ex).SpikeTimes),100)
    xline(sAP.cellBlock{2}.vecStimOnTime(1),'r-')
    xline(sAP.cellBlock{2}.vecStimOffTime(end),'r-')
end

% B. Inhibited Cell

cell_sel_Inh = [(InhTab.Subject + "_" + InhTab.RecDate + "_AP.mat") InhTab.ClusterN]

for i = 1:length(cell_sel_Inh)
    sAP = load(cell_sel_Inh(i,1)); sAP = sAP.sAP;
    Inh_Ex = str2double(cell_sel_Inh(i,2));
    
    st = sAP.sCluster(Inh_Ex).SpikeTimes;
    stim = sAP.cellBlock{1, 3}.vecLaserOnTime(sAP.cellBlock{1, 3}.PulseDurTrial == pulseSelect) - 0.1;
    maxDur = 0.2;
    
    figure; subplot(3,1,1); hold on;
    plotRaster(st,stim,maxDur); xline(0.1, 'r--'); xline(0.101, 'b--');
    xline(0.110, 'g--'); xline(0.130, 'g--');
    title(InhTab(i,:).Subject + "; Cl: " + cell_sel_Inh(i,2) + " Pulse: " + pulseSelect);
    fixfig; hold off;
    % IFR
    subplot(3,1,2); hold on;
    plot(InhTab(i,:).IFR_Time{:}, InhTab(i,:).IFR_Rate{:});
    fixfig; hold off;
    % Histogram Over Whole Rec.
    subplot(3,1,3); hold on;
    histogram((sAP.sCluster(Inh_Ex).SpikeTimes),100)
    xline(sAP.cellBlock{2}.vecStimOnTime(1),'r-')
    xline(sAP.cellBlock{2}.vecStimOffTime(end),'r-')
end

%% Fig 3 D-E: Overall Plots

% Heatmap (-0.05s - 0.1s); Organized by 1. GAD2, 2. Inhibited, 3.

BinCenters = DOT.PSTHBinCenters(1,:);
excl_idx = (BinCenters > -0.001 & BinCenters < 0.006) | (BinCenters > 0.047 & BinCenters < 0.053);

GadTab_PSTH_z = ((GadTab.PSTHMean_20ms - GadTab.SpontRate)./GadTab.SpontRate_STD);
ActTab_PSTH_z = ((ActTab.PSTHMean_20ms - ActTab.SpontRate)./ActTab.SpontRate_STD);
InhTab_PSTH_z = ((InhTab.PSTHMean_20ms - InhTab.SpontRate)./InhTab.SpontRate_STD);

for i = 1:8
    figure; hold on
    bar(InhTab_PSTH_z(i,:))
    yline(-1,'r')
end

ReorderedPSTH = [GadTab_PSTH_z; InhTab_PSTH_z; ActTab_PSTH_z];
ReorderedPSTH(:,excl_idx) = NaN;
figure; imagesc(ReorderedPSTH); yline(2)

% PSTH Overlay Gad2, Inhibited, Activated response dynamics
% figure; hold on;
% plot(DOT.Overall.PSTHBinCenters, DOT.Overall.PSTHMean_20ms_GAD2);
% plot(DOT.Overall.PSTHBinCenters, DOT.Overall.PSTHMean_20ms_Act);
% plot(DOT.Overall.PSTHBinCenters, DOT.Overall.PSTHMean_20ms_Inh);
% legend(["GAD2+","Activated","Inhibited"]);
% ylabel('Mean Response (Z-Score)');
% xlabel('Time from laser onset (ms)');
% title('z-scored');
% hold off;

% Z-Scored
figure; hold on;

PSTHMean_GAD2_z = DOT.PSTHMean_20ms_GAD2_z;
PSTHMean_GAD2_z(excl_idx) = NaN;
PSTHSEM_GAD2_z = DOT.PSTHSEM_20ms_GAD2_z;
PSTHSEM_GAD2_z(excl_idx) = NaN;

PSTHMean_Act_z = DOT.PSTHMean_20ms_Act_z;
PSTHMean_Act_z(excl_idx) = NaN;
PSTHSEM_Act_z = DOT.PSTHSEM_20ms_Act_z;
PSTHSEM_Act_z(excl_idx) = NaN;

PSTHMean_Inh_z = DOT.PSTHMean_20ms_Inh_z;
PSTHMean_Inh_z(excl_idx) = NaN;
PSTHSEM_Inh_z = DOT.PSTHSEM_20ms_Inh_z;
PSTHSEM_Inh_z(excl_idx) = NaN;
shadedErrorBar(DOT.PSTHBinCenters, PSTHMean_GAD2_z, PSTHSEM_GAD2_z, 'lineProps', {'Color',[0 0.4470 0.7410]});
shadedErrorBar(DOT.PSTHBinCenters, PSTHMean_Act_z, PSTHSEM_Act_z, 'lineProps', {'Color',[0.8500 0.3250 0.0980]});
shadedErrorBar(DOT.PSTHBinCenters, PSTHMean_Inh_z, PSTHSEM_Inh_z, 'lineProps', {'Color',[0.9290 0.6940 0.1250]});


% shadedErrorBar(DOT.PSTHBinCenters, DOT.PSTHMean_20ms_GAD2_z, DOT.PSTHSEM_20ms_GAD2_z, 'lineProps', {'Color',[0 0.4470 0.7410]});
% shadedErrorBar(DOT.PSTHBinCenters, DOT.PSTHMean_20ms_Act_z, DOT.PSTHSEM_20ms_Act_z, 'lineProps', {'Color',[0.8500 0.3250 0.0980]});
% shadedErrorBar(DOT.PSTHBinCenters, DOT.PSTHMean_20ms_Inh_z, DOT.PSTHSEM_20ms_Inh_z, 'lineProps', {'Color',[0.9290 0.6940 0.1250]});
xlim([-0.1 0.1]);
% ylim([-0.2 1.2])
legend(["GAD2+","Activated","Inhibited"]);
% legend(["Activated","Inhibited"]);
ylabel('Mean Response (Z-Score)');
xlabel('Time from laser onset (ms)');
title('Z-Scored');
fixfig; hold off;

% Baseline-Corrected
figure; hold on;
plot(DOT.PSTHBinCenters, DOT.PSTHMean_20ms_GAD2_norm);
plot(DOT.PSTHBinCenters, DOT.PSTHMean_20ms_Act_norm);
plot(DOT.PSTHBinCenters, DOT.PSTHMean_20ms_Inh_norm);
xlim([-0.1 0.1]);
ylim([-0.2 1.2])
legend(["GAD2+","Activated","Inhibited"]);
ylabel('Mean Response (Normalized to Peak)');
xlabel('Time from laser onset (ms)');
title('BL-Subtracted and Normalized to Peak');
fixfig; hold off;

%% Suppl. Figure 1

% Maybe some plots comparing different pulse durations?