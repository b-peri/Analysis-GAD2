% Plots for Optotagging 
DOT = DataOut_OT_50ms.Overall;

GadTab_50ms = DataOut_OT_50ms.ClusterData(DataOut_OT_50ms.ClusterData.CellClass == "GAD2+",:);
GadTab_20ms = DataOut_OT_20ms.ClusterData(DataOut_OT_20ms.ClusterData.CellClass == "GAD2+",:);
GadTab_10ms = DataOut_OT_10ms.ClusterData(DataOut_OT_10ms.ClusterData.CellClass == "GAD2+",:);

% ActTab = DataOut_OT_50ms.ClusterData(DataOut_OT_50ms.ClusterData.CellClass == "Activated",:);
% InhTab = DataOut_OT_50ms.ClusterData(DataOut_OT_50ms.ClusterData.CellClass == "Inhibited",:);

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

%% Test Fig 3 A-C: Individual Plots - Raster & Instantaneous FR Plot



% -- A. GAD2 Cell Raster & Instantaneous FR (50 ms) --
% RasterPlot (79159_20230426: 446, 473)
pulseSelect = 0.05;
cell_sel_gad_50ms = [(GadTab_50ms.Subject + "_" + GadTab_50ms.RecDate + "_AP.mat") GadTab_50ms.ClusterN];

figure; 
for i = 1:length(cell_sel_gad_50ms)
    sAP = load(cell_sel_gad_50ms(i,1)); sAP = sAP.sAP;
    GAD2_Ex = str2double(cell_sel_gad_50ms(i,2));
    
    st = sAP.sCluster(GAD2_Ex).SpikeTimes;
    stim = sAP.cellBlock{1, 3}.vecLaserOnTime(sAP.cellBlock{1, 3}.PulseDurTrial == pulseSelect) - 0.05;
    maxDur = 0.15;
    
    % RasterPlot
    subplot(3,length(cell_sel_gad_50ms),i); hold on;
    plotRaster(st,stim,maxDur); xline(0.1, 'r--'); xline(0.101, 'b--');
    xline(0.110, 'g--'); xline(0.130, 'g--');
    title(GadTab_50ms(i,:).Subject + "; Cl: " + cell_sel_gad_50ms(i,2) + " Pulse: " + pulseSelect);
    fixfig; hold off;
    % IFR
    subplot(3,length(cell_sel_gad_50ms),i+length(cell_sel_gad_50ms)); hold on;
    plot(GadTab_50ms(i,:).IFR_Time{:}, GadTab_50ms(i,:).IFR_Rate{:});
    fixfig; hold off;
    % Average Waveform
    % subplot(3,length(cell_sel_gad_50ms),i+6); hold on;
    % plot(sAP.sCluster(GAD2_Ex).Waveform);
    % fixfig; hold off;
    subplot(3,length(cell_sel_gad_50ms),i+(length(cell_sel_gad_50ms)*2)); hold on;
    histogram((sAP.sCluster(GAD2_Ex).SpikeTimes),100)
    xline(sAP.cellBlock{2}.vecStimOnTime(1),'r-')
    xline(sAP.cellBlock{2}.vecStimOffTime(end),'r-')
end

% -- B. GAD2 Cell Raster & Instantaneous FR (20 ms) --
% RasterPlot (79159_20230426: 446, 473)
pulseSelect = 0.02;
cell_sel_gad_20ms = [(GadTab_20ms.Subject + "_" + GadTab_20ms.RecDate + "_AP.mat") GadTab_20ms.ClusterN];

figure; 
for i = 1:length(cell_sel_gad_20ms)
    sAP = load(cell_sel_gad_20ms(i,1)); sAP = sAP.sAP;
    GAD2_Ex = str2double(cell_sel_gad_20ms(i,2));
    
    st = sAP.sCluster(GAD2_Ex).SpikeTimes;
    stim = sAP.cellBlock{1, 3}.vecLaserOnTime(sAP.cellBlock{1, 3}.PulseDurTrial == pulseSelect) - 0.05;
    maxDur = 0.15;
    
    % RasterPlot
    subplot(3,length(cell_sel_gad_20ms),i); hold on;
    plotRaster(st,stim,maxDur); xline(0.05, 'r--');
    xline(0.06, 'g--'); xline(0.08, 'g--');
    title(GadTab_20ms(i,:).Subject + "; Cl: " + cell_sel_gad_20ms(i,2) + " Pulse: " + pulseSelect);
    fixfig; hold off;
    % IFR
    subplot(3,length(cell_sel_gad_20ms),i+length(cell_sel_gad_20ms)); hold on;
    plot(GadTab_20ms(i,:).IFR_Time{:}, GadTab_20ms(i,:).IFR_Rate{:});
    fixfig; hold off;
    % Average Waveform
    % subplot(3,length(cell_sel_gad_20ms),i+6); hold on;
    % plot(sAP.sCluster(GAD2_Ex).Waveform);
    % fixfig; hold off;
    subplot(3,length(cell_sel_gad_20ms),i+(length(cell_sel_gad_20ms)*2)); hold on;
    histogram((sAP.sCluster(GAD2_Ex).SpikeTimes),100)
    xline(sAP.cellBlock{2}.vecStimOnTime(1),'r-')
    xline(sAP.cellBlock{2}.vecStimOffTime(end),'r-')
end

% -- B. GAD2 Cell Raster & Instantaneous FR (10 ms) --
% RasterPlot (79159_20230426: 446, 473)
pulseSelect = 0.01;
cell_sel_gad_10ms = [(GadTab_10ms.Subject + "_" + GadTab_10ms.RecDate + "_AP.mat") GadTab_10ms.ClusterN];

figure; 
for i = 1:length(cell_sel_gad_10ms)
    sAP = load(cell_sel_gad_10ms(i,1)); sAP = sAP.sAP;
    GAD2_Ex = str2double(cell_sel_gad_10ms(i,2));
    
    st = sAP.sCluster(GAD2_Ex).SpikeTimes;
    stim = sAP.cellBlock{1, 3}.vecLaserOnTime(sAP.cellBlock{1, 3}.PulseDurTrial == pulseSelect) - 0.05;
    maxDur = 0.15;
    
    % RasterPlot
    subplot(3,length(cell_sel_gad_10ms),i); hold on;
    plotRaster(st,stim,maxDur); xline(0.05, 'r--');
    xline(0.06, 'g--'); xline(0.08, 'g--');
    title(GadTab_10ms(i,:).Subject + "; Cl: " + cell_sel_gad_10ms(i,2) + " Pulse: " + pulseSelect);
    fixfig; hold off;
    % IFR
    subplot(3,length(cell_sel_gad_10ms),i+length(cell_sel_gad_10ms)); hold on;
    plot(GadTab_10ms(i,:).IFR_Time{:}, GadTab_10ms(i,:).IFR_Rate{:});
    fixfig; hold off;
    % Average Waveform
    % subplot(3,length(cell_sel_gad_10ms),i+6); hold on;
    % plot(sAP.sCluster(GAD2_Ex).Waveform);
    % fixfig; hold off;
    subplot(3,length(cell_sel_gad_10ms),i+(length(cell_sel_gad_10ms)*2)); hold on;
    histogram((sAP.sCluster(GAD2_Ex).SpikeTimes),100)
    xline(sAP.cellBlock{2}.vecStimOnTime(1),'r-')
    xline(sAP.cellBlock{2}.vecStimOffTime(end),'r-')
end

%% 

% B. Activated/Disinhibited Cell
cell_sel_Act = [(ActTab.Subject + "_" + ActTab.RecDate + "_AP.mat") ActTab.ClusterN];

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

% C. Inhibited Cell

cell_sel_Inh = [(InhTab.Subject + "_" + InhTab.RecDate + "_AP.mat") InhTab.ClusterN];

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
    plot(InhTab(i,:).PSTHBinCenters, InhTab(i,:).PSTHMean);
    fixfig; hold off;
    % Histogram Over Whole Rec.
    subplot(3,1,3); hold on;
    histogram((sAP.sCluster(Inh_Ex).SpikeTimes),100)
    xline(sAP.cellBlock{2}.vecStimOnTime(1),'r-')
    xline(sAP.cellBlock{2}.vecStimOffTime(end),'r-')
end

% D. Other Cell 

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
xregion(0.05, pulseSelect+0.05, 'FaceColor', cols.xreg)
% title(string(sAP.sJson.subject) + " Cl: " + string(Act_Ex) + " Pulse: " + pulseSelect);
title("Activated");
fixfig; hold off;

%% Fig 3 D-E: Overall Plots

% Heatmap (-0.05s - 0.1s); Organized by 1. GAD2, 2. Inhibited, 3.
figure; BinCenters = DOT.PSTHBinCenters(1,:);
excl_idx = (BinCenters > -0.001 & BinCenters < 0.002) | (BinCenters > 0.049 & BinCenters < 0.052);

ReorderedTab = [GadTab; InhTab; ActTab];
PSTH_BL = ReorderedTab.PSTHMean - ReorderedTab.SpontRate;
PSTH_norm = PSTH_BL./max(PSTH_BL,[],2);
PSTH_norm(:,excl_idx) = NaN;
colormap('bone');
colorbar;
subplot(2,3,[2 3]); imagesc(PSTH_norm);
ylabel('Neuron #');
annotation('rectangle',[0.4 0.88 0.0075 0.045], 'FaceColor', cols.GAD2, 'LineStyle', 'none')
annotation('rectangle',[0.4 0.675 0.0075 0.202], 'FaceColor', cols.Inh, 'LineStyle', 'none')
annotation('rectangle',[0.4 0.585 0.0075 0.088], 'FaceColor', cols.Act, 'LineStyle', 'none')


% GadTab_PSTH_z = ((GadTab.PSTHMean - GadTab.SpontRate)./GadTab.SpontRate_STD);
% ActTab_PSTH_z = ((ActTab.PSTHMean - ActTab.SpontRate)./ActTab.SpontRate_STD);
% InhTab_PSTH_z = ((InhTab.PSTHMean - InhTab.SpontRate)./InhTab.SpontRate_STD);
% ReorderedPSTH = [GadTab_PSTH_z; InhTab_PSTH_z; ActTab_PSTH_z];

% 
% figure; hold on;
% plot(DOT.Overall.PSTHBinCenters, DOT.Overall.PSTHMean_GAD2);
% plot(DOT.Overall.PSTHBinCenters, DOT.Overall.PSTHMean_Act);
% plot(DOT.Overall.PSTHBinCenters, DOT.Overall.PSTHMean_Inh);
% legend(["GAD2+","Activated","Inhibited"]);
% ylabel('Mean Response (Z-Score)');
% xlabel('Time from laser onset (ms)');
% title('z-scored');
% hold off;

% % PSTH Gad2, Inhibited, Activated response dynamics: Z-Scored
% subplot(2,3,5); hold on;
% 
% PSTHMean_GAD2_z = DOT.PSTHMean_GAD2_z;
% PSTHMean_GAD2_z(excl_idx) = NaN;
% PSTHSEM_GAD2_z = DOT.PSTHSEM_GAD2_z;
% PSTHSEM_GAD2_z(excl_idx) = NaN;
% 
% PSTHMean_Act_z = DOT.PSTHMean_Act_z;
% PSTHMean_Act_z(excl_idx) = NaN;
% PSTHSEM_Act_z = DOT.PSTHSEM_Act_z;
% PSTHSEM_Act_z(excl_idx) = NaN;
% 
% PSTHMean_Inh_z = DOT.PSTHMean_Inh_z;
% PSTHMean_Inh_z(excl_idx) = NaN;
% PSTHSEM_Inh_z = DOT.PSTHSEM_Inh_z;
% PSTHSEM_Inh_z(excl_idx) = NaN;
% 
% shadedErrorBar(DOT.PSTHBinCenters, PSTHMean_GAD2_z, PSTHSEM_GAD2_z, 'lineProps', {'Color',cols.GAD2});
% shadedErrorBar(DOT.PSTHBinCenters, PSTHMean_Act_z, PSTHSEM_Act_z, 'lineProps', {'Color',cols.Act});
% shadedErrorBar(DOT.PSTHBinCenters, PSTHMean_Inh_z, PSTHSEM_Inh_z, 'lineProps', {'Color', cols.Inh});
% xlim([-0.1 0.1]);
% % ylim([-0.2 1.2])
% legend(["GAD2+","Activated","Inhibited"]);
% % legend(["Activated","Inhibited"]);
% ylabel('Mean Response (Z-Score)');
% xlabel('Time from laser onset (ms)');
% % title('Z-Scored');
% fixfig; hold off;

% Baseline-Corrected
subplot(2,3,[5 6]); hold on;
shadedErrorBar(DOT.PSTHBinCenters, DOT.PSTHMean_GAD2_norm, DOT.PSTHSEM_GAD2_norm, 'lineProps', {'Color',cols.GAD2});
shadedErrorBar(DOT.PSTHBinCenters, DOT.PSTHMean_Act_norm, DOT.PSTHSEM_Act_norm, 'lineProps', {'Color',cols.Act});
shadedErrorBar(DOT.PSTHBinCenters, DOT.PSTHMean_Inh_norm, DOT.PSTHSEM_Inh_norm, 'lineProps', {'Color',cols.Inh});
xlim([-0.1 0.1]);
ylim([-0.2 1.2])
legend(["GAD2+","Activated","Inhibited"]);
ylabel('Mean Response (Normalized to Peak)');
xlabel('Time from laser onset (ms)');
% title('BL-Subtracted and Normalized to Peak');
fixfig; hold off;

%% Boxplots Latency

% [] -> Ask Robin about whether boxplot is the most appropriate way to
% demonstrate this point given that N is so small

figure;hold on;

boxplot(DataOut_OT_50ms.ClusterData.Peak_Lat,DataOut_OT_50ms.ClusterData.CellClass);
xlabel("Peak Latency (s)");

%% Suppl. Figure 1

% Maybe some plots comparing different pulse durations?