%% Figure 1: Plots Opto Gratings
Data = DataOut.AllMice;

%% Example Neuron Plots
NumPlots_Down = 3 % Number of Downregulated Example Cells to Plot
NumPlots_Up = 1 % Number of Upregulated Example Cells to Plot

ClIdx_Down = randi(height(Data.ClusterDataDown), [1 NumPlots_Down]);
figure; hold on;
for n = 1:NumPlots_Down % Should probably normalize 

    i = ClIdx_Down(n)
    data_tmp = Data.ClusterDataDown(i,:);
    plot(data_tmp.PSTHBinCenters_Off, data_tmp.PSTHMean_Off, 'k');
    plot(data_tmp.PSTHBinCenters_On, data_tmp.PSTHMean_On, 'b');
    title(["Cl: " + string(data_tmp.ClusterN)]);
    xline(0,'r--');
    ylabel('Normalized Response');
    xlabel('Time from stimulus onset (s)');
    xlim([min(Data.Overall.PSTHtime) max(Data.Overall.PSTHtime)]);
    legend('No Opto', 'Opto');
end
hold off

ClIdx_Up = randi(height(Data.ClusterDataUp), [1 NumPlots_Up]);
figure; hold on
for i = ClIdx_Up % Should probably normalize 
    data_tmp = Data.ClusterDataUp(i,:);
    figure; hold on;
    plot(data_tmp.PSTHBinCenters_Off, data_tmp.PSTHMean_Off, 'k');
    plot(data_tmp.PSTHBinCenters_On, data_tmp.PSTHMean_On, 'b');
    title(["Cl: " + string(data_tmp.ClusterN)]);
    xline(0,'r--');
    ylabel('Normalized Response');
    xlabel('Time from stimulus onset (s)');
    xlim([min(Data.Overall.PSTHtime) max(Data.Overall.PSTHtime)]);
    legend('No Opto', 'Opto');
end
hold off

%% PSTH

% Mean Response Overall
figure; hold on;
% plot(Data.Overall.PSTHBinCenters, Data.Overall.PSTHMean_Off_norm,'k');
% plot(Data.Overall.PSTHBinCenters, Data.Overall.PSTHMean_On_norm,'b');
shadedErrorBar(Data.Overall.PSTHBinCenters, Data.Overall.PSTHMean_Off_norm, Data.Overall.PSTHSEM_Off_norm, 'lineProps', 'k');
shadedErrorBar(Data.Overall.PSTHBinCenters, Data.Overall.PSTHMean_On_norm, Data.Overall.PSTHSEM_On_norm, 'lineProps', 'b');
xline(0,'r--');
ylabel('Normalized Response');
xlabel('Time from stimulus onset (s)');
xlim([min(Data.Overall.PSTHtime) max(Data.Overall.PSTHtime)]);
ylim([-0.1 0.8]);
yline(0);
text(0.94, 0.62, sprintf('%g units', Data.Overall.NCells), 'FontSize',15);
legend('No Opto', 'Opto');
title('Mean Response Overall');
fixfig;
drawnow;
hold off;

% Mean Response Upregulated Cells -> Normalized
norm_up = max(Data.ClusterDataUp.PSTHMean_Off,[],2);
figure; hold on;
% plot(Data.Overall.PSTHBinCenters, mean(Data.ClusterDataUp.PSTHMean_Off./norm_up, 1),'k');
% plot(Data.Overall.PSTHBinCenters, mean(Data.ClusterDataUp.PSTHMean_On./norm_up, 1),'b');
shadedErrorBar(Data.Overall.PSTHBinCenters, mean(Data.ClusterDataUp.PSTHMean_Off./norm_up, 1), ...
    std(Data.ClusterDataUp.PSTHMean_Off./norm_up)/sqrt(Data.Overall.NCells_Increased), 'lineProps', 'k');
shadedErrorBar(Data.Overall.PSTHBinCenters, mean(Data.ClusterDataUp.PSTHMean_On./norm_up, 1), ..., 'lineProps', 'b');
    std(Data.ClusterDataUp.PSTHMean_On./norm_up)/sqrt(Data.Overall.NCells_Increased), 'lineProps', 'b');
xline(0,'r--');
ylabel('Normalized Response');
xlabel('Time from stimulus onset (s)');
xlim([min(Data.Overall.PSTHtime) max(Data.Overall.PSTHtime)]);
% ylim([-0.1 1]);
yline(0);
text(0.94, 2.2, sprintf('%g units', Data.Overall.NCells_Increased), 'FontSize',15);
legend('No Opto', 'Opto');
title('Upregulated Cells');
fixfig;
drawnow;
hold off;

% Mean Response Downregulated Cells
norm_down = max(Data.ClusterDataDown.PSTHMean_Off,[],2);
figure; hold on;
% plot(Data.Overall.PSTHBinCenters, mean(Data.ClusterDataDown.PSTHMean_Off./norm_down, 1),'k');
% plot(Data.Overall.PSTHBinCenters, mean(Data.ClusterDataDown.PSTHMean_On./norm_down, 1),'b');
shadedErrorBar(Data.Overall.PSTHBinCenters, mean(Data.ClusterDataDown.PSTHMean_Off./norm_down, 1), ...
    std(Data.ClusterDataDown.PSTHMean_Off./norm_down)/sqrt(Data.Overall.NCells_Reduced), 'lineProps', 'k');
shadedErrorBar(Data.Overall.PSTHBinCenters, mean(Data.ClusterDataDown.PSTHMean_On./norm_down, 1), ...
    std(Data.ClusterDataDown.PSTHMean_On./norm_down)/sqrt(Data.Overall.NCells_Reduced), 'lineProps', 'b');
xline(0,'r--');
ylabel('Normalized Response');
xlabel('Time from stimulus onset (s)');
xlim([min(Data.Overall.PSTHtime) max(Data.Overall.PSTHtime)]);
% ylim([-0.1 1]);
yline(0);
text(0.94, 0.62, sprintf('%g units', Data.Overall.NCells_Reduced), 'FontSize',15);
legend('No Opto', 'Opto');
title('Downregulated Cells');
fixfig;
drawnow;
hold off;

% Scatter Plot SR Opto vs SR No-Opto -> Split axis to include outliers
figure; hold on;
scatter(Data.ClusterDataUp.ER_OptoOff, Data.ClusterDataUp.ER_OptoOn, 'r', 'filled'); % Scatter SigUp
scatter(Data.ClusterDataDown.ER_OptoOff, Data.ClusterDataDown.ER_OptoOn, 'b', 'filled'); % Scatter SigDown
scatter(Data.ClusterDataNonSig.ER_OptoOff, Data.ClusterDataNonSig.ER_OptoOn, 'k'); % Scatter NonSig
plot([-20:1:100],[-20:1:100],'--'); % Reference Line
% ylim([-15 100]);
xlim([-15 100]);
ylabel('Evoked Rate Opto ON (spks/s)');
xlabel('Evoked Rate Opto OFF (spks/s)');
legend('Upregulated (p<0.01)', 'Downregulated (p<0.01)', 'No Effect');
grid on
axis square
% fixfig;
hold off

% Barplots
figure; hold on;
bar_pl = bar([Data.Overall.ER_OptoOff Data.Overall.ER_OptoOn]);
% dot_pl = plot([1 2],[Data.ClusterData.ER_OptoOff Data.ClusterData.ER_OptoOn],'-','Color',[0, 0, 0, 0.3]);
% scatter([1 2],[Data.ClusterData.ER_OptoOff Data.ClusterData.ER_OptoOn],'MarkerEdgeColor', [0 0 0], 'MarkerEdgeAlpha',0.5);
bar_err = errorbar([1,2],[Data.Overall.ER_OptoOff Data.Overall.ER_OptoOn], [Data.Overall.SE_OptoOff Data.Overall.SE_OptoOn]);
bar_err.Color = [0.8 0 0];
bar_err.LineStyle = 'none';
bar_err.LineWidth = 1.5;
ylabel('Evoked Rate (spks/s)');
xticks([1 2]);
xticklabels({"Laser OFF", "Laser ON"});
fixfig;
hold off

%% Supplementary Plots for Leonie


%% Plots Receptive Field Mapping

% Location
% Preferred size!

