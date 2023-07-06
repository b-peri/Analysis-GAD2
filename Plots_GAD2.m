%% Plots Opto Gratings
Data = DataOut.AllMice;

% Barplots
figure; hold on;
bar_pl = bar([Data.Overall.ER_OptoOff Data.Overall.ER_OptoOn]);
dot_pl = plot([1 2],[Data.ClusterData.ER_OptoOff Data.ClusterData.ER_OptoOn],'-','Color',[0, 0, 0, 0.3]);
scatter([1 2],[Data.ClusterData.ER_OptoOff Data.ClusterData.ER_OptoOn],'MarkerEdgeColor', [0 0 0], 'MarkerEdgeAlpha',0.5);
bar_err = errorbar([1,2],[Data.Overall.ER_OptoOff Data.Overall.ER_OptoOn], [Data.Overall.SE_OptoOff Data.Overall.SE_OptoOn]);
bar_err.Color = [0.8 0 0];
bar_err.LineStyle = 'none';
bar_err.LineWidth = 1.5;
ylabel('Evoked Rate (spks/s)');
xticks([1 2]);
xticklabels({"Laser OFF", "Laser ON"});
hold off

% Mean Response
figure; hold on;
plot(Data.Overall.PSTHBinCenters, Data.Overall.PSTHMean_Off_norm,'k');
plot(Data.Overall.PSTHBinCenters, Data.Overall.PSTHMean_On_norm,'b');
xline(0,'r--');
ylabel('Spiking Rate (% of Peak)');
xlabel('Time from stimulus onset (s)');
xlim([min(Data.Overall.PSTHtime) max(Data.Overall.PSTHtime)]);
text(0.94, 0.62, sprintf('%g units', Data.Overall.NCells), 'FontSize',15);
legend('No Opto', 'Opto');
fixfig;
drawnow;
hold off;

%% Plots Receptive Field Mapping

% Location
% Preferred size!

