% Plot Untagged

figure; hold on;
shadedErrorBar(test_Un_PSTH.PSTHBinCenters, test_Un_PSTH.PSTHMean_Off_norm, test_Un_PSTH.PSTHSEM_Off_norm, 'lineProps', 'k');
shadedErrorBar(test_Un_PSTH.PSTHBinCenters, test_Un_PSTH.PSTHMean_On_norm, test_Un_PSTH.PSTHSEM_On_norm, 'lineProps', 'b');
xline(0,'r--');
ylabel('Normalized Response');
xlabel('Time from stimulus onset (s)');
% xlim([min(test_Un_PSTH.PSTHtime) max(test_Un_PSTH.PSTHtime)]);
% ylim([-0.1 1]);
yline(0);
text(0.94, 0.62, sprintf('%g units', height(test_Un)), 'FontSize',15);
legend('No Opto', 'Opto');
title('Mean Response Untagged');
fixfig;
drawnow;
hold off;

% Plot Activated

figure; hold on;
shadedErrorBar(test_Act_PSTH.PSTHBinCenters, test_Act_PSTH.PSTHMean_Off_norm, test_Act_PSTH.PSTHSEM_Off_norm, 'lineProps', 'k');
shadedErrorBar(test_Act_PSTH.PSTHBinCenters, test_Act_PSTH.PSTHMean_On_norm, test_Act_PSTH.PSTHSEM_On_norm, 'lineProps', 'b');
xline(0,'r--');
ylabel('Normalized Response');
xlabel('Time from stimulus onset (s)');
% xlim([min(test_Act_PSTH.PSTHtime) max(test_Act_PSTH.PSTHtime)]);
% ylim([-0.1 1]);
yline(0);
% text(0.94, 0.62, sprintf('%g units', height(test_Act), 'FontSize',15);
legend('No Opto', 'Opto');
title('Mean Response Activated');
fixfig;
drawnow;
hold off;

% Plot Inhibited

figure; hold on;
shadedErrorBar(test_Inh_PSTH.PSTHBinCenters, test_Inh_PSTH.PSTHMean_Off_norm, test_Inh_PSTH.PSTHSEM_Off_norm, 'lineProps', 'k');
shadedErrorBar(test_Inh_PSTH.PSTHBinCenters, test_Inh_PSTH.PSTHMean_On_norm, test_Inh_PSTH.PSTHSEM_On_norm, 'lineProps', 'b');
xline(0,'r--');
ylabel('Normalized Response');
xlabel('Time from stimulus onset (s)');
% xlim([min(test_Inh_PSTH.PSTHtime) max(test_Inh_PSTH.PSTHtime)]);
% ylim([-0.1 1]);
yline(0);
% text(0.94, 0.62, sprintf('%g units', height(test_Inh), 'FontSize',15);
legend('No Opto', 'Opto');
title('Mean Response Inhibited');
fixfig;
drawnow;
hold off;