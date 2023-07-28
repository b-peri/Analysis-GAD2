% Plot Untagged

figure; hold on;
shadedErrorBar(Un_PSTH.PSTHBinCenters, Un_PSTH.PSTHMean_Off_z, Un_PSTH.PSTHSEM_Off_z, 'lineProps', 'k');
shadedErrorBar(Un_PSTH.PSTHBinCenters, Un_PSTH.PSTHMean_On_z, Un_PSTH.PSTHSEM_On_z, 'lineProps', 'b');
xline(0,'r--');
ylabel('Normalized Response (Z-Score)');
xlabel('Time from stimulus onset (s)');
xlim([-0.2 1.2]);
% ylim([-0.1 1]);
yline(0);
text(1, 2.25, sprintf('%g units', height(test_Un)), 'FontSize',15);
legend('No Opto', 'Opto');
title('Mean Response Untagged');
fixfig;
drawnow;
hold off;

% Plot Activated

figure; hold on;
shadedErrorBar(Act_PSTH.PSTHBinCenters, Act_PSTH.PSTHMean_Off_z, Act_PSTH.PSTHSEM_Off_z, 'lineProps', 'k');
shadedErrorBar(Act_PSTH.PSTHBinCenters, Act_PSTH.PSTHMean_On_z, Act_PSTH.PSTHSEM_On_z, 'lineProps', 'b');
xline(0,'r--');
ylabel('Normalized Response (Z-Score)');
xlabel('Time from stimulus onset (s)');
xlim([-0.2 1.2]);
% ylim([-0.1 1]);
yline(0);
text(1, 23, sprintf('%g units', height(test_Act)), 'FontSize',15);
legend('No Opto', 'Opto');
title('Mean Response Activated');
fixfig;
drawnow;
hold off;

% Plot Inhibited

figure; hold on;
shadedErrorBar(Inh_PSTH.PSTHBinCenters, Inh_PSTH.PSTHMean_Off_z, Inh_PSTH.PSTHSEM_Off_z, 'lineProps', 'k');
shadedErrorBar(Inh_PSTH.PSTHBinCenters, Inh_PSTH.PSTHMean_On_z, Inh_PSTH.PSTHSEM_On_z, 'lineProps', 'b');
xline(0,'r--');
ylabel('Normalized Response (Z-Score)');
xlabel('Time from stimulus onset (s)');
xlim([-0.2 1.2]);
% ylim([-0.1 1]);
yline(0);
text(1, 9, sprintf('%g units', height(test_Inh)), 'FontSize',15);
legend('No Opto', 'Opto');
title('Mean Response Inhibited');
fixfig;
drawnow;
hold off;