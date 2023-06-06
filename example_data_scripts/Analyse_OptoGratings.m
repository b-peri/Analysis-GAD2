%% OptoGratings Analysis
% 
% [] Why is vecStimOffSecs 1 second BEFORE stim onsets?
% 
% Robin's example script makes PLOTS, but what we want is an OUTPUT FILE
% that should give:
% [] Per Cluster: Optograting PSTH
% [] Per Mouse: Some metric of OVERALL reduction in sSC specifically (in firing rate?)
% [] OPTIONAL: Conversion to CSV

% [] Figure out where/whether it says somewhere if clusters are MU or SU
% [] Separate Figs for SINGLE-UNIT and MULTI-UNIT data?


%% 
% Get stimulus info
sAP = sSynthData;
intNumClu = length(sAP.sCluster); % Grabs number of clusters/putative cells (Idk where I can see whether this is single or MU?)
structEP = sAP.cellStim{1,1}.structEP;  % structEP contains data on a single recording block (incl. stim onset/offset times)
% -> This index is hardcoded. Not a huge deal if OptoGratings is always the
% first recording, but not very pretty/generalizable
vecStimOnSecs = structEP.vecStimOnTime; % For OptoGratings, synchronized stimulus onset times
vecStimOffSecs = structEP.vecStimOffTime; % For Optogratings, synchronized stimulus OFFset times
vecLaserOn = structEP.vecOptoOn; % Logical array: Tells whether opto was on for each trial!

%% Prepare Output Table
% 

%% loop through clusters
sOptions.handleFig = -1; % COME BACK TO THIS
dblBinDur = 5e-3; % Binsize of PSTH
vecTime = -0.2:dblBinDur:1.2; % PSTH X-Axis range
%-0.2:dblBinDur:0.5; %
indExcludeOn = vecTime > -2* dblBinDur & vecTime < 2*dblBinDur; % NOT COMPLETELY UNDERSTANDING WHAT THE POINT IS HERE...
indExcludeOff = vecTime > (1+ -2* dblBinDur) & vecTime < (1+ 2*dblBinDur);
indExclude = [find(indExcludeOn) find(indExcludeOff)];

%vecUnique = unique(vecSpikeCh);
%for intCh = 50:250%:length(vecUnique)
for intCh = 1:length(sAP.sCluster) % For each cluster:
    % vecSpikes = vecSpikeSecs(vecSpikeCh == vecUnique(intCh));
    vecSpikes = sAP.sCluster(intCh).SpikeTimes;
    dblZetaP = getZeta(vecSpikes,vecStimOnSecs(~vecLaserOn),0.9); % Compute zetatest for Stimuli w/o Opto
    if dblZetaP < 0.01 %&& sCluster(intCh).Violations1ms < 0.25 && abs(sCluster(intCh).NonStationarity) < 0.25 % ONLY plots figures for units that are VISUALLY RESPONSIVE (according to Zeta)
        figure; hold on;
        [vecMean,vecSEM,vecWindowBinCenters,~] = doPEP(vecSpikes,vecTime,vecStimOnSecs(~vecLaserOn),sOptions); %FROM JORRIT'S GENERALANALYSIS repo
        vecMean(indExclude) = NaN;
        plot(vecWindowBinCenters,vecMean,'k');
            % plot(vecWindowBinCenters,vecMean-vecSEM,'k--');
            % plot(vecWindowBinCenters,vecMean+vecSEM,'k--');
        [vecMean,vecSEM,vecWindowBinCenters,~] = doPEP(vecSpikes,vecTime,vecStimOnSecs(vecLaserOn),sOptions);
        vecMean(indExclude) = NaN;
        plot(vecWindowBinCenters,vecMean,'b');
            % plot(vecWindowBinCenters,vecMean-vecSEM,'b--');
            % plot(vecWindowBinCenters,vecMean+vecSEM,'b--');
        
        % title(['Channel: ' num2str(vecUnique(intCh))]);
        title(num2str(intCh));
        xline(0,'r--')
        ylabel('Rate (spks/s)')
        xlabel('Time relative to onset (s)')
        xlim([min(vecTime) max(vecTime)])
        legend('No Opto', 'Opto');
        fixfig;
        drawnow;
    else
        continue
    end
end

 % [~,vecUnique] = val2idx(vecSpikeCh);
 % for intCh = 100:200; %length(vecUnique)
 %     vecSpikes = vecSpikeSecs(vecSpikeCh == vecUnique(intCh));
 %     dblZetaP = zetatest(vecSpikes,vecStimOnSecs(~vecLaserOn),0.9);
 %     if dblZetaP < 0.05
 %         figure;subplot(2,1,1); hold on
 %        plotRaster(vecSpikes,vecStimOnSecs(~vecLaserOn),0.9,[],[0 0 0])
 % 
 % subplot(2,1,2);
 %         plotRaster(vecSpikes,vecStimOnSecs(vecLaserOn),0.9,[],[0 0 1])
 % 
 % 
 %         sgtitle(num2str(vecUnique(intCh)));
 % 
 % 
 %         drawnow;
 % 
 %        else
 %            continue
 %     end
 % end

 %% Export
