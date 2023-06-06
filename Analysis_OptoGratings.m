%% OptoGratings Analysis
% 
% Robin's example script makes PLOTS, but what we want is an OUTPUT FILE
% that should give:
% [] Per Cluster: Optograting PSTH
% [] Per Mouse: Some metric of OVERALL reduction in sSC specifically (in firing rate?)
% [] OPTIONAL: Conversion to CSV

% OUTPUT -> Currently, Output is a struct w/ output data for a single 
% mouse. Ultimately want this to be a cell array/table where each entry 
% corresponds to an individual subject/mouse
% 
% For each MouseN = struct;
%   [] Struct 1: Multi-Unit (Per Cluster, per Mouse)
%       [X] Zeta-Test value per cluster
%       [] Spontaneous Rate, Evoked Rate (Opto), Evoked Rate (No-Opto)
%       [] Repeated measures T-test value PER CLUSTER
%           [] Also want to know HOW MANY clusters are significantly
%           _reduced_ by Opto
%           [] Ask Robin about what test to use and what p-threshold should
%           be... Do we want a directional test (and then test separately
%           for reduction and increase?)
%       [] ULTIMATELY (Low-Priority), would maybe also be nice to have
%       clusters sorted on anatomical location, and have some value for
%       where they are exactly (from alignment)
%   [] ?? -> Struct 2: Single-Unit (Per Cluster, per Mouse)
%   [] Struct 2/3: Overall
%       [] T-test value for MEAN REDUCTION over ALL (sSC)
%       channels
%   [] Struct 3 (Overall over ALL mice included in analysis)
%       Should report: 1. Number of mice 2. Total number of neurons
%       included
%       [] Should this be mean of means per mouse or per cluster?
% [] Figure out where/whether it says somewhere if clusters are MU or SU
% [] Separate Figs for SINGLE-UNIT and MULTI-UNIT data?


%% 

% Load in Correct Data File
DataFiles = dir(fullfile("C:\Software and Code\Analysis-GAD2\example_data_scripts\*.mat"));
FileNames = {DataFiles.name};
sAP_Files = FileNames(endsWith(FileNames, '_Synthesis.mat'));
load(fullfile(DataFiles(1).folder,sAP_Files{1})); % [] This currently just grabs the first folder in DataFile. May want to change this to be more specific later!

%%
% NOTE: SCRIPT IS CURRENTLY WRITTEN TO HANDLE ONE MOUSE AT A TIME, CHANGE
% THIS TO LOOP ONCE BASIC ANALYSIS HAS BEEN IMPLEMENTED

% Get stimulus info
sAP = sSynthData;
intNumClu = length(sAP.sCluster); % Grabs number of clusters/putative cells
structEP = sAP.cellStim{1,1}.structEP;  % structEP contains data on a single recording block (incl. stim onset/offset times)
% -> This index is hardcoded. Not a huge deal if OptoGratings is always the
% first block in the recording, but not very pretty/generalizable
vecStimOnSecs = structEP.vecStimOnTime; % For OptoGratings, synchronized stimulus onset times
vecStimOffSecs = structEP.vecStimOffTime; % For Optogratings, synchronized stimulus OFFset times
vecLaserOn = structEP.vecOptoOn; % Logical array: Tells whether opto was on for each trial!

%% Prepare Output Table

Subject = sAP.sJson.subject
% Output = table();
DataOut = struct;
DataOut.Cluster = table;
% Output.Mouse = struct;

%% Analysis 1: Opto vs No-Opto

% Spontaneous Rate = Mean rate in -1s - -100ms before stim onset
% Visually Evoked Rate = Mean rate in 0-1s after stim onset (entire stim
% presentation) - Spontaneous Rate

sOptions.handleFig = -1; % COME BACK TO THIS
vecTime = [-1 -0.1 0 1] % Bin edges for Spike Counting w/ Histcounts
binDur = [vecTime(2) - vecTime(1), vecTime(4) - vecTime(3)] % Binsize [Spontaneous, Evoked] 
ExclWindowSRate = [vecStimOnSecs'-1 vecStimOnSecs'-0.1] % Spontaneous Rate

for intCh = 1:length(sAP.sCluster) % For each cluster:
    vecSpikes = sAP.sCluster(intCh).SpikeTimes;
    dblZetaP = zetatest(vecSpikes,vecStimOnSecs(~vecLaserOn),0.9); % Compute zetatest for Stimuli w/o Opto -> Visually responsive neurons
    if dblZetaP < 0.01 %&& sCluster(intCh).Violations1ms < 0.25 && abs(sCluster(intCh).NonStationarity) < 0.25 % ONLY plots figures for units that are VISUALLY RESPONSIVE (according to Zeta)
        % Opto Inhibition
        sCounts_Opto = zeros(numel(vecStimOnSecs),2);
        for intTrial=1:numel(vecStimOnSecs)
	        vecTheseEdges = vecTime + vecStimOnSecs(intTrial); % Add stim onset time for this trial to rel. bin edges to get absolute bin edges
	        [vecCounts,edges] = histcounts(vecSpikes,vecTheseEdges);
	        sCounts_Opto(intTrial,1) = vecCounts(1); % Counts Spontaneous Rate
            sCounts_Opto(intTrial,2) = vecCounts(3); % Count for Visual Response
        end
        
        sRate_Opto = [sCounts_Opto(:,1)/binDur(1) sCounts_Opto(:,2)/binDur(2)];
        sRate_OptoOn = sRate_Opto(vecLaserOn, 2); % Evoked rates for Opto trials
        sRate_OptoOff = sRate_Opto(~vecLaserOn, 2); % Evoked rates non-Opto trials
        
        SpontRate = mean(sCounts_Opto(:,1)); 
        EvokedRate_OptoOn = mean(sRate_OptoOn) - SpontRate;
        EvokedRate_OptoOff = mean(sRate_OptoOff) - SpontRate;

        SEM_OptoOn = std(sRate_OptoOn, [], 1)
        SEM_OptoOff = std(sRate_OptoOff, [], 1)

        % Paired-Sample T-test
        

        % Write to Table
        % ...{'Zeta'} = 
    else
        continue
    end
end        
   
%% Plots

% dblBinDur = 5e-3; % Binsize of PSTH
% vecTime = -0.2:dblBinDur:1.2; % PSTH X-Axis range/binsize
% indExcludeOn = vecTime > -2* dblBinDur & vecTime < 2*dblBinDur; % First millisecond before and after stim onset
% indExcludeOff = vecTime > (1+ -2* dblBinDur) & vecTime < (1+ 2*dblBinDur); % Ms before and after stim offset
% indExclude = [find(indExcludeOn) find(indExcludeOff)];

% This was originally in cluster loop above
        % % PSTH Laser Off
        % [vecMean_Off,vecSEM_Off,vecWindowBinCenters_Off,~] = doPEP(vecSpikes,vecTime,vecStimOnSecs(~vecLaserOn),sOptions); %FROM JORRIT'S GENERALANALYSIS repo
        % vecMean_Off(indExclude) = NaN;
        % % PSTH Laser On
        % [vecMean_On,vecSEM_On,vecWindowBinCenters_On,~] = doPEP(vecSpikes,vecTime,vecStimOnSecs(vecLaserOn),sOptions);
        % vecMean_On(indExclude) = NaN;        
% plot(vecWindowBinCenters,vecMean,'k');
        %     % plot(vecWindowBinCenters,vecMean-vecSEM,'k--');
        %     % plot(vecWindowBinCenters,vecMean+vecSEM,'k--');
        
        % plot(vecWindowBinCenters,vecMean,'b');
        %     % plot(vecWindowBinCenters,vecMean-vecSEM,'b--');
        %     % plot(vecWindowBinCenters,vecMean+vecSEM,'b--');
        % 
        % % title(['Channel: ' num2str(vecUnique(intCh))]);
        % title(num2str(intCh));
        % xline(0,'r--')
        % ylabel('Rate (spks/s)')
        % xlabel('Time relative to onset (s)')
        % xlim([min(vecTime) max(vecTime)])
        % legend('No Opto', 'Opto');
        % fixfig;
        % drawnow;

 %% Export
