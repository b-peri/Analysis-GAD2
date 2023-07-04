%% OptoGratings Analysis
% 
% Output should give:
% [] Per Cluster: Opto vs No-Opto
%       - Single Values (Spontaneous, Evoked Firing Rate (Opto & No-Opto)
%       - PSTH Values
% [] Per Mouse: NCells_SigReduced, AverageReduction_sig (In sSC), p_val 
% (for Opto vs No-Opto), vecMean_mean, vecSEM_mean (HOW? -> Currently SEM over means), 
% vecBinCenters, n_trials, n_clusters

% OUTPUT -> Currently, Output is a struct w/ output data for a single 
% mouse. Ultimately want this to be a cell array/table where each entry 
% corresponds to an individual subject/mouse
% 
% For each SubjectN = Cell? Table? Struct?;
%   [] Struct 1: Single-Unit (Per Cluster)
%       [] Find a way to pre-allocate vector size
%       [] ULTIMATELY alignment data -> Need to know whether clusters are
%       located in sSC
%   [X] Struct 2: Overall
  
% [] LATER -> Struct 3 (Overall over ALL mice included in analysis)
%       Should report: 1. Number of mice 2. Total number of
%       neurons/clusters included
%       [] Should this be mean of means per mouse or per cluster?


%% Load in Data

% Load in Correct Data File
% DataDirectory = dir(fullfile("C:\Software and Code\Analysis-GAD2\example_data_scripts\*.mat"));
% FileNames = {DataDirectory.name};
% sAP_Files = FileNames(endsWith(FileNames, '_Synthesis.mat'));
% load(fullfile(DataDirectory(1).folder,sAP_Files{1})); % [] This currently just grabs the first folder in DataFile. May want to change this to be more specific later!

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

% Info to Have Per Mouse: 

DataOut.Subject = sAP.sJson.subject;

% Table Headers: ClusterN, zeta, SpontRate, EvokedRate_OptoOn,
% EvokedRate_OptoOff, SEM_On, SEM_Off
ClusterN = [];
zeta_p = [];
SpontRate = [];
ER_OptoOn = [];
ER_OptoOff = [];
SE_OptoOn = [];
SE_OptoOff = [];
PctChange = [];
p_val = [];
PSTHMean_On = [];
PSTHSEM_On = [];
PSTHBinCenters_On = [];
PSTHMean_Off = [];
PSTHSEM_Off = [];
PSTHBinCenters_Off = [];

%% Analysis

% -- Prep Analysis pt 1. --
% Spontaneous Rate = Mean rate in -1s - -100ms before stim onset
% Visually Evoked Rate = Mean rate in 10-250ms after stim onset - Spontaneous Rate
BinEdge = [-1 -0.1 0.01 0.250]; % Bin edges for Spike Counting w/ Histcounts
binDur = [BinEdge(2) - BinEdge(1), BinEdge(4) - BinEdge(3)]; % Bin Duration [Spontaneous, Evoked] -> Divisor when calculating Firing Rate later (Spikes/Period)
ExclWindowSRate = [vecStimOnSecs'-1 vecStimOnSecs'-0.1]; % Excl. 100ms before stimulus onset -> Remove artefact
n_trials = numel(vecStimOnSecs); % Number of Trials

% -- Prep Analysis pt 2. --
dblBinDur = 5e-3; % Binsize of PSTH
vecTime = -0.2:dblBinDur:1.2; % PSTH X-Axis range/binsize
indExcludeOn = vecTime > -2* dblBinDur & vecTime < 2*dblBinDur; % First millisecond before and after stim onset
indExcludeOff = vecTime > (1+ -2* dblBinDur) & vecTime < (1+ 2*dblBinDur); % Ms before and after stim offset
indExclude = [find(indExcludeOn) find(indExcludeOff)];
sOptions = -1;

% -- Start Loop --
for intCh = 1:length(sAP.sCluster) % For each cluster:
    vecSpikes = sAP.sCluster(intCh).SpikeTimes;
    dblZetaP = zetatest(vecSpikes,vecStimOnSecs(~vecLaserOn),0.9); % Compute zetatest for Stimuli w/o Opto -> Visually responsive neurons
    if dblZetaP < 0.01 %&& sCluster(intCh).Violations1ms < 0.25 && abs(sCluster(intCh).NonStationarity) < 0.25 % ONLY plots figures for units that are VISUALLY RESPONSIVE (according to Zeta)
        % -- Analysis pt 1.  Opto vs No-Opto --
        sCounts_Opto = zeros(numel(vecStimOnSecs),2);
        for intTrial=1:n_trials
	        vecTheseEdges = BinEdge + vecStimOnSecs(intTrial); % Add stim onset time for this trial to rel. bin edges to get absolute bin edges
	        [vecCounts,edges] = histcounts(vecSpikes,vecTheseEdges);
	        sCounts_Opto(intTrial,1) = vecCounts(1); % Counts Spontaneous Rate
            sCounts_Opto(intTrial,2) = vecCounts(3); % Count for Visual Response
        end
        
        sRate_Opto = [sCounts_Opto(:,1)/binDur(1) sCounts_Opto(:,2)/binDur(2)];
        sRate_OptoOn = sRate_Opto(vecLaserOn, 2); % Rates during stim for Opto trials
        sRate_OptoOff = sRate_Opto(~vecLaserOn, 2); % Rates during stim non-Opto trials
        
        % Spontaneous Rate
        SpontRate_Cl = mean(sCounts_Opto(:,1)); 
    
        % Overall Visually Evoked FRs (Channel)
        EvokedRate_OptoOn = mean(sRate_OptoOn) - SpontRate_Cl;
        EvokedRate_OptoOff = mean(sRate_OptoOff) - SpontRate_Cl;

        % Standard Error
        SEM_OptoOn = std(sRate_OptoOn, [], 1)/sqrt(n_trials);
        SEM_OptoOff = std(sRate_OptoOff, [], 1)/sqrt(n_trials);

        % Magnitude of Reduction per Channel
        PctChange_Cl = (EvokedRate_OptoOn - EvokedRate_OptoOff)/EvokedRate_OptoOff;

        % EvokedRate per Trial
        EvokedRate_OptoOn_Cl = sRate_OptoOn - SpontRate_Cl;
        EvokedRate_OptoOff_Cl = sRate_OptoOff - SpontRate_Cl;

        % Paired-Sample T-test (Channel) -> THIS IS DEFINITELY WEIRD ?
        [~,p_val_Cl] = ttest(EvokedRate_OptoOn_Cl, EvokedRate_OptoOff_Cl, 'Alpha', 0.01);
        
        % -- Analysis pt. 2: PSTH --
        % PSTH Laser Off
        [vecMean_Off,vecSEM_Off,vecWindowBinCenters_Off,~] = doPEP(vecSpikes,vecTime,vecStimOnSecs(~vecLaserOn),sOptions);
        vecMean_Off(indExclude) = NaN;
        
        % PSTH Laser On 
        [vecMean_On,vecSEM_On,vecWindowBinCenters_On,~] = doPEP(vecSpikes,vecTime,vecStimOnSecs(vecLaserOn),sOptions);
        vecMean_On(indExclude) = NaN;
        
        % -- Write to Vectors --
        ClusterN = [ClusterN; intCh];
        zeta_p = [zeta_p; dblZetaP];
        SpontRate = [SpontRate; SpontRate_Cl];
        ER_OptoOn = [ER_OptoOn; EvokedRate_OptoOn];
        ER_OptoOff = [ER_OptoOff; EvokedRate_OptoOff];
        SE_OptoOn = [SE_OptoOn; SEM_OptoOn];
        SE_OptoOff = [SE_OptoOff; SEM_OptoOff];
        PctChange = [PctChange; PctChange_Cl];
        p_val = [p_val; p_val_Cl];
        PSTHMean_On = [PSTHMean_On; vecMean_On];
        PSTHSEM_On = [PSTHSEM_On; vecSEM_On];
        PSTHBinCenters_On = [PSTHBinCenters_On; vecWindowBinCenters_On];
        PSTHMean_Off = [PSTHMean_Off; vecMean_Off];
        PSTHSEM_Off = [PSTHSEM_Off; vecSEM_Off];
        PSTHBinCenters_Off = [PSTHBinCenters_Off; vecWindowBinCenters_Off];
    else
        continue
    end
end        

%% Analysis pt 3: Overall Values

n_clusters = numel(ClusterN);

Overall.NCells_Reduced = sum((p_val < 0.01) & PctChange < 0);
Overall.MeanPctReduction = mean(PctChange(PctChange < 0));
Overall.NCells_Increased = sum((p_val < 0.01) & PctChange > 0);
Overall.MeanPctIncrease = mean(PctChange((p_val<0.01) & (PctChange > 0)));
Overall.ER_OptoOn = mean(ER_OptoOn);
Overall.ER_OptoOff = mean(ER_OptoOff);
Overall.SE_OptoOn = std(ER_OptoOn, [], 1)/sqrt(n_clusters);
Overall.SE_OptoOff = std(ER_OptoOff, [], 1)/sqrt(n_clusters);
[~, Overall.p_val] = ttest(ER_OptoOn, ER_OptoOff, 'Alpha', 0.01); % Inclusion criteria? Only negatively modulated cells?

norm = max(PSTHMean_Off,[],2);
Overall.PSTHMean_Off_norm = mean(PSTHMean_Off./norm, 1);
Overall.PSTHSEM_Off_norm = std(PSTHMean_Off./norm, 1)/sqrt(n_clusters);
Overall.PSTHMean_On_norm = mean(PSTHMean_On./norm, 1);
Overall.PSTHSEM_On_norm = std(PSTHMean_On./norm, 1)/sqrt(n_clusters);
Overall.PSTHBinSize = 5e-3; % Binsize of PSTH
Overall.PSTHtime = vecTime; % PSTH X-Axis range/binsize
Overall.PSTHBinCenters = PSTHBinCenters_On(1,:);

Overall.n_trials = n_trials;
Overall.n_clusters = n_clusters;


%% Write Output

% Write Table
DataOut.ClusterData = table(ClusterN, zeta_p, SpontRate, ER_OptoOff, ...
    ER_OptoOn, SE_OptoOff, SE_OptoOn, PctChange, p_val, PSTHMean_Off, ...
    PSTHSEM_Off, PSTHBinCenters_Off, PSTHMean_On, PSTHSEM_On, ...
    PSTHBinCenters_On);

% Overall Subject Data
DataOut.OverallData = Overall;

%% Test Plots

DataOut.ClusterSig = DataOut.ClusterData(p_val < 0.01,:);

PlotClusters = 0;

if PlotClusters
    for intCh = DataOut.ClusterSig.ClusterN'
        row = DataOut.ClusterSig(DataOut.ClusterSig.ClusterN == intCh, :);
        vecSpikes = sAP.sCluster(intCh).SpikeTimes;
        figure; hold on;
        plot(row.PSTHBinCenters_Off,row.PSTHMean_Off,'k');
            % plot(vecWindowBinCenters,vecMean-vecSEM,'k--');
            % plot(vecWindowBinCenters,vecMean+vecSEM,'k--');
        plot(row.PSTHBinCenters_On,row.PSTHMean_On,'b');
            % plot(vecWindowBinCenters,vecMean-vecSEM,'b--');
            % plot(vecWindowBinCenters,vecMean+vecSEM,'b--');
        
        % title(['Channel: ' num2str(vecUnique(intCh))]);
        title(num2str(intCh));
        xline(0,'r--')
        ylabel('Spiking Rate (spks/s)')
        xlabel('Time from stimulus onset (s)')
        xlim([min(vecTime) max(vecTime)])
        legend('No Opto', 'Opto');
        fixfig;
        drawnow;
        hold off;
    end
end

% Barplots
figure; hold on;
bar_pl = bar([Overall.ER_OptoOff Overall.ER_OptoOn]);
dot_pl = plot([1 2],[ER_OptoOff ER_OptoOn],'-','Color',[0, 0, 0, 0.3]);
scatter([1 2],[ER_OptoOff ER_OptoOn],'MarkerEdgeColor', [0 0 0], 'MarkerEdgeAlpha',0.5);
bar_err = errorbar([1,2],[Overall.ER_OptoOff Overall.ER_OptoOn], [Overall.SE_OptoOff Overall.SE_OptoOn]);
bar_err.Color = [0.8 0 0];
bar_err.LineStyle = 'none';
bar_err.LineWidth = 1.5;
ylabel('Evoked Rate (spks/s)');
xticks([1 2]);
xticklabels({"Laser OFF", "Laser ON"});
hold off

% Mean Response
figure; hold on;
plot(DataOut.OverallData.PSTHBinCenters, DataOut.OverallData.PSTHMean_Off_norm,'k');
plot(DataOut.OverallData.PSTHBinCenters, DataOut.OverallData.PSTHMean_On_norm,'b');
xline(0,'r--');
ylabel('Spiking Rate (% of Peak)');
xlabel('Time from stimulus onset (s)');
xlim([min(vecTime) max(vecTime)]);
text(0.94, 0.6, sprintf('%g units', n_clusters), 'FontSize',15);
legend('No Opto', 'Opto');
fixfig;
drawnow;
hold off;

 %% Export
