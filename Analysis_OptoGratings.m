%% OptoGratings Analysis

% Output should give:
% [X] Per Cluster: Opto vs No-Opto
%       - Single Values (Spontaneous, Evoked Firing Rate (Opto & No-Opto)
%       - PSTH Values
% [X] Per Mouse: NCells_Reduced, AverageReduction_sig (In sSC), NCells_Increased, AverageIncrease, p_val 
% (for Opto vs No-Opto), vecMean_mean, vecSEM_mean (HOW? -> Currently SEM over means), 
% vecBinCenters, n_trials, n_clusters
% 
% For each SubjectN = Cell? Table? Struct?;
%   [X] Struct 1: Single-Unit (Per Cluster)
%       [] ULTIMATELY alignment data -> Need to know whether clusters are
%       located in sSC
%   [X] Struct 2: Overall
  
% [X] Struct 3 (Overall over ALL mice included in analysis)


%% Load in Data

% Load in Correct Data File
[FileNames, PathName] = uigetfile('*.mat', 'MultiSelect', 'on');
% DataDirectory = dir(fullfile("C:\Software and Code\Analysis-GAD2\example_data_scripts\*.mat"));
% FileNames = {DataDirectory.name};
% sAP_Files = FileNames(endsWith(FileNames, '_Synthesis.mat'));
% load(fullfile(DataDirectory(1).folder,sAP_Files{1})); % [] This currently just grabs the first folder in DataFile. May want to change this to be more specific later!

%% Set Up Struct

DataOut.AllMice.ClusterData = cell2table(cell(0,17), 'VariableNames', {'Subject', 'RecDate', 'ClusterN', 'zeta_p', 'SpontRate', 'ER_OptoOff', ...
    'ER_OptoOn', 'SE_OptoOff', 'SE_OptoOn', 'PctChange', 'p_val', 'PSTHMean_Off', ...
    'PSTHSEM_Off', 'PSTHBinCenters_Off', 'PSTHMean_On', 'PSTHSEM_On', ...
    'PSTHBinCenters_On'});

%% Start Loop

if isa(FileNames, 'cell')
    NumFiles = numel(FileNames);
else
    NumFiles = 1;
end

for idx = 1:NumFiles % For each recording
if isa(FileNames, 'cell')
    load(fullfile(PathName, FileNames{idx}));
else
    load(fullfile(PathName, FileNames));
end

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
RecData.Subject = sAP.sJson.subject;
RecData.SubjectType = sAP.sJson.subjecttype;

% Table Headers: ClusterN, zeta_p, SpontRate, EvokedRate_OptoOn,
% EvokedRate_OptoOff, SEM_On, SEM_Off, PctChange, p_val, PSTH stuff ...
ClusterN = [];
zeta_p = [];
SpontRate = [];
ER_OptoOn = [];
ER_OptoOff = [];
SE_OptoOn = [];
SE_OptoOff = [];
PctChange = [];
p_val = [];
PSTHMean_Off = [];
PSTHSEM_Off = [];
PSTHBinCenters_Off = [];
PSTHMean_On = [];
PSTHSEM_On = [];
PSTHBinCenters_On = [];

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
        % SpontRate_Cl = mean(sCounts_Opto(:,1));
        SpontRate_Cl = mean(sRate_Opto(:,1));
    
        % Overall Visually Evoked FRs (Channel)
        EvokedRate_OptoOn = mean(sRate_OptoOn) - SpontRate_Cl;
        EvokedRate_OptoOff = mean(sRate_OptoOff) - SpontRate_Cl;

        % Standard Error
        SEM_OptoOn = std(sRate_OptoOn, [], 1)/sqrt(n_trials);
        SEM_OptoOff = std(sRate_OptoOff, [], 1)/sqrt(n_trials);

        % Magnitude of Reduction per Channel
        PctChange_Cl = (EvokedRate_OptoOn - EvokedRate_OptoOff)/abs(EvokedRate_OptoOff);

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
        PSTHMean_On = [PSTHMean_On; vecMean_On - SpontRate_Cl]; % Note that all PSTH values are now baseline-subtracted!
        PSTHSEM_On = [PSTHSEM_On; vecSEM_On];
        PSTHBinCenters_On = [PSTHBinCenters_On; vecWindowBinCenters_On];
        PSTHMean_Off = [PSTHMean_Off; vecMean_Off - SpontRate_Cl];
        PSTHSEM_Off = [PSTHSEM_Off; vecSEM_Off];
        PSTHBinCenters_Off = [PSTHBinCenters_Off; vecWindowBinCenters_Off];
    else
        continue
    end
end

%% Analysis pt 3: Overall Values per Recording

n_clusters = numel(ClusterN);

RecOverall.NCells_Reduced = sum((p_val < 0.01) & PctChange < 0);
RecOverall.MeanPctReduction = mean(PctChange(PctChange < 0)); % 
RecOverall.NCells_Increased = sum((p_val < 0.01) & PctChange > 0);
RecOverall.MeanPctIncrease = mean(PctChange((p_val<0.01) & (PctChange > 0)));
RecOverall.ER_OptoOn = mean(ER_OptoOn);
RecOverall.ER_OptoOff = mean(ER_OptoOff);
RecOverall.SE_OptoOn = std(ER_OptoOn, [], 1)/sqrt(n_clusters);
RecOverall.SE_OptoOff = std(ER_OptoOff, [], 1)/sqrt(n_clusters);
[~, RecOverall.p_val] = ttest(ER_OptoOn, ER_OptoOff, 'Alpha', 0.01);

% Non-Normalized PSTH Values
norm = max(PSTHMean_Off,[],2);
RecOverall.PSTHMean_Off = mean(PSTHMean_Off./norm, 1);
RecOverall.PSTHSEM_Off = std(PSTHMean_Off./norm, 1)/sqrt(n_clusters);
RecOverall.PSTHMean_On = mean(PSTHMean_On./norm, 1);
RecOverall.PSTHSEM_On = std(PSTHMean_On./norm, 1)/sqrt(n_clusters);

%Normalized PSTH Values
norm = max(PSTHMean_Off,[],2);
RecOverall.PSTHMean_Off_norm = mean(PSTHMean_Off./norm, 1);
RecOverall.PSTHSEM_Off_norm = std(PSTHMean_Off./norm, 1)/sqrt(n_clusters);
RecOverall.PSTHMean_On_norm = mean(PSTHMean_On./norm, 1);
RecOverall.PSTHSEM_On_norm = std(PSTHMean_On./norm, 1)/sqrt(n_clusters);
RecOverall.PSTHBinSize = 5e-3; % Binsize of PSTH
RecOverall.PSTHtime = vecTime; % PSTH X-Axis range/binsize
RecOverall.PSTHBinCenters = PSTHBinCenters_On(1,:);

RecOverall.n_trials = n_trials;
RecOverall.n_clusters = n_clusters;


%% Write Output

% Write Table (Subject)
RecData.ClusterData = table(ClusterN, zeta_p, SpontRate, ER_OptoOff, ...
    ER_OptoOn, SE_OptoOff, SE_OptoOn, PctChange, p_val, PSTHMean_Off, ...
    PSTHSEM_Off, PSTHBinCenters_Off, PSTHMean_On, PSTHSEM_On, ...
    PSTHBinCenters_On);

% Overall Subject Data
RecData.OverallData = RecOverall;

% Add to DataOut
SubjectN = table(repmat(RecData.Subject, [n_clusters 1]), 'VariableNames', {'Subject'});
RecDate = cell2table(repmat({sAP.sJson.date}', [n_clusters 1]), 'VariableNames', {'RecDate'});
% ADD CHANNEL DEPTH!
DataOut.AllMice.ClusterData = [DataOut.AllMice.ClusterData; [SubjectN RecDate RecData.ClusterData]];

RecordingName = [replace(sAP.sJson.experiment(1:end-6),'-','_')];
DataOut.(RecordingName) = RecData;

end

 %% Overall Values

AllM_Overall = struct;

n_clusters_overall = height(DataOut.AllMice.ClusterData);
AllM_Overall.NMice = numel(unique(DataOut.AllMice.ClusterData(:,1)));
AllM_Overall.NRecs = numel(fieldnames(DataOut)) - 1;
AllM_Overall.NCells = n_clusters_overall; % Probably will need to tweak this!
AllM_Overall.NCells_Reduced = sum((DataOut.AllMice.ClusterData.p_val < 0.01) & DataOut.AllMice.ClusterData.PctChange < 0);
AllM_Overall.MeanPctReduction = mean(DataOut.AllMice.ClusterData.PctChange(DataOut.AllMice.ClusterData.PctChange < 0));
AllM_Overall.NCells_Increased = sum((DataOut.AllMice.ClusterData.p_val < 0.01) & DataOut.AllMice.ClusterData.PctChange > 0);
AllM_Overall.MeanPctIncrease = mean(DataOut.AllMice.ClusterData.PctChange((DataOut.AllMice.ClusterData.p_val<0.01) & (DataOut.AllMice.ClusterData.PctChange > 0)));
AllM_Overall.ER_OptoOn = mean(DataOut.AllMice.ClusterData.ER_OptoOn);
AllM_Overall.ER_OptoOff = mean(DataOut.AllMice.ClusterData.ER_OptoOff);
AllM_Overall.SE_OptoOn = std(DataOut.AllMice.ClusterData.ER_OptoOn, [], 1)/sqrt(n_clusters_overall);
AllM_Overall.SE_OptoOff = std(DataOut.AllMice.ClusterData.ER_OptoOff, [], 1)/sqrt(n_clusters_overall);
[~, AllM_Overall.p_val] = ttest(DataOut.AllMice.ClusterData.ER_OptoOn, DataOut.AllMice.ClusterData.ER_OptoOff, 'Alpha', 0.01);

% Non-Normalized 
AllM_Overall.PSTHMean_Off = mean(DataOut.AllMice.ClusterData.PSTHMean_Off, 1);
AllM_Overall.PSTHSEM_Off = std(DataOut.AllMice.ClusterData.PSTHMean_Off, 1)/sqrt(n_clusters_overall);
AllM_Overall.PSTHMean_On = mean(DataOut.AllMice.ClusterData.PSTHMean_On, 1);
AllM_Overall.PSTHSEM_On = std(DataOut.AllMice.ClusterData.PSTHMean_On, 1)/sqrt(n_clusters_overall);

% Normalized Overall PSTH Values
norm = max(DataOut.AllMice.ClusterData.PSTHMean_Off,[],2);
AllM_Overall.PSTHMean_Off_norm = mean(DataOut.AllMice.ClusterData.PSTHMean_Off./norm, 1);
AllM_Overall.PSTHSEM_Off_norm = std(DataOut.AllMice.ClusterData.PSTHMean_Off./norm, 1)/sqrt(n_clusters_overall);
AllM_Overall.PSTHMean_On_norm = mean(DataOut.AllMice.ClusterData.PSTHMean_On./norm, 1);
AllM_Overall.PSTHSEM_On_norm = std(DataOut.AllMice.ClusterData.PSTHMean_On./norm, 1)/sqrt(n_clusters_overall);

AllM_Overall.PSTHBinSize = 5e-3; % Binsize of PSTH
AllM_Overall.PSTHtime = vecTime; % PSTH X-Axis range/binsize
AllM_Overall.PSTHBinCenters = DataOut.AllMice.ClusterData.PSTHBinCenters_On(1,:);

DataOut.AllMice.ClusterDataUp = DataOut.AllMice.ClusterData((DataOut.AllMice.ClusterData.p_val < 0.01) & (DataOut.AllMice.ClusterData.PctChange > 0),:);
DataOut.AllMice.ClusterDataDown = DataOut.AllMice.ClusterData((DataOut.AllMice.ClusterData.p_val < 0.01) & (DataOut.AllMice.ClusterData.PctChange < 0),:);
DataOut.AllMice.ClusterDataNonSig = DataOut.AllMice.ClusterData(~(DataOut.AllMice.ClusterData.p_val < 0.01),:);

DataOut.AllMice.Overall = AllM_Overall;

%% Save Output
if any(startsWith(string(fieldnames(DataOut)),'Rec7'))
    SaveFile = ['DataOut_OptoGratings_' datestr(datetime("today"),"dd-mm-yy") '_GAD2' '.mat'];
elseif any(startsWith(string(fieldnames(DataOut)),'Rec8'))
    SaveFile = ['DataOut_OptoGratings_' datestr(datetime("today"),"dd-mm-yy") '_Control' '.mat'];
end
save(SaveFile, 'DataOut');