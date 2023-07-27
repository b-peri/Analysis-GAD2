%% OptoGratings Analysis

% To-DO:
%   [] Separate by Orientation? -> Not sure I'm seeing this one... >
%   Supplementary fig at most
%   [] Check SEM calculation...
%   [] Adapt script to extract LaserOn times from those rec that have it
%   but to skip for those who don't

%% Load in Data

% Load in Correct Data File
[FileNames, PathName] = uigetfile('*.mat', 'MultiSelect', 'on');

%% Set Up Struct

DataOut_OG.ClusterData = cell2table(cell(0,18), 'VariableNames', {'Subject', 'RecDate', 'ClusterN', 'zeta_p', 'Area', 'SpontRate', 'ER_OptoOff', ...
    'ER_OptoOn', 'SE_OptoOff', 'SE_OptoOn', 'PctChange', 'p_val', 'PSTHMean_Off', ...
    'PSTHSEM_Off', 'PSTHBinCenters_Off', 'PSTHMean_On', 'PSTHSEM_On', ...
    'PSTHBinCenters_On'});

DiffTimes_Onset = [];
DiffTimes_Offset = [];

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
intNumClu = length(sAP.sCluster); % Grabs number of clusters/putative cells
structEP = sAP.cellBlock{1,1};  % structEP contains data on a single recording block (incl. stim onset/offset times)
vecStimOnSecs = structEP.vecStimOnTime; % For OptoGratings, synchronized stimulus onset times
vecStimOffSecs = structEP.vecStimOffTime; % For Optogratings, synchronized stimulus OFFset times
vecLaserOn = structEP.vecOptoOn; % Logical array: Tells whether opto was on for each trial!
if isfield(structEP,"vecLaserOnTime")
    vecLaserOnTime = structEP.vecLaserOnTime; % Laser on times
    vecLaserOffTime = structEP.vecLaserOffTime; % Laser off times
    
    % Compute Differences & Append to Main Vec
    DiffTimes_Onset = [DiffTimes_Onset; (vecLaserOnTime - vecStimOnSecs(vecLaserOn))'];
    DiffTimes_Offset = [DiffTimes_Offset; (vecLaserOffTime - vecStimOffSecs(vecLaserOn))'];
end

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
Area = [];

%% Analysis

% vecROI = ["Superior colliculus zonal layer" "Superior colliculus superficial gray layer" "Superior colliculus optic layer"];

vecROI = ["Superior colliculus zonal layer" "Superior colliculus" + ...
    " superficial gray layer" "Superior colliculus optic layer" ...
    "Superior colliculus motor related intermediate gray layer"];

% -- Prep Analysis pt 1. --
% Spontaneous Rate = Mean rate in -1s - -100ms before stim onset
% Visually Evoked Rate = (Mean rate in 10-250ms after stim onset) - Spontaneous Rate
BinEdge = [-1 -0.1 0 0.2]; % Bin edges for Spike Counting w/ Histcounts
binDur = [BinEdge(2) - BinEdge(1), BinEdge(4) - BinEdge(3)]; % Bin Duration [Spontaneous, Evoked] -> Divisor when calculating Firing Rate later (Spikes/Period)
n_trials = numel(vecStimOnSecs); % Number of Trials

% -- Prep Analysis pt 2. --
dblBinDur = 5e-3; % Binsize of PSTH
vecTime = -0.2:dblBinDur:1.2; % PSTH X-Axis range/binsize
indExcludeOn = vecTime > -0.045 & vecTime < -0.015; % Exclude between -0.045 and -0.015 before stim onset (15ms before and after laser onset)
indExcludeOff = vecTime > 0.95 & vecTime < 0.98; % Ms before and after stim offset
indExclude = [find(indExcludeOn) find(indExcludeOff)];
sOptions.handleFig = -1;

% -- Start Loop --
for intCl = 1:length(sAP.sCluster) % For each cluster:
    vecSpikes = sAP.sCluster(intCl).SpikeTimes;
    dblZetaP = zetatest(vecSpikes,vecStimOnSecs(~vecLaserOn),0.9); % Compute zetatest for Stimuli w/o Opto -> Visually responsive neurons
    if dblZetaP < 0.01 && ismember(sAP.sCluster(intCl).Area, vecROI) && sAP.sCluster(intCl).Violations1ms < 0.25 %&& abs(sCluster(intCh).NonStationarity) < 0.25
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

        % T-test (Channel)
        [~,p_val_Cl] = ttest(EvokedRate_OptoOn_Cl, EvokedRate_OptoOff_Cl, 'Alpha', 0.01);
        
        % -- Analysis pt. 2: PSTH --
        % PSTH Laser Off
        [vecMean_Off,vecSEM_Off,vecWindowBinCenters_Off,~] = doPEP(vecSpikes,vecTime,vecStimOnSecs(~vecLaserOn),sOptions);
        %vecMean_Off(indExclude) = NaN;
        
        % PSTH Laser On 
        [vecMean_On,vecSEM_On,vecWindowBinCenters_On,~] = doPEP(vecSpikes,vecTime,vecStimOnSecs(vecLaserOn),sOptions);
        %vecMean_On(indExclude) = NaN;
        
        % -- Write to Vectors --
        ClusterN = [ClusterN; intCl];
        zeta_p = [zeta_p; dblZetaP];
        SpontRate = [SpontRate; SpontRate_Cl];
        ER_OptoOn = [ER_OptoOn; EvokedRate_OptoOn];
        ER_OptoOff = [ER_OptoOff; EvokedRate_OptoOff];
        SE_OptoOn = [SE_OptoOn; SEM_OptoOn];
        SE_OptoOff = [SE_OptoOff; SEM_OptoOff];
        PctChange = [PctChange; PctChange_Cl];
        p_val = [p_val; p_val_Cl];
        PSTHMean_On = [PSTHMean_On; vecMean_On]; % Note that all PSTH values are now baseline-subtracted!
        PSTHSEM_On = [PSTHSEM_On; vecSEM_On];
        PSTHBinCenters_On = [PSTHBinCenters_On; vecWindowBinCenters_On];
        PSTHMean_Off = [PSTHMean_Off; vecMean_Off];
        PSTHSEM_Off = [PSTHSEM_Off; vecSEM_Off];
        PSTHBinCenters_Off = [PSTHBinCenters_Off; vecWindowBinCenters_Off];
        Area = [Area; string(sAP.sCluster(intCl).Area)];
    else
        continue
    end
end

%% Write Output

n_clusters = numel(ClusterN);

% Write Table (Subject)
ClusterData_Rec = table(ClusterN, zeta_p, Area, SpontRate, ER_OptoOff, ...
    ER_OptoOn, SE_OptoOff, SE_OptoOn, PctChange, p_val, PSTHMean_Off, ...
    PSTHSEM_Off, PSTHBinCenters_Off, PSTHMean_On, PSTHSEM_On, ...
    PSTHBinCenters_On);

% Add to DataOut
SubjectN = table(repmat(string(sAP.sJson.subject), [n_clusters 1]), 'VariableNames', {'Subject'});
RecDate = table(repmat(string(sAP.sJson.date), [n_clusters 1]), 'VariableNames', {'RecDate'});
DataOut_OG.ClusterData = [DataOut_OG.ClusterData; [SubjectN RecDate ClusterData_Rec]];

end

 %% Overall Values

AllM_Overall = struct;

n_clusters_overall = height(DataOut_OG.ClusterData);
AllM_Overall.NMice = numel(unique(DataOut_OG.ClusterData.Subject));
AllM_Overall.NRecs = numel(unique([string(DataOut_OG.ClusterData.Subject) DataOut_OG.ClusterData.RecDate],"rows"));
AllM_Overall.NCells = n_clusters_overall; % Probably will need to tweak this!
AllM_Overall.NCells_Reduced = sum((DataOut_OG.ClusterData.p_val < 0.01) & DataOut_OG.ClusterData.PctChange < 0);
AllM_Overall.MeanPctReduction = mean(DataOut_OG.ClusterData.PctChange(DataOut_OG.ClusterData.PctChange < 0));
AllM_Overall.NCells_Increased = sum((DataOut_OG.ClusterData.p_val < 0.01) & DataOut_OG.ClusterData.PctChange > 0);
AllM_Overall.MeanPctIncrease = mean(DataOut_OG.ClusterData.PctChange((DataOut_OG.ClusterData.p_val<0.01) & (DataOut_OG.ClusterData.PctChange > 0)));
AllM_Overall.ER_OptoOn = mean(DataOut_OG.ClusterData.ER_OptoOn);
AllM_Overall.ER_OptoOff = mean(DataOut_OG.ClusterData.ER_OptoOff);
AllM_Overall.SE_OptoOn = std(DataOut_OG.ClusterData.ER_OptoOn, [], 1)/sqrt(n_clusters_overall);
AllM_Overall.SE_OptoOff = std(DataOut_OG.ClusterData.ER_OptoOff, [], 1)/sqrt(n_clusters_overall);
[~, AllM_Overall.p_val] = ttest(DataOut_OG.ClusterData.ER_OptoOn, DataOut_OG.ClusterData.ER_OptoOff, 'Alpha', 0.01);

% Filter PSTH Values on Min and Max Laser Times
PSTHBinCenters = PSTHBinCenters_Off(1,:);
exclidx = (PSTHBinCenters > (mean(DiffTimes_Onset) - 2*std(DiffTimes_Onset)) & PSTHBinCenters < (mean(DiffTimes_Onset) + 2*std(DiffTimes_Onset))) | ...
    (PSTHBinCenters > (1 + mean(DiffTimes_Offset) - 2*std(DiffTimes_Offset)) & PSTHBinCenters < (1 + mean(DiffTimes_Offset) + 2*std(DiffTimes_Offset)));
DataOut_OG.ClusterData.PSTHMean_Off(:,exclidx) = NaN;
DataOut_OG.ClusterData.PSTHMean_On(:,exclidx) = NaN;

% Z-Scored PSTH
PSTHoff_z = (DataOut_OG.ClusterData.PSTHMean_Off - mean(DataOut_OG.ClusterData.SpontRate))./std(DataOut_OG.ClusterData.SpontRate);
PSTHon_z = (DataOut_OG.ClusterData.PSTHMean_On - mean(DataOut_OG.ClusterData.SpontRate))./std(DataOut_OG.ClusterData.SpontRate);

AllM_Overall.PSTHMean_Off_z = mean(PSTHoff_z, 1);
AllM_Overall.PSTHSEM_Off_z = std(PSTHoff_z, 1)/sqrt(n_clusters_overall);
AllM_Overall.PSTHMean_On_z = mean(PSTHon_z, 1);
AllM_Overall.PSTHSEM_On_z = std(PSTHon_z, 1)/sqrt(n_clusters_overall);

% Normalized Overall PSTH Values
PSTHMean_Off_BL = DataOut_OG.ClusterData.PSTHMean_Off - DataOut_OG.ClusterData.SpontRate;
PSTHMean_On_BL = DataOut_OG.ClusterData.PSTHMean_On - DataOut_OG.ClusterData.SpontRate;
norm = max(PSTHMean_Off_BL,[],2);
AllM_Overall.PSTHMean_Off_norm = mean(PSTHMean_Off_BL./norm, 1);
AllM_Overall.PSTHSEM_Off_norm = std(PSTHMean_Off_BL./norm, 1)/sqrt(n_clusters_overall);
AllM_Overall.PSTHMean_On_norm = mean(PSTHMean_On_BL./norm, 1);
AllM_Overall.PSTHSEM_On_norm = std(PSTHMean_On_BL./norm, 1)/sqrt(n_clusters_overall);

AllM_Overall.PSTHBinSize = 5e-3; % Binsize of PSTH
AllM_Overall.PSTHtime = vecTime; % PSTH X-Axis range/binsize
AllM_Overall.PSTHBinCenters = DataOut_OG.ClusterData.PSTHBinCenters_On(1,:);

DataOut_OG.ClusterDataUp = DataOut_OG.ClusterData((DataOut_OG.ClusterData.p_val < 0.01) & (DataOut_OG.ClusterData.PctChange > 0),:);
DataOut_OG.ClusterDataDown = DataOut_OG.ClusterData((DataOut_OG.ClusterData.p_val < 0.01) & (DataOut_OG.ClusterData.PctChange < 0),:);
DataOut_OG.ClusterDataNonSig = DataOut_OG.ClusterData(~(DataOut_OG.ClusterData.p_val < 0.01),:);

DataOut_OG.Overall = AllM_Overall;
% DataOut_OG.LatencyOpto_On = DiffTimes_Onset;
% DataOut_OG.LatencyOpto_Off = DiffTimes_Offset;

%% Save Output

if all(startsWith(string(DataOut_OG.ClusterData.Subject),'7'))
    SaveFile = ['DataOut_OptoGratings_' datestr(datetime("today"),"dd-mm-yy") '_GAD2' '.mat'];
    DataOut_OG_GAD2 = DataOut_OG;
    save(SaveFile, 'DataOut_OG_GAD2');
elseif all(startsWith(string(DataOut_OG.ClusterData.Subject),'8'))
    SaveFile = ['DataOut_OptoGratings_' datestr(datetime("today"),"dd-mm-yy") '_Control' '.mat'];
    DataOut_OG_Cont = DataOut_OG;
    save(SaveFile, 'DataOut_OG_Cont');
end

clearvars -except DataOut_OT DataOut_OG_GAD2 DataOut_OG_Cont DataOut_OGt