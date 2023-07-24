% Analysis Optotagging GAD2

% Note for resp dynamics: We want to see 1. Optotagged GAD2 units, 2.
% Inhibited units, 3. (Putative) disinhibited units! (and maybe 4.
% Other/visually responsive units?)

% To-Do

%%

%Prep Table
DataOut_OT.AllMice.ClusterData = cell2table(cell(0,13), 'VariableNames', ...
    {'Subject', 'RecDate', 'ClusterN', 'CellClass', 'Area', 'zeta_p_20ms', ...
    'Peak_Lat', 'IFR_Rate', 'IFR_Time', 'SpontRate', 'SpontRate_STD', ...
    'PSTHMean_20ms', 'PSTHSEM_20ms', 'PSTHBinCenters'});

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
%% Grab Stimulus Data

intNumClu = length(sAP.sCluster);
structEP = sAP.cellBlock{1,3};
vecLaserOnSecs = structEP.vecLaserOnTime;
pulseDurTrial = structEP.PulseDurTrial;
% pulseDurations = structEP.vecPulseDuration;
% vecStimOffSecs = structEP.vecLaserOffTime;

%% Prep Output Table

% Initial Data RecData
RecData.SubjectType = sAP.sJson.subjecttype;

% MouseN, RecN, ClusterN, zeta_p, peak latency, Instaneous FR, mean PSTH values (probably
% w/ 2 bin sizes?)
MouseN = sAP.sJson.subject;
RecN = sAP.sJson.date;
ClusterN = [];
CellClass = []
Area = [];
zeta_p_20ms = [];
Peak_Lat = [];
IFR_Rate = [];
IFR_Time = [];
SpontRate = [];
SpontRate_STD = [];

PSTHMean_20ms = [];
PSTHSEM_20ms = [];
PSTHBinCenters = [];

%% Compute Zeta, Get Latencies, and Get Average Responses

% vecROI = ["Superior colliculus zonal layer" "Superior colliculus" + ...
%     " superficial gray layer" "Superior colliculus optic layer"];

vecROI = ["Superior colliculus zonal layer" "Superior colliculus" + ...
    " superficial gray layer" "Superior colliculus optic layer" ...
    "Superior colliculus motor related intermediate gray layer" ...
    "Superior colliculus motor related intermediate white layer"];

% --- Select Pulse Duration ---
% pulseDur = 0.02;
% n_trials_20ms = sum(pulseDurTrial == 0.02);
% vecLaserOnSecs_dur = vecLaserOnSecs(pulseDurTrial == pulseDurTrial);

% --- Prep Overall Modulation (Baseline vs. Evoked) ---
% Spontaneous Rate = Mean rate in -1s - -100ms before stim onset
% Visually Evoked Rate = Mean rate in 10-250ms after stim onset - Spontaneous Rate
BinEdge = [-0.05 0 0.01 0.03]; % Bin edges for Spike Counting w/ Histcounts
binDur = [BinEdge(2) - BinEdge(1), BinEdge(4) - BinEdge(3)]; % Bin Duration [Spontaneous, Evoked] -> Divisor when calculating Firing Rate later (Spikes/Period)

% --- Prep PSTH ---
dblBinDur = 5e-3; % Binsize of PSTH
vecTime = -0.1:dblBinDur:0.1; % PSTH X-Axis range/binsize
% indExcludeOn = vecTime > -2* dblBinDur & vecTime < 2*dblBinDur; % First millisecond before and after stim onset
% indExcludeOff = vecTime > (1+ -2* dblBinDur) & vecTime < (1+ 2*dblBinDur); % Ms before and after stim offset
% indExclude = [find(indExcludeOn) find(indExcludeOff)];
sOptions = -1;

for intCl = 1:intNumClu
    if ismember(sAP.sCluster(intCl).Area, vecROI) %selection criteria: Area, quality criteria?
        % Compute Zeta and Inst. Firing Rate of for pulseDur = 0.02
        vecSpikeTimes = sAP.sCluster(intCl).SpikeTimes;
        [dblZetaP,~,sRate,vecLatency] = zetatest(vecSpikeTimes,vecLaserOnSecs(pulseDurTrial == 0.02)-0.5,1,[],[],3);
        PeakLat_Cl = vecLatency(3);
        % --- Classify Clusters ---
        if (dblZetaP > 0.05) || (PeakLat_Cl < 0.501)
            continue;
        elseif (PeakLat_Cl < 0.51);
            % --- Get SpontRate_Cl for GAD2+ neuron ---
            sCounts_ER = zeros(numel(vecLaserOnSecs),1);
            for intTrial=1:structEP.intTrialNum
	            vecTrialEdges = BinEdge + vecLaserOnSecs(intTrial); % Add stim onset time for this trial to rel. bin edges to get absolute bin edges
	            [vecCounts,edges] = histcounts(vecSpikeTimes,vecTrialEdges);
	            sCounts_ER(intTrial) = vecCounts(1); % Counts Spontaneous Rate
            end
            
            sRate_ER = mean(sCounts_ER(:,1)/abs(binDur(1))); % Spike rates for each trial [Rate during -0.05s-0s, Rate during 0.01-0.03s]
            SpontRate_Cl = mean(sRate_ER(:,1));
            SpontRate_STD_Cl = std(sRate_ER(:,1));
            
            CellClass_Cl = "GAD2+";
        else
            % -- Calculate if sig difference between Spont Rate and FR (10-30 ms after optogenetic onset) ---
            sCounts_ER = zeros(numel(vecLaserOnSecs),2);
            for intTrial=1:structEP.intTrialNum
	            vecTrialEdges = BinEdge + vecLaserOnSecs(intTrial); % Add stim onset time for this trial to rel. bin edges to get absolute bin edges
	            [vecCounts,edges] = histcounts(vecSpikeTimes,vecTrialEdges);
	            sCounts_ER(intTrial,1) = vecCounts(1); % Counts Spontaneous Rate
                sCounts_ER(intTrial,2) = vecCounts(3); % Count for Response After Stimulation
            end
            
            sRate_ER = [sCounts_ER(:,1)/abs(binDur(1)) sCounts_ER(:,2)/abs(binDur(2))]; % Spike rates for each trial [Rate during -0.05s-0s, Rate during 0.01-0.03s]
            
            % Mean Spontaneous & Evoked Rates
            SpontRate_Cl = mean(sRate_ER(:,1));
            SpontRate_STD_Cl = std(sRate_ER(:,1));
            EvokedRate_20ms_Cl = mean(sRate_ER(pulseDurTrial == 0.02,2));
            
            % T-tests
            [~,p_val_20ms_Cl] = ttest(sRate_ER(pulseDurTrial == 0.02,2), sRate_ER(pulseDurTrial == 0.02,1), 'Alpha', 0.01);
            
            % Simple Classification for Now
            if (p_val_20ms_Cl < 0.01) && (EvokedRate_20ms_Cl > SpontRate_Cl)
                CellClass_Cl = "Activated";
            elseif (p_val_20ms_Cl < 0.01) && (EvokedRate_20ms_Cl < SpontRate_Cl);
                CellClass_Cl = "Inhibited";
            else
                CellClass_Cl = "Other";
            end
        end

        % --- PSTH ---
        [PSTHMean_20ms_Cl,PSTHSEM_20ms_Cl,PSTHBinCenters_20ms_Cl,~] = doPEP(vecSpikeTimes,vecTime,vecLaserOnSecs(pulseDurTrial == 0.02),sOptions);

        % --- Export Cluster Data ---

        ClusterN = [ClusterN; intCl];
        Area = [Area; string(sAP.sCluster(intCl).Area)];
        zeta_p_20ms = [zeta_p_20ms; dblZetaP];
        CellClass = [CellClass; CellClass_Cl];
        Peak_Lat = [Peak_Lat; PeakLat_Cl - 0.5];
        IFR_Rate = [IFR_Rate; {sRate.vecRate}];
        IFR_Time = [IFR_Time; {sRate.vecT}];
        SpontRate = [SpontRate; SpontRate_Cl];
        SpontRate_STD = [SpontRate_STD; SpontRate_STD_Cl];
        
        % PSTH Vals 20ms
        PSTHMean_20ms = [PSTHMean_20ms; PSTHMean_20ms_Cl];
        PSTHSEM_20ms = [PSTHSEM_20ms; PSTHSEM_20ms_Cl];
        PSTHBinCenters = [PSTHBinCenters; PSTHBinCenters_20ms_Cl];

    else
        continue;
    end
end

%% RecOverall

RecOverall.NCells = numel(ClusterN);
RecOverall.NCells_GAD2 = sum(CellClass == "GAD2+");
% RecOverall.PLatency_GAD2 = % look for function that gives values for boxplot
RecOverall.NCells_Inhibited = sum(CellClass == "Inhibited");
% RecOverall.PLatency_Inhibited = %
RecOverall.NCells_Activated = sum(CellClass == "Activated");
% RecOverall.PLatency_Activated = %
RecOverall.NCells_Other = sum(CellClass == "Other");

% Z-Scored PSTH Values
zscored_PSTH = ((PSTHMean_20ms - SpontRate)./SpontRate_STD);

RecData.Overall.PSTHMean_20ms_GAD2 = mean(zscored_PSTH(CellClass == "GAD2+",:), 1);
RecData.Overall.PSTHSEM_20ms_GAD2 = std(zscored_PSTH(CellClass == "GAD2+",:))/sum(CellClass == "GAD2+");

RecData.Overall.PSTHMEAN_20ms_Inhibited = mean(zscored_PSTH(CellClass == "Inhibited",:), 1);
RecData.Overall.PSTHSEM_20ms_Inhibited = std(zscored_PSTH(CellClass == "Inhibited",:))/sum(CellClass == "Inhibited");

RecData.Overall.PSTHMEAN_20ms_Activated = mean(zscored_PSTH(CellClass == "Activated",:), 1);
RecData.Overall.PSTHSEM_20ms_Activated = std(zscored_PSTH(CellClass == "Activated",:))/sum(CellClass == "Activated");

RecData.Overall.PSTHMEAN_20ms_Other = mean(zscored_PSTH(CellClass == "Other",:), 1);
RecData.Overall.PSTHSEM_20ms_Other = std(zscored_PSTH(CellClass == "Other",:))/sum(CellClass == "Other");

%Normalized PSTH Values
PSTHMean_20ms_BL = PSTHMean_20ms - SpontRate;
norm = max(PSTHMean_20ms_BL,[],2);

RecOverall.PSTHMean_20ms_GAD2_norm = mean(PSTHMean_20ms_BL(CellClass == "GAD2+",:)./norm(CellClass == "GAD2+"), 1);
RecOverall.PSTHSEM_20ms_GAD2_norm = std(PSTHMean_20ms_BL./norm(CellClass == "GAD2+"), 1)/sqrt(sum(CellClass == "GAD2+"));

RecOverall.PSTHMean_20ms_Inhibited = mean(PSTHMean_20ms_BL(CellClass == "Inhibited",:)./norm(CellClass == "Inhibited"), 1);
RecOverall.PSTHSEM_20ms_Inhibited = std(PSTHMean_20ms_BL./norm(CellClass == "Inhibited"), 1)/sqrt(sum(CellClass == "Inhibited"));

RecOverall.PSTHMean_20ms_Activated = mean(PSTHMean_20ms_BL(CellClass == "Activated",:)./norm(CellClass == "Activated"), 1);
RecOverall.PSTHSEM_20ms_Activated = std(PSTHMean_20ms_BL./norm(CellClass == "Activated"), 1)/sqrt(sum(CellClass == "Activated"));

RecOverall.PSTHMean_20ms_Other = mean(PSTHMean_20ms_BL(CellClass == "Other",:)./norm(CellClass == "Other"), 1);
RecOverall.PSTHSEM_20ms_Other = std(PSTHMean_20ms_BL./norm(CellClass == "Other"), 1)/sqrt(sum(CellClass == "Other"));

RecOverall.PSTHBinSize = dblBinDur; % Binsize of PSTH
RecOverall.PSTHtime = vecTime; % PSTH X-Axis range/binsize
RecOverall.PSTHBinCenters = PSTHBinCenters(1,:);

%% Write RecData to Overall Table

% Write Table (Subject)
RecData.ClusterData = table(ClusterN, CellClass, Area, zeta_p_20ms, ...
    Peak_Lat, IFR_Rate, IFR_Time, SpontRate, SpontRate_STD, PSTHMean_20ms, ...
    PSTHSEM_20ms, PSTHBinCenters);

% Overall Subject Data
RecData.OverallData = RecOverall;

% Add to DataOut
SubjectN = table(repmat(sAP.sJson.subject, [numel(ClusterN) 1]), 'VariableNames', {'Subject'});
RecDate = cell2table(repmat({sAP.sJson.date}', [numel(ClusterN) 1]), 'VariableNames', {'RecDate'});
DataOut_OT.AllMice.ClusterData = [DataOut_OT.AllMice.ClusterData; [SubjectN RecDate RecData.ClusterData]];

RecordingName = [replace(sAP.sJson.experiment(1:end-6),'-','_')];
DataOut_OT.(RecordingName) = RecData;

end

%%

%Create (additional) separate structs for optotagged and non-optotagged
%units!

AllM_Overall = struct;

n_clusters_overall = height(DataOut_OG.AllMice.ClusterData);
AllM_Overall.NMice = numel(unique(DataOut_OG.AllMice.ClusterData(:,1)));
AllM_Overall.NRecs = numel(fieldnames(DataOut_OG)) - 1;
AllM_Overall.NCells = n_clusters_overall; % Probably will need to tweak this!
AllM_Overall.NCells_Reduced = sum((DataOut_OG.AllMice.ClusterData.p_val < 0.01) & DataOut_OG.AllMice.ClusterData.PctChange < 0);
AllM_Overall.MeanPctReduction = mean(DataOut_OG.AllMice.ClusterData.PctChange(DataOut_OG.AllMice.ClusterData.PctChange < 0));
AllM_Overall.NCells_Increased = sum((DataOut_OG.AllMice.ClusterData.p_val < 0.01) & DataOut_OG.AllMice.ClusterData.PctChange > 0);
AllM_Overall.MeanPctIncrease = mean(DataOut_OG.AllMice.ClusterData.PctChange((DataOut_OG.AllMice.ClusterData.p_val<0.01) & (DataOut_OG.AllMice.ClusterData.PctChange > 0)));
AllM_Overall.ER_OptoOn = mean(DataOut_OG.AllMice.ClusterData.ER_OptoOn);
AllM_Overall.ER_OptoOff = mean(DataOut_OG.AllMice.ClusterData.ER_OptoOff);
AllM_Overall.SE_OptoOn = std(DataOut_OG.AllMice.ClusterData.ER_OptoOn, [], 1)/sqrt(n_clusters_overall);
AllM_Overall.SE_OptoOff = std(DataOut_OG.AllMice.ClusterData.ER_OptoOff, [], 1)/sqrt(n_clusters_overall);
[~, AllM_Overall.p_val] = ttest(DataOut_OG.AllMice.ClusterData.ER_OptoOn, DataOut_OG.AllMice.ClusterData.ER_OptoOff, 'Alpha', 0.01);

% Non-Normalized 
AllM_Overall.PSTHMean_Off = mean(DataOut_OG.AllMice.ClusterData.PSTHMean_Off, 1);
AllM_Overall.PSTHSEM_Off = std(DataOut_OG.AllMice.ClusterData.PSTHMean_Off, 1)/sqrt(n_clusters_overall);
AllM_Overall.PSTHMean_On = mean(DataOut_OG.AllMice.ClusterData.PSTHMean_On, 1);
AllM_Overall.PSTHSEM_On = std(DataOut_OG.AllMice.ClusterData.PSTHMean_On, 1)/sqrt(n_clusters_overall);

% Normalized Overall PSTH Values
norm = max(DataOut_OG.AllMice.ClusterData.PSTHMean_Off,[],2);
AllM_Overall.PSTHMean_Off_norm = mean(DataOut_OG.AllMice.ClusterData.PSTHMean_Off./norm, 1);
AllM_Overall.PSTHSEM_Off_norm = std(DataOut_OG.AllMice.ClusterData.PSTHMean_Off./norm, 1)/sqrt(n_clusters_overall);
AllM_Overall.PSTHMean_On_norm = mean(DataOut_OG.AllMice.ClusterData.PSTHMean_On./norm, 1);
AllM_Overall.PSTHSEM_On_norm = std(DataOut_OG.AllMice.ClusterData.PSTHMean_On./norm, 1)/sqrt(n_clusters_overall);

AllM_Overall.PSTHBinSize = 5e-3; % Binsize of PSTH
AllM_Overall.PSTHtime = vecTime; % PSTH X-Axis range/binsize
AllM_Overall.PSTHBinCenters = DataOut_OG.AllMice.ClusterData.PSTHBinCenters_On(1,:);

DataOut_OG.AllMice.ClusterDataUp = DataOut_OG.AllMice.ClusterData((DataOut_OG.AllMice.ClusterData.p_val < 0.01) & (DataOut_OG.AllMice.ClusterData.PctChange > 0),:);
DataOut_OG.AllMice.ClusterDataDown = DataOut_OG.AllMice.ClusterData((DataOut_OG.AllMice.ClusterData.p_val < 0.01) & (DataOut_OG.AllMice.ClusterData.PctChange < 0),:);
DataOut_OG.AllMice.ClusterDataNonSig = DataOut_OG.AllMice.ClusterData(~(DataOut_OG.AllMice.ClusterData.p_val < 0.01),:);

DataOut_OG.AllMice.Overall = AllM_Overall;
DataOut_OG.AllMice.LatencyOpto = mean(DiffTimes);
DataOut_OG.AllMice.LatencyOptoSEM = std(DiffTimes)/sqrt(length(DiffTimes));

%%

st = sAP.sCluster(488).SpikeTimes;

stim = vecLaserOnSecs(pulseDurTrial == 0.02) - 0.1;

maxDur = 0.2;

figure; hold on
plotRaster(st,stim,maxDur); xline(0, 'r'); fixfig