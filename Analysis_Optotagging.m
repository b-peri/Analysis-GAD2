% Analysis Optotagging GAD2

% Note for resp dynamics: We want to see 1. Optotagged GAD2 units, 2.
% Inhibited units, 3. (Putative) disinhibited units! (and maybe 4.
% Other/visually responsive units?)

% To-Do

%%

%Prep Table
DataOut_OT.AllMice.ClusterData = cell2table(cell(0,14), 'VariableNames', ...
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
CellClass = [];
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
        [dblZetaP,~,sRate] = zetatest(vecSpikeTimes,vecLaserOnSecs(pulseDurTrial == 0.02)-0.5,1);
        if dblZetaP >= 0.01
            continue;
        end
        [~,indPeak] = max(sRate.vecRate);
        PeakLat_Cl = sRate.vecT(indPeak);
        % --- Classify Clusters ---
        if (PeakLat_Cl < 0.501)
            continue;
        elseif (PeakLat_Cl < 0.51)
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
DOT = DataOut_OT.AllMice.ClusterData;
All_M = struct;

n_clusters_overall = height(DOT);
All_M.NMice = numel(unique(DOT.Subject));
All_M.NRecs = numel(fieldnames(DataOut_OT)) - 1;
All_M.NCells = n_clusters_overall; % Probably will need to tweak this!
All_M.NCells_GAD2 = sum(DOT.CellClass == "GAD2+");
All_M.NCells_Inh = sum(DOT.CellClass == "Inhibited");
All_M.NCells_Act = sum(DOT.CellClass == "Activated");
All_M.NCells_Oth = sum(DOT.CellClass == "Other");

% Z-Scored PSTH Values
zscored_PSTH = ((DOT.PSTHMean_20ms - DOT.SpontRate)./DOT.SpontRate_STD);

All_M.PSTHMean_20ms_GAD2 = mean(zscored_PSTH(DOT.CellClass == "GAD2+",:), 1);
All_M.PSTHSEM_20ms_GAD2 = std(zscored_PSTH(DOT.CellClass == "GAD2+",:))/sum(DOT.CellClass == "GAD2+");

All_M.PSTHMean_20ms_Inh = mean(zscored_PSTH(DOT.CellClass == "Inhibited",:), 1);
All_M.PSTHSEM_20ms_Inh = std(zscored_PSTH(DOT.CellClass == "Inhibited",:))/sum(DOT.CellClass == "Inhibited");

All_M.PSTHMean_20ms_Act = mean(zscored_PSTH(DOT.CellClass == "Activated",:), 1);
All_M.PSTHSEM_20ms_Act = std(zscored_PSTH(DOT.CellClass == "Activated",:))/sum(DOT.CellClass == "Activated");

All_M.PSTHMean_20ms_Oth = mean(zscored_PSTH(DOT.CellClass == "Other",:), 1);
All_M.PSTHSEM_20ms_Oth = std(zscored_PSTH(DOT.CellClass == "Other",:))/sum(DOT.CellClass == "Other");

%Normalized PSTH Values
PSTHMean_20ms_BL = DOT.PSTHMean_20ms - DOT.SpontRate;
norm = max(PSTHMean_20ms_BL,[],2);

All_M.PSTHMean_20ms_GAD2_norm = mean(PSTHMean_20ms_BL(DOT.CellClass == "GAD2+",:)./norm(DOT.CellClass == "GAD2+"), 1);
All_M.PSTHSEM_20ms_GAD2_norm = std(PSTHMean_20ms_BL(DOT.CellClass == "GAD2+",:)./norm(DOT.CellClass == "GAD2+"), 1)/sqrt(sum(DOT.CellClass == "GAD2+"));

All_M.PSTHMean_20ms_Inh_norm = mean(PSTHMean_20ms_BL(DOT.CellClass == "Inhibited",:)./norm(DOT.CellClass == "Inhibited"), 1);
All_M.PSTHSEM_20ms_Inh_norm = std(PSTHMean_20ms_BL(DOT.CellClass == "Inhibited",:)./norm(DOT.CellClass == "Inhibited"), 1)/sqrt(sum(DOT.CellClass == "Inhibited"));

All_M.PSTHMean_20ms_Act_norm = mean(PSTHMean_20ms_BL(DOT.CellClass == "Activated",:)./norm(DOT.CellClass == "Activated"), 1);
All_M.PSTHSEM_20ms_Act_norm = std(PSTHMean_20ms_BL(DOT.CellClass == "Activated",:)./norm(DOT.CellClass == "Activated"), 1)/sqrt(sum(DOT.CellClass == "Activated"));

All_M.PSTHMean_20ms_Oth_norm = mean(PSTHMean_20ms_BL(DOT.CellClass == "Other",:)./norm(DOT.CellClass == "Other"), 1);
All_M.PSTHSEM_20ms_Oth_norm = std(PSTHMean_20ms_BL(DOT.CellClass == "Other",:)./norm(DOT.CellClass == "Other"), 1)/sqrt(sum(DOT.CellClass == "Other"));

All_M.PSTHBinSize = dblBinDur; % Binsize of PSTH
All_M.PSTHtime = vecTime; % PSTH X-Axis range/binsize
All_M.PSTHBinCenters = PSTHBinCenters(1,:);

DataOut_OT.AllMice.Overall = All_M;

%%

st = sAP.sCluster(488).SpikeTimes;

stim = vecLaserOnSecs(pulseDurTrial == 0.02) - 0.1;

maxDur = 0.2;

figure; hold on
plotRaster(st,stim,maxDur); xline(0, 'r'); fixfig