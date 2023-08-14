% Analysis Optotagging GAD2

% Note for resp dynamics: We want to see 1. Optotagged GAD2 units, 2.
% Inhibited units, 3. (Putative) disinhibited units! (and maybe 4.
% Other/visually responsive units?)

%% Load in Data

% Load in Correct Data File
[FileNames, PathName] = uigetfile('*.mat', 'MultiSelect', 'on');

%% Set Up Struct

DataOut_OT.Overall = struct;
DataOut_OT.ClusterData = cell2table(cell(0,15), 'VariableNames', ...
    {'Subject', 'RecDate', 'ClusterN', 'CellClass', 'Area', 'zeta_p', ...
    'Peak_Lat', 'sRate_Early', 'sRate_Late', 'SpontRate', 'SpontRate_STD', ...
    'PSTHMean', 'PSTHSEM', 'PSTHBinCenters', 'NonStationarity'});

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
zeta_p = [];
Peak_Lat = [];
sRate_Early = [];
sRate_Late = [];
SpontRate = [];
SpontRate_STD = [];

PSTHMean = [];
PSTHSEM = [];
PSTHBinCenters = [];

NonStationarity = [];

%% Compute Zeta, Get Latencies, and Get Average Responses

% vecROI = ["Superior colliculus zonal layer" "Superior colliculus" + ...
%     " superficial gray layer" "Superior colliculus optic layer"];

vecROI = ["Superior colliculus zonal layer" "Superior colliculus" + ...
    " superficial gray layer" "Superior colliculus optic layer" ...
    "Superior colliculus motor related intermediate gray layer"];

% --- Select Pulse Duration ---
pulseSelect = 0.01;

% --- Prep Overall Modulation (Baseline vs. Evoked) ---
% Spontaneous Rate = Mean rate in -1s - -100ms before stim onset
% Visually Evoked Rate = Mean rate in 10-250ms after stim onset - Spontaneous Rate
BinEdge = [-0.5 -0.005 0.002 0.01 0.03]; % Bin edges for Spike Counting w/ Histcounts
binDur = [BinEdge(2) - BinEdge(1), BinEdge(4) - BinEdge(3), BinEdge(5) - BinEdge(4)]; % Bin Duration [Spontaneous, Evoked] -> Divisor when calculating Firing Rate later (Spikes/Period)

% --- Prep PSTH ---
dblBinDur = 2e-3; % Binsize of PSTH
vecTime = -0.1:dblBinDur:0.1; % PSTH X-Axis range/binsize
sOptions.handleFig = -1;

for intCl = 1:intNumClu
    if ismember(sAP.sCluster(intCl).Area, vecROI) && sAP.sCluster(intCl).Violations1ms < 0.25 %&& abs(sAP.sCluster(intCl).NonStationarity) < 0.25 %selection criteria
        % Compute Zeta and Inst. Firing Rate of for pulseDur = 0.02
        vecSpikeTimes = sAP.sCluster(intCl).SpikeTimes;
        if size(vecSpikeTimes, 1) < 1000
            continue;
        end

        % --- Compute Zeta ---
        [dblZetaP,~,structRate] = zetatest(vecSpikeTimes,vecLaserOnSecs(pulseDurTrial == pulseSelect)-0.5,1);
        if dblZetaP >= 0.01
            continue;
        end

        % --- Find Peak in Inst. FR Plot (Outside of excl. window) ---
        peakWin_on = [0.5 + BinEdge(2) 0.5 + BinEdge(3)];
        peakWin_off = [0.499+pulseSelect 0.502+pulseSelect];
        vecRate_Filt = structRate.vecRate((structRate.vecT < peakWin_on(1) | structRate.vecT > peakWin_on(2)) ...
            & (structRate.vecT < peakWin_off(1) | structRate.vecT > peakWin_off(2)));
        [~,indPeak] = max(vecRate_Filt);
        vecT_short = structRate.vecT((structRate.vecT < peakWin_on(1) | structRate.vecT > peakWin_on(2)) ...
            & (structRate.vecT < peakWin_off(1) | structRate.vecT > peakWin_off(2)));
        PeakLat_Cl = vecT_short(indPeak); 

        RFmap_times = [sAP.cellBlock{2}.vecStimOnTime(1) sAP.cellBlock{2}.vecStimOffTime(end)];
        spikesRF = sum(vecSpikeTimes > RFmap_times(1) & vecSpikeTimes < RFmap_times(2));

        if isempty(PeakLat_Cl) || spikesRF == 0 % Exclusion: If no spikes during RF mapper or if no spikes outside of artefact period
            continue;
        end

        % --- Get SpikeRates ---

        trial_sel = vecLaserOnSecs(pulseDurTrial == pulseSelect);
        sCounts = zeros(numel(trial_sel),3);
        for intTrial=1:numel(trial_sel)
            vecTrialEdges = BinEdge + trial_sel(intTrial); % Add stim onset time for this trial to rel. bin edges to get absolute bin edges
            intLaserOn = trial_sel(intTrial);
            vecSpikeShort = vecSpikeTimes(vecSpikeTimes > vecTrialEdges(1) & vecSpikeTimes < vecTrialEdges(end));
            vecSpikeShort(vecSpikeShort > (intLaserOn+pulseSelect - 0.001) & vecSpikeShort < (intLaserOn + pulseSelect + 0.002)) = []; %Remove spikes around pulse offset

            [vecCounts,edges] = histcounts(vecSpikeShort,vecTrialEdges);

            sCounts(intTrial,1) = vecCounts(1); % Counts Spontaneous Rate
            sCounts(intTrial,2) = vecCounts(3); % Counts spikes during first 10ms
            sCounts(intTrial,3) = vecCounts(4); % Counts spikes during 10-30ms
        end

        sRate = [sCounts(:,1)/abs(binDur(1)), sCounts(:,2)/abs(binDur(2)), sCounts(:,3)/abs(binDur(3))]; % Spike rates for each trial [Rate during -0.05s-0s, Rate during 0.01-0.03s]
        
        SpontRate_Cl = mean(sRate(:,1));
        SpontRate_STD_Cl = std(sRate(:,1));
        sRate_Early_Cl = mean(sRate(:,2));
        sRate_Late_Cl = mean(sRate(:,3));

        % --- Classify Clusters ---
        if (ttest(sRate(:,1),sRate(:,2), "Alpha", 0.01) == 1) && (sRate_Early_Cl > SpontRate_Cl) && (PeakLat_Cl > 0.5) && (PeakLat_Cl <= 0.510) && (sRate_Early_Cl >= sRate_Late_Cl)
            CellClass_Cl = "GAD2+";
        elseif (ttest(sRate(:,1), sRate(:,3), "Alpha", 0.01) == 1)
            if sRate_Late_Cl > SpontRate_Cl
                CellClass_Cl = "Activated";
            elseif sRate_Late_Cl < SpontRate_Cl
                CellClass_Cl = "Inhibited";
            end
        else
            CellClass_Cl = "Other";
        end

        % --- PSTH ---
        [PSTHMean_Cl,PSTHSEM_Cl,PSTHBinCenters_Cl,~] = doPEP({vecSpikeTimes},vecTime,vecLaserOnSecs(pulseDurTrial == pulseSelect),sOptions);

        % --- Export Cluster Data ---

        ClusterN = [ClusterN; intCl];
        Area = [Area; string(sAP.sCluster(intCl).Area)];
        zeta_p = [zeta_p; dblZetaP];
        CellClass = [CellClass; CellClass_Cl];
        Peak_Lat = [Peak_Lat; PeakLat_Cl - 0.5];
        sRate_Early = [sRate_Early; sRate_Early_Cl];
        sRate_Late = [sRate_Late; sRate_Late_Cl];
        SpontRate = [SpontRate; SpontRate_Cl];
        SpontRate_STD = [SpontRate_STD; SpontRate_STD_Cl];
        
        % PSTH Vals
        PSTHMean = [PSTHMean; PSTHMean_Cl];
        PSTHSEM = [PSTHSEM; PSTHSEM_Cl];
        PSTHBinCenters = [PSTHBinCenters; PSTHBinCenters_Cl];

        NonStationarity = [NonStationarity; sAP.sCluster(intCl).NonStationarity];
    else
        continue;
    end
end

if isempty(ClusterN); % If no cells found, skip rest of analysis
    continue;
end

%% Write RecData to Overall Table

% Add to DataOut
SubjectN = table(repmat(string(sAP.sJson.subject), [numel(ClusterN) 1]), 'VariableNames', {'Subject'});
RecDate = table(repmat(string(sAP.sJson.date), [numel(ClusterN) 1]), 'VariableNames', {'RecDate'});
DataOut_OT.ClusterData = [DataOut_OT.ClusterData; [SubjectN ...
    RecDate table(ClusterN, CellClass, Area, zeta_p, ...
    Peak_Lat, sRate_Early, sRate_Late, SpontRate, SpontRate_STD, PSTHMean, ...
    PSTHSEM, PSTHBinCenters, NonStationarity);]];

end

%% Overall Data

DOT = DataOut_OT.ClusterData;

n_clusters_overall = height(DOT);
DataOut_OT.Overall.NMice = numel(unique(DOT.Subject));
DataOut_OT.Overall.NRecs = numel(unique([DataOut_OT.ClusterData.Subject DataOut_OT.ClusterData.RecDate], 'rows'));
DataOut_OT.Overall.NCells = n_clusters_overall; % Probably will need to tweak this!
DataOut_OT.Overall.NCells_GAD2 = sum(DOT.CellClass == "GAD2+");
DataOut_OT.Overall.NCells_Inh = sum(DOT.CellClass == "Inhibited");
DataOut_OT.Overall.NCells_Act = sum(DOT.CellClass == "Activated");
DataOut_OT.Overall.NCells_Oth = sum(DOT.CellClass == "Other");

% Z-Scored PSTH Values
zscored_PSTH = ((DOT.PSTHMean - DOT.SpontRate)./DOT.SpontRate_STD);

DataOut_OT.Overall.PSTHMean_GAD2_z = mean(zscored_PSTH(DOT.CellClass == "GAD2+",:), 1);
DataOut_OT.Overall.PSTHSEM_GAD2_z = std(zscored_PSTH(DOT.CellClass == "GAD2+",:))/sum(DOT.CellClass == "GAD2+");

DataOut_OT.Overall.PSTHMean_Inh_z = mean(zscored_PSTH(DOT.CellClass == "Inhibited",:), 1);
DataOut_OT.Overall.PSTHSEM_Inh_z = std(zscored_PSTH(DOT.CellClass == "Inhibited",:))/sum(DOT.CellClass == "Inhibited");

DataOut_OT.Overall.PSTHMean_Act_z = mean(zscored_PSTH(DOT.CellClass == "Activated",:), 1);
DataOut_OT.Overall.PSTHSEM_Act_z = std(zscored_PSTH(DOT.CellClass == "Activated",:))/sum(DOT.CellClass == "Activated");

DataOut_OT.Overall.PSTHMean_Oth_z = mean(zscored_PSTH(DOT.CellClass == "Other",:), 1);
DataOut_OT.Overall.PSTHSEM_Oth_z = std(zscored_PSTH(DOT.CellClass == "Other",:))/sum(DOT.CellClass == "Other");

%Normalized PSTH Values
PSTHMean_BL = DOT.PSTHMean - DOT.SpontRate;
norm = max(PSTHMean_BL,[],2);

DataOut_OT.Overall.PSTHMean_GAD2_norm = mean(PSTHMean_BL(DOT.CellClass == "GAD2+",:)./norm(DOT.CellClass == "GAD2+"), 1);
DataOut_OT.Overall.PSTHSEM_GAD2_norm = std(PSTHMean_BL(DOT.CellClass == "GAD2+",:)./norm(DOT.CellClass == "GAD2+"), 1)/sqrt(sum(DOT.CellClass == "GAD2+"));

DataOut_OT.Overall.PSTHMean_Inh_norm = mean(PSTHMean_BL(DOT.CellClass == "Inhibited",:)./norm(DOT.CellClass == "Inhibited"), 1);
DataOut_OT.Overall.PSTHSEM_Inh_norm = std(PSTHMean_BL(DOT.CellClass == "Inhibited",:)./norm(DOT.CellClass == "Inhibited"), 1)/sqrt(sum(DOT.CellClass == "Inhibited"));

DataOut_OT.Overall.PSTHMean_Act_norm = mean(PSTHMean_BL(DOT.CellClass == "Activated",:)./norm(DOT.CellClass == "Activated"), 1);
DataOut_OT.Overall.PSTHSEM_Act_norm = std(PSTHMean_BL(DOT.CellClass == "Activated",:)./norm(DOT.CellClass == "Activated"), 1)/sqrt(sum(DOT.CellClass == "Activated"));

DataOut_OT.Overall.PSTHMean_Oth_norm = mean(PSTHMean_BL(DOT.CellClass == "Other",:)./norm(DOT.CellClass == "Other"), 1);
DataOut_OT.Overall.PSTHSEM_Oth_norm = std(PSTHMean_BL(DOT.CellClass == "Other",:)./norm(DOT.CellClass == "Other"), 1)/sqrt(sum(DOT.CellClass == "Other"));

DataOut_OT.Overall.PSTHBinSize = dblBinDur; % Binsize of PSTH
DataOut_OT.Overall.PSTHtime = vecTime; % PSTH X-Axis range/binsize
DataOut_OT.Overall.PSTHBinCenters = PSTHBinCenters(1,:);

%% Save Output

if all(startsWith(string(DataOut_OT.ClusterData.Subject),'7'))
    if pulseSelect == 0.05
        SaveFile = ['DataOut_Optotagging50ms_' datestr(datetime("today"),"dd-mm-yy") '.mat'];
        DataOut_OT_50ms = DataOut_OT;
        save(SaveFile, 'DataOut_OT_50ms');
    elseif pulseSelect == 0.02
        SaveFile = ['DataOut_Optotagging20ms_' datestr(datetime("today"),"dd-mm-yy") '.mat'];
        DataOut_OT_20ms = DataOut_OT;
        save(SaveFile, 'DataOut_OT_20ms');
    elseif pulseSelect == 0.01
        SaveFile = ['DataOut_Optotagging10ms_' datestr(datetime("today"),"dd-mm-yy") '.mat'];
        DataOut_OT_10ms = DataOut_OT;
        save(SaveFile, 'DataOut_OT_10ms');
    elseif pulseSelect == 0.005
        SaveFile = ['DataOut_Optotagging5ms_' datestr(datetime("today"),"dd-mm-yy") '.mat'];
        DataOut_OT_5ms = DataOut_OT;
        save(SaveFile, 'DataOut_OT_5ms');
    elseif pulseSelect == 0.002
        SaveFile = ['DataOut_Optotagging2ms_' datestr(datetime("today"),"dd-mm-yy") '.mat'];
        DataOut_OT_2ms = DataOut_OT;
        save(SaveFile, 'DataOut_OT_2ms');
    end
end

clearvars -except DataOut_OT_50ms DataOut_OT_20ms DataOut_OT_10ms DataOut_OG_GAD2 DataOut_OG_Cont DataOut_OGt