% Analysis Optotagging GAD2

% Note for resp dynamics: We want to see 1. Optotagged GAD2 units, 2.
% Inhibited units, 3. Potential disinhibited units!

%Prep Table
DataOut_OT.AllMice.ClusterData = cell2table(cell(0,13), 'VariableNames', ...
    {'Subject', 'RecDate', 'ClusterN', 'Area', 'zeta_p', 'Peak_lat', ...
    'Inst_FR', 'PSTHMean', 'PSTHSEM_Off', 'PSTHBinCenters_Off', ...
    'PSTHMean_On', 'PSTHSEM_On','PSTHBinCenters_On'});

%START LOOP OVER RECS
%% Grab Stimulus Data

intNumClu = length(sAP.sCluster);
structEP = sAP.cellBlock{1,3};
vecLaserOnSecs = structEP.vecLaserOnTime;
pulseDurTrial = structEP.PulseDurTrial;
pulseDurs = structEP.vecPulseDuration;
% vecStimOffSecs = structEP.vecLaserOffTime;

%% Prep Output Table

% Initial Data RecData
RecData.Subject = sAP.sJson.subject;
RecData.SubjectType = sAP.sJson.subjecttype;

% MouseN, RecN, ClusterN, zeta_p, peak latency, Instaneous FR, mean PSTH values (probably
% w/ 2 bin sizes?)
MouseN = sAP.sJson.subject;
RecN = sAP.sJson.date;
ClusterN = [];
Area = [];
zeta_p = [];
Peak_Lat = [];
Inst_FR = [];
PSTHMean = [];
PSTHSEM = [];
PSTHBinCenters = [];

%% Compute Zeta, Get Latencies, and Get Average Responses

% vecROI = ["Superior colliculus zonal layer" "Superior colliculus" + ...
%     " superficial gray layer" "Superior colliculus optic layer"];

vecROI = ["Superior colliculus zonal layer" "Superior colliculus" + ...
    " superficial gray layer" "Superior colliculus optic layer" ...
    "Superior colliculus motor related intermediate gray layer" ...
    "Superior colliculus motor related intermediate white layer"];

% --- Select Pulse Duration ---
pulseDur = 0.02;
n_trials_dur = sum(pulseDurTrial == pulseDur);
% vecLaserOnSecs_dur = vecLaserOnSecs(pulseDurTrial == pulseDurTrial);

% --- Prep Overall Modulation (Baseline vs. Evoked) ---
% Spontaneous Rate = Mean rate in -1s - -100ms before stim onset
% Visually Evoked Rate = Mean rate in 10-250ms after stim onset - Spontaneous Rate
BinEdge = [-0.05 0 0.01 0.03]; % Bin edges for Spike Counting w/ Histcounts
binDur = [BinEdge(2) - BinEdge(1), BinEdge(4) - BinEdge(3)]; % Bin Duration [Spontaneous, Evoked] -> Divisor when calculating Firing Rate later (Spikes/Period)

% --- Prep PSTH ---
dblBinDur = 5e-3; % Binsize of PSTH
vecTime = -0.1:dblBinDur:0.1; % PSTH X-Axis range/binsize
indExcludeOn = vecTime > -2* dblBinDur & vecTime < 2*dblBinDur; % First millisecond before and after stim onset
indExcludeOff = vecTime > (1+ -2* dblBinDur) & vecTime < (1+ 2*dblBinDur); % Ms before and after stim offset
indExclude = [find(indExcludeOn) find(indExcludeOff)];
sOptions = -1;

for intCl = 1:intNumClu
    if ismember(sAP.sCluster(intCl).Area, vecROI) %selection criteria: Area, quality criteria?
        % Compute Zeta and Inst. Firing Rate
        vecSpikeTimes = sAP.sCluster(intCl).SpikeTimes;
        [dblZetaP,~,sRate,vecLatency] = zetatest(vecSpikeTimes,vecLaserOnSecs(pulseDurTrial == 0.01)-0.5,2,[],4,3);
        % dblZetaP = zetatest(vecSpikeTimes,vecLaserOnSecs,0.5); % -> 1. Should I be tested for Zeta over all stimDurs? 2. Is this window large enough to capture visual response?
        
        %if peak latency within 1-10ms -> optotagged; else continue with Sig Difference btwn spont rate and fr  
        % Check If Unit is optotagged
        % if ~(peakLatency < 0.001) && (peakLatency < 0.01);
            % CellResp = "GAD2+";
        % else
        % -- Analysis pt 2. Calculate if sig difference between Spont Rate and FR
            % (10-30 ms after optogenetic onset) ---
    
            sCounts_ER = zeros(numel(vecLaserOnSecs),2);
            for intTrial=1:n_trials_dur
	            vecTrialEdges = BinEdge + vecLaserOnSecs_dur(intTrial); % Add stim onset time for this trial to rel. bin edges to get absolute bin edges
	            [vecCounts,edges] = histcounts(vecSpikeTimes,vecTrialEdges);
	            sCounts_ER(intTrial,1) = vecCounts(1); % Counts Spontaneous Rate
                sCounts_ER(intTrial,2) = vecCounts(3); % Count for Visual Response
            end
            
            sRate_ER = [sCounts_ER(:,1)/binDur(1) sCounts_ER(:,2)/binDur(2)]; % Spike rates for each trial [Rate during -0.05s-0s, Rate during 0.01-0.03s]
            
            % Spontaneous Rate Overall
            SpontRate_Cl = mean(sRate_ER(:,1));
        
            % Mean Evoked FRs (Cluster)
            EvokedRate_Overall_Cl = mean(EvokedRate_Trial);
            % EvokedRate_2ms_Cl = mean(sRate_ER(pulseDurTrial == 0.002,2));
            % EvokedRate_5ms_Cl = mean(sRate_ER(pulseDurTrial == 0.005,2));
            % EvokedRate_10ms_Cl = mean(sRate_ER(pulseDurTrial == 0.01,2));
            % EvokedRate_20ms_Cl = mean(sRate_ER(pulseDurTrial == 0.02,2));
            % EvokedRate_50ms_Cl = mean(sRate_ER(pulseDurTrial == 0.05,2));
            
            % T-tests -> P-threshold?
            [~,p_val_Cl] = ttest(sRate_ER, sRate_ER(:,1), 'Alpha', 0.01);
            [~,p_val_2ms_Cl] = ttest(EvokedRate_2ms_Cl, sRate_ER(pulseDurTrial == 0.002,1), 'Alpha', 0.01);
            [~,p_val_5ms_Cl] = ttest(EvokedRate_5ms_Cl, sRate_ER(pulseDurTrial == 0.005,1), 'Alpha', 0.01);
            [~,p_val_10ms_Cl] = ttest(EvokedRate_10ms_Cl, sRate_ER(pulseDurTrial == 0.01,1), 'Alpha', 0.01);
            [~,p_val_20ms_Cl] = ttest(EvokedRate_20ms_Cl, sRate_ER(pulseDurTrial == 0.02,1), 'Alpha', 0.01);
            [~,p_val_50ms_Cl] = ttest(EvokedRate_50ms_Cl, sRate_ER(pulseDurTrial == 0.05,1), 'Alpha', 0.01);
    
            % Simple Classification for Now
            if (p_val_Cl) EvokedRate_Overall_Cl > SpontRate_Cl;
                CellResp_Cl = "Activated";
            elseif EvokedRate_Overall_Cl =< SpontRate_Cl;
                CellResp_Cl = "Inhibited";
            else
                CellResp_Cl = "Other";
        % end

        % --- Analysis pt 3. PSTH ---

        [PSTHMean_Cl,PSTHSEM_Cl,PSTHBinCenters_Cl,~] = doPEP(vecSpikeTimes,vecTime,vecLaserOnSecs,sOptions);
        [PSTHMean_2ms_Cl,PSTHSEM_2ms_Cl,PSTHBinCenters_2ms_Cl,~] = doPEP(vecSpikes,vecTime,vecLaserOnSecs(pulseDurTrial == 0.002),sOptions);
        [PSTHMean_5ms_Cl,PSTHSEM_5ms_Cl,PSTHBinCenters_5ms_Cl,~] = doPEP(vecSpikes,vecTime,vecLaserOnSecs(pulseDurTrial == 0.005),sOptions);
        [PSTHMean_10ms_Cl,PSTHSEM_10ms_Cl,PSTHBinCenters_10ms_Cl,~] = doPEP(vecSpikes,vecTime,vecLaserOnSecs(pulseDurTrial == 0.01),sOptions);
        [PSTHMean_20ms_Cl,PSTHSEM_20ms_Cl,PSTHBinCenters_20ms_Cl,~] = doPEP(vecSpikes,vecTime,vecLaserOnSecs(pulseDurTrial == 0.02),sOptions);
        [PSTHMean_50ms_Cl,PSTHSEM_50ms_Cl,PSTHBinCenters_50ms_Cl,~] = doPEP(vecSpikes,vecTime,vecLaserOnSecs(pulseDurTrial == 0.05),sOptions);

        % --- Export Cluster Data ---

        ClusterN = [ClusterN; intCl];
        Area = [Area; string(sAP.sCluster(intCl).Area)];
        zeta_p = [zeta_p; dblZetaP];
        CellResp = [CellResp; CellResp_Cl];
        % Peak_Lat = [];
        % Inst_FR = [];

        %PSTH VALS -> MAKE BASELINE-CORRECTED
        
        % PSTH Vals Overall
        PSTHMean = [PSTHMean; PSTHMean_Cl - SpontRate_Cl];
        PSTHSEM = [PSTHSEM; PSTHSEM_Cl];
        PSTHBinCenters = [PSTHBinCenters; PSTHBinCenters_Cl];
        % PSTH Vals 2ms
        PSTHMean_2ms = [PSTHMean_2ms; PSTHMean_2ms_Cl - SpontRate_Cl];
        PSTHSEM_2ms = [PSTHSEM_2ms; PSTHSEM_2ms_Cl];
        PSTHBinCenters_2ms = [PSTHBinCenters_2ms; PSTHBinCenters_2ms_Cl];
        % PSTH Vals 5ms
        PSTHMean_5ms = [PSTHMean_5ms; PSTHMean_5ms_Cl - SpontRate_Cl];
        PSTHSEM_5ms = [PSTHSEM_5ms; PSTHSEM_5ms_Cl];
        PSTHBinCenters_5ms = [PSTHBinCenters_5ms; PSTHBinCenters_5ms_Cl];
        % PSTH Vals 10ms
        PSTHMean_10ms = [PSTHMean_10ms; PSTHMean_10ms_Cl - SpontRate_Cl];
        PSTHSEM_10ms = [PSTHSEM_10ms; PSTHSEM_10ms_Cl];
        PSTHBinCenters_10ms = [PSTHBinCenters_10ms; PSTHBinCenters_10ms_Cl];
        % PSTH Vals 20ms
        PSTHMean_20ms = [PSTHMean_20ms; PSTHMean_20ms_Cl - SpontRate_Cl];
        PSTHSEM_20ms = [PSTHSEM_20ms; PSTHSEM_20ms_Cl];
        PSTHBinCenters_20ms = [PSTHBinCenters_20ms; PSTHBinCenters_20ms_Cl];
        % PSTH Vals 50ms
        PSTHMean_50ms = [PSTHMean_50ms; PSTHMean_50ms_Cl - SpontRate_Cl];
        PSTHSEM_50ms = [PSTHSEM_50ms; PSTHSEM_50ms_Cl];
        PSTHBinCenters_50ms = [PSTHBinCenters_50ms; PSTHBinCenters_50ms_Cl];
    else
        continue;
    end
end

%% RecOverall

%mean response PSTH 20ms

RecData.Overall.PSTHMean_20ms_GAD2 = mean(PSTHMean_20ms(CellResp == "GAD2+",:),2);
RecData.Overall.PSTHSEM_20ms_GAD2 = mean

%% Output RecData

%END LOOP
%%

%Create (additional) separate structs for optotagged and non-optotagged
%units!

%%

st = sAP.sCluster(488).SpikeTimes;

stim = vecLaserOnSecs(pulseDurTrial == 0.02) - 0.1;

maxDur = 0.2;

figure; hold on
plotRaster(st,stim,maxDur); xline(0, 'r'); fixfig