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
pulseDurs = structEP.PulseDurTrial
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

vecROI = ["Superior colliculus zonal layer" "Superior colliculus" + ...
    " superficial gray layer" "Superior colliculus optic layer" ...
    "Superior colliculus motor related intermediate gray layer" ...
    "Superior colliculus motor related intermediate white layer"];

% Prep bin stuff for PSTH

for intCl = 1:intNumClu
    if ismember(sAP.sCluster(intCl).Area, vecROI) %selection criteria: Area, quality criteria?
        vecSpikeTimes = sAP.sCluster(intCl).SpikeTimes;
        [dblZetaP,~,sRate,vecLatencies] = zetatest(vecSpikeTimes,vecLaserOnSecs(,0.5);

        ClusterN = [ClusterN; intCl];
        Area = [Area; string(sAP.sCluster(intCl).Area)];
        zeta_p = [zeta_p; dblZetaP];
        % Peak_Lat = [];
        % Inst_FR = [];
        % PSTHMean = [];
        % PSTHSEM = [];
        % PSTHBinCenters = [];
    else
        continue;
    end
end

%% Output RecData

%END LOOP
%%

%Create (additional) separate structs for optotagged and non-optotagged
%units!

%%

st = sSynthData.sCluster(488).SpikeTimes;

stim = sSynthData.cellStim{1,3}.structEP.vecStimOnTime - 0.1;
stim = stim(sSynthData.cellStim{1,3}.structEP.PulseDurTrial == 0.02);

maxDur = 0.2;

figure; hold on
plotRaster(st,stim,maxDur); xline(0.1, 'r'); fixfig