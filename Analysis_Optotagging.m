% Analysis Optotagging GAD2

intNumClu = length(sAP.sCluster);
structEP = sAP.cellStim{1,2}.structEP;
vecStimOnSecs = structEP.vecStimOnTime;
vecStimOffSecs = structEP.vecStimOffTime;

%% Prep Output Table

% MouseN, RecN, ClusterN, zeta_p, Instaneous FR, mean PSTH values (probably
% w/ 2 bin sizes?)
MouseN = [];
RecN = [];
ClusterN = [];


%% Compute Zeta, Get Latencies, and Get Average Responses

% Prep bin stuff for PSTH

for intCl = 1:intNumClu
    if (sAP.sCluster(intCl).Area == "Superior colliculus superificial grey layer") %selection criteria: Area, Visual responsiveness?, quality criteria?
    end
    vecSpikeTimes = 
end

%%

st = sSynthData.sCluster(488).SpikeTimes;

stim = sSynthData.cellStim{1,3}.structEP.vecStimOnTime - 0.1;
stim = stim(sSynthData.cellStim{1,3}.structEP.PulseDurTrial == 0.02);

maxDur = 0.2;

figure; hold on
plotRaster(st,stim,maxDur); xline(0.1, 'r'); fixfig