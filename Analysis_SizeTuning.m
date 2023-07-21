%% Size Tuning

% Get Stimulus Onset Times

sAP = sSynthData;
intNumClu = length(sAP.sCluster);
structEP = sAP.cellStim{1,2}.structEP;  
vecStimOnSecs = structEP.vecStimOnTime;
vecStimOffSecs = structEP.vecStimOffTime;