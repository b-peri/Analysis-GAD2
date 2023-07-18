function dblRateSpontaneous = computeRateSpontaneous(vecSpikeTimes,vecStimOnSecs,vecStimOffSecs,sParams)
%compute unit's spontaneous/baseline rate
intCount = 0;
dblPeriod = 0;
for intTrial=1:length(vecStimOffSecs)-1
    intCount = intCount + ...
        length(find(vecSpikeTimes>(vecStimOffSecs(intTrial)+sParams.dblSecsFromPrevStimOff) & ...
        vecSpikeTimes<vecStimOnSecs(intTrial+1)));
    dblPeriod = dblPeriod + vecStimOnSecs(intTrial+1) - (vecStimOffSecs(intTrial)+sParams.dblSecsFromPrevStimOff);
end
dblRateSpontaneous = intCount / dblPeriod;
if dblPeriod<5
    fprintf('Less than 5s to compute spontaneous rate.\n')
end
if intCount<10
    fprintf('Less than 10 spikes to compute spontaneous rate.\n')
end
end