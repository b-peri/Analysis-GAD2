
%GratingPatchesAnalysis

%% get stimulus info
sAP = sSynthData;
intNumClu = length(sAP.sCluster);
structEP = sAP.cellStim{1,2}.structEP;  
vecStimOnSecs = structEP.vecStimOnTime;
vecStimOffSecs = structEP.vecStimOffTime;

%% get grid data
vecUniqueRects = unique(structEP.vecDstRect','rows'); %unique dst rects
vecUniqueStims = 1:length(vecUniqueRects);
vecStimIdx = zeros(size(structEP.vecDstRect,2),1);
for intStim = 1:length(vecUniqueRects)
    vecStimIdx(ismember(structEP.vecDstRect',vecUniqueRects(intStim,:),'rows')) = vecUniqueStims(intStim);
end
vecX_pix = unique(vecUniqueRects(:,1))+(vecUniqueRects(1,3)-unique(vecUniqueRects(1,1)))/2;
vecY_pix = unique(vecUniqueRects(:,2))+(vecUniqueRects(1,4)-unique(vecUniqueRects(1,2)))/2;

%% loop through data
matAvgRespAll = NaN(numel(vecY_pix),numel(vecX_pix),intNumClu);
for intClu = 1:intNumClu
    vecSpikes = sAP.sCluster(intClu).SpikeTimes; 
    vecRate = zeros(1,structEP.intTrialNum);
    for intTrial = 1:structEP.intTrialNum
        vecSpikeT = vecSpikes(vecSpikes>vecStimOnSecs(intTrial)&vecSpikes<vecStimOffSecs(intTrial));
        vecRate(intTrial) = numel(vecSpikeT)/(vecStimOffSecs(intTrial)-vecStimOnSecs(intTrial));
    end
    matAvgResp = NaN(numel(vecY_pix),numel(vecX_pix));

    for intLoc = vecUniqueStims
        matAvgResp(intLoc) = mean(vecRate(vecStimIdx==intLoc));
    end
    sParams.dblSecsFromPrevStimOff = 0.1; %s, for computing unit's baseline rate
    dblRateSpontaneous = computeRateSpontaneous(vecSpikeT,vecStimOnSecs,vecStimOffSecs,sParams);
    matAvgRespAll(:,:,intClu) = matAvgResp-dblRateSpontaneous;
end

%% plot data
%interpolate
vecX_pix_interp = linspace(vecX_pix(1),vecX_pix(end),16);
vecY_pix_interp = linspace(vecY_pix(1),vecY_pix(end),9);

% %get colormap(s)
% cellColorMaps = RH_ColorMaps;
%%
%loop through channels
set(0,'DefaultFigureWindowStyle','docked')
for intClu = 1:intNumClu
    matAvgRespAll_interp = matAvgRespAll(:,:,intClu);
    if max(matAvgRespAll_interp,[],'all') > 10
    figure; hold on;
    title(['Cluster: ' num2str(intClu)]);
    imagesc(vecX_pix_interp,vecY_pix_interp,matAvgRespAll_interp);
    set(gca, 'YDir','reverse');
%     colormap(cellColorMaps{2});
    cb=colorbar;
    cb.Label.String='spks/s';
    axis image
    fixfig;
    else
        continue
    end
end

%%
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
