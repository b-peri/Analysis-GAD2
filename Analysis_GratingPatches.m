%% Visual Receptive Field (Grating Patches) Analysis

% [] Think about selection criteria for neurons to be analyzed... Can start
% by just using same channels that passed Zeta threshold in OptoGratings?

%% get stimulus onset times

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

% Stim Center Locations (?) -> Come back to confirm this
vecX_pix = unique(vecUniqueRects(:,1))+(vecUniqueRects(1,3)-unique(vecUniqueRects(1,1)))/2;
vecY_pix = unique(vecUniqueRects(:,2))+(vecUniqueRects(1,4)-unique(vecUniqueRects(1,2)))/2;

%% loop through data

%Check cells are visually responsive!

matAvgRespAll = NaN(numel(vecY_pix),numel(vecX_pix),intNumClu);
for intCh = 1:intNumClu
    vecSpikesCh = sAP.sCluster(intCh).SpikeTimes; % Load in spike times
    vecRate = zeros(1,structEP.intTrialNum); % Initialize horizontal vector, each element contains total # of spikes during stimulus presentation on that trial
    for intTrial = 1:structEP.intTrialNum
        vecSpikeT = vecSpikesCh(vecSpikesCh>vecStimOnSecs(intTrial)&vecSpikesCh<vecStimOffSecs(intTrial)); % Grabs spike times that lie within stim presentation time for this trial
        vecRate(intTrial) = numel(vecSpikeT)/(vecStimOffSecs(intTrial)-vecStimOnSecs(intTrial)); % Counts number of spikes
    end
    matAvgResp = NaN(numel(vecY_pix),numel(vecX_pix)); % NaN Array where each element corresponds to a unique stim center; Will contain average number of spikes per stim location

    for intLoc = vecUniqueStims
        matAvgResp(intLoc) = mean(vecRate(vecStimIdx==intLoc));
    end
    sParams.dblSecsFromPrevStimOff = 0.1; %s, for computing unit's baseline rate
    dblRateSpontaneous = computeRateSpontaneous(vecSpikeT,vecStimOnSecs,vecStimOffSecs,sParams);
    matAvgRespAll(:,:,intCh) = matAvgResp-dblRateSpontaneous; % Mean response minus baseline firing rate for that neuron/cluster!
end

%% plot data

% COMMENTED OUT FOR NOW -> We should in theory have sufficient spatial
% resolution atm to not need interpolation

%interpolate
% vecX_pix_interp = linspace(vecX_pix(1),vecX_pix(end),16);
% vecY_pix_interp = linspace(vecY_pix(1),vecY_pix(end),9);

%get colormap(s)
% cellColorMaps = RH_ColorMaps;
%%
%loop through channels
set(0,'DefaultFigureWindowStyle','docked')
for intCh = 1:intNumClu
    matAvgRespAll_interp = matAvgRespAll(:,:,intCh);
    figure; hold on;
    title(['Channel: ' num2str(intCh)]);
    imagesc(vecX_pix,vecY_pix,matAvgRespAll_interp);
    set(gca, 'YDir','reverse');
    % colormap(cellColorMaps{2});
    cb=colorbar;
    cb.Label.String='spks/s';
    axis image
    fixfig;
end

%%