function plotRasterSplit(vecSpikes,vecTrialStarts,dblTrialDur,vecTrialType,matColor,intPlotMaxSpikes)
%plotRaster Makes raster plot
%syntax: plotRaster(vecSpikes,vecTrialStarts,dblTrialDur,intPlotMaxSpikes)
%	input:
%	- vecSpikes; spike times (s)
%	- vecTrialStarts: trial start times (s)
%	- dblTrialDur: max trial dur
%   - vecTrialType
%   - matColor
%	- intPlotMaxSpikes: maximum number of spikes to plot (default inf)
%
%Version history:
%1.0 - 18 June 2019
%	Created by Jorrit Montijn
%1.1 - 2 Dec 2019
%	Added max spike number to plot [by JM]
%2.1 - 27 July 2023
%   Added the possibility to split trials in different types [by RH]

%get inputs
if ~exist('intPlotMaxSpikes','var') || isempty(intPlotMaxSpikes)
    intPlotMaxSpikes = inf;
end
if ~exist('dblTrialDur','var') || isempty(dblTrialDur)
    dblTrialDur = median(diff(vecTrialStarts));
end
if ~exist('vecTrialType','var') || isempty(vecTrialType)
    vecTrialType = ones(size(vecTrialStarts));
end
% if ~isempty(matColor)% && (size(matColor,1) ~= 3 || size(matColor,2) ~= 3)
%     error('Check matColor!')
% end

[vecIdx,vecUnique,~,~,~] = val2idx(vecTrialType);
intType = numel(vecUnique);
if ~exist('matColor','var') || isempty(matColor) || size(matColor,1) ~= intType
    matColor = [0 0 0; lines(intType-1)];
end

%subselect
if numel(vecSpikes) > intPlotMaxSpikes
    vecSpikes = vecSpikes(sort(randperm(numel(vecSpikes),intPlotMaxSpikes)));
end

cla;
hold all;
intOffset = 0;
for intTrialType = 1:numel(vecUnique)

    %get trial onsets
    vecThisTrialStarts = vecTrialStarts(vecIdx==intTrialType);

    %get spike times in subset of trials
    [vecThisTrialPerSpike,vecThisTimePerSpike] = getSpikesInTrial(vecSpikes,vecThisTrialStarts,dblTrialDur);

    %plot per trial
    for intTrial=1:numel(vecThisTrialStarts)
        vecTimes = vecThisTimePerSpike(vecThisTrialPerSpike==intTrial);
        vecTimes(vecTimes>dblTrialDur)=[];
        line([vecTimes(:)';vecTimes(:)'],[(intTrial+intOffset)*ones(1,numel(vecTimes))-0.5;(intTrial+intOffset)*ones(1,numel(vecTimes))+0.5],...
            'Color',matColor(intTrialType,:),'LineWidth',1.5);
    end
    intOffset = intOffset + numel(vecThisTrialStarts);
end
hold off

%set fig props
ylim([0.5 numel(vecTrialStarts)+0.5]);
xlim([0 dblTrialDur]);
xlabel('Time after trial start (s)');
ylabel('Trial #');
fixfig(gca);
end

