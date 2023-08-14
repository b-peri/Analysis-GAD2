function plotRaster(vecSpikes,vecTrialStarts,dblTrialDur,intPlotMaxSpikes, Color)
	%plotRaster Makes raster plot
	%syntax: plotRaster(vecSpikes,vecTrialStarts,dblTrialDur,intPlotMaxSpikes)
	%	input:
	%	- vecSpikes; spike times (s)
	%	- vecTrialStarts: trial start times (s)
	%	- dblTrialDur: max trial dur
	%	- intPlotMaxSpikes: maximum number of spikes to plot (default inf)
	%
	%Version history:
	%1.0 - 18 June 2019
	%	Created by Jorrit Montijn
	%1.1 - 2 Dec 2019
	%	Added max spike number to plot [by JM]
	
	%get inputs
	if ~exist('intPlotMaxSpikes','var') || isempty(intPlotMaxSpikes)
		intPlotMaxSpikes = inf;
	end
	if ~exist('dblTrialDur','var') || isempty(dblTrialDur)
		dblTrialDur = median(diff(vecTrialStarts));
    end
    if ~exist('Color','var') || isempty(Color)
        Color = 'k'
    end
	
	%subselect
	if numel(vecSpikes) > intPlotMaxSpikes
		vecSpikes = vecSpikes(sort(randperm(numel(vecSpikes),intPlotMaxSpikes)));
	end
	
	%get spike times in trials
	[vecTrialPerSpike,vecTimePerSpike] = getSpikesInTrial(vecSpikes,vecTrialStarts,dblTrialDur);
	
	%plot per trial
	cla;
	hold all;
	for intTrial=1:numel(vecTrialStarts)
		vecTimes = vecTimePerSpike(vecTrialPerSpike==intTrial);
		vecTimes(vecTimes>dblTrialDur)=[];
        vecTimes(vecTimes > 0.049 & vecTimes < 0.051) = [];
		line([vecTimes(:)';vecTimes(:)'],[intTrial*ones(1,numel(vecTimes))-0.5;intTrial*ones(1,numel(vecTimes))+0.5],'Color',Color,'LineWidth',1.5);
	end
	hold off
	
	%set fig props
	ylim([0.5 numel(vecTrialStarts)+0.5]);
	xlim([0 dblTrialDur]);
	xlabel('Time after trial start (s)');
	ylabel('Trial #');
	fixfig(gca);
end

