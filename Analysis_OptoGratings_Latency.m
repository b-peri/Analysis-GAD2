[FileNames, PathName] = uigetfile('*.mat', 'MultiSelect', 'on');

NumFiles = numel(FileNames);

DiffTimes = [];

for idx = 1:NumFiles % For each recording
    load(fullfile(PathName, FileNames{idx}));
    % Load in Rec Data
    structEP = sAP.cellBlock{1,1};  % structEP contains data on a single recording block (incl. stim onset/offset times)
    vecStimOnSecs = structEP.vecStimOnTime; % For OptoGratings, synchronized stimulus onset times
    vecStimOffSecs = structEP.vecStimOffTime; % For Optogratings, synchronized stimulus OFFset times
    vecLaserOn = structEP.vecOptoOn; % Logical array: Tells whether opto was on for each trial!
    vecLaserOnSecs = structEP.vecLaserOnTime; % Laser on times
    vecLaserOffSecs = structEP.vecLaserOffTime; % Laser off times
    
    % Compute Differences & Append to Main Vec
    DiffTimes = [DiffTimes; (vecLaserOnSecs - vecStimOnSecs(vecLaserOn))'];
end

DataOut.AllMice.Overall.LatencyOpto = mean(DiffTimes);
DataOut.AllMice.Overall.LatencyOptoSEM = std(DiffTimes)/sqrt(length(DiffTimes));