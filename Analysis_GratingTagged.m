% Analysis GratingTagged

% %% Prep Output
% 
% DataOut_OGt.ClusterData = cell2table(cell(0,15), 'VariableNames', ...
%     {'Subject', 'RecDate', 'ClusterN', 'CellClass', 'Area', 'zeta_p_20ms', ...
%     'Peak_Lat', 'IFR_Rate', 'IFR_Time', 'SpontRate', 'SpontRate_STD', ...
%     'PSTHMean_20ms', 'PSTHSEM_20ms', 'PSTHBinCenters', 'NonStationarity'});
% 
% %% Select DataFiles
% 
% RecID = unique([DataOut_OT.ClusterData.Subject DataOut_OT.ClusterData.RecDate],'rows');
% FileNames = [RecID(:,1) + "_" + RecID(:,2) + "_AP.mat"];
% 
% for idx = 1:length(FileNames)
%     load(fullfile("C:\Software and Code\Data GAD2 Analysis", FileNames(idx)))
% 
%     % Get stimulus info
%     TaggedCl = DataOut_OT.ClusterData.ClusterN(DataOut_OT.ClusterData.Subject == RecID(idx,1) & DataOut_OT.ClusterData.RecDate == RecID(idx,2));
% 
%     % intNumClu = length(sAP.sCluster); % Grabs number of clusters/putative cells
%     structEP = sAP.cellBlock{1,1};  % structEP contains data on a single recording block (incl. stim onset/offset times)
%     vecStimOnSecs = structEP.vecStimOnTime; % For OptoGratings, synchronized stimulus onset times
%     vecStimOffSecs = structEP.vecStimOffTime; % For Optogratings, synchronized stimulus OFFset times
%     vecLaserOn = structEP.vecOptoOn; % Logical array: Tells whether opto was on for each trial!
%     vecLaserOnTime = structEP.vecLaserOnTime; % Laser on times
%     vecLaserOffTime = structEP.vecLaserOffTime; % Laser off times
% 
% end

%% Begin loop


%% TEST

% test = cell2table(cell(0,15), 'VariableNames', ...
%     {'Subject', 'RecDate', 'ClusterN', 'CellClass', 'Area', 'zeta_p_OG', ...
%     'SpontRate', 'ER_Off', 'ER_On', 'PSTHMean_Off', 'PSTHSEM_Off', ...
%     'PSTHMean_On', 'PSTHSEM_On', 'PSTHBinCenters', 'NonStationarity'});

Subject = [];
RecDate = [];
ClusterN = [];
CellClass = [];
zeta_p_OG = [];
SpontRate = [];
ER_Off = [];
ER_On = [];
PSTHMean_Off = [];
PSTHSEM_Off = [];
PSTHMean_On = [];
PSTHSEM_On = [];
PSTHBinCenters = [];

for i = 1:height(DataOut_OG.AllMice.ClusterData)
    row = DataOut_OG.AllMice.ClusterData(i,:);
    TagCell = DataOut_OT.ClusterData(Tags.Subject == row.Subject & Tags.RecDate == row.RecDate & Tags.ClusterN == row.ClusterN,:);
    if ~isempty(TagCell)
        CellClass = [CellClass; TagCell.CellClass];
    else
        CellClass = [CellClass; "Untagged"];
    end
    
    Subject = [Subject; row.Subject];
    RecDate = [RecDate; row.RecDate];
    ClusterN = [ClusterN; row.ClusterN];
    zeta_p_OG = [zeta_p_OG; row.zeta_p];
    SpontRate = [SpontRate; row.SpontRate];
    ER_Off = [ER_Off; row.ER_OptoOff];
    ER_On = [ER_On; row.ER_OptoOn];
    PSTHMean_Off = [PSTHMean_Off; row.PSTHMean_Off];
    PSTHSEM_Off = [PSTHSEM_Off; row.PSTHSEM_Off];
    PSTHMean_On = [PSTHMean_On; row.PSTHMean_On];
    PSTHSEM_On = [PSTHSEM_On; row.PSTHSEM_On];
    PSTHBinCenters = [PSTHBinCenters; row.PSTHBinCenters_Off];
end

test = table(Subject,RecDate,ClusterN,CellClass,zeta_p_OG,PSTHMean_Off, ...
    PSTHSEM_Off, PSTHMean_On, PSTHSEM_On, PSTHBinCenters);
test_Un = test(test.CellClass == "Untagged",:);
test_Act = test(test.CellClass == "Activated",:);
test_Inh = test(test.CellClass == "Inhibited",:);
test_Oth = test(test.CellClass == "Other",:);

% Normalized PSTH Values Untagged
nClust_Un = height(test_Un);
norm_Un = max(test_Un.PSTHMean_Off,[],2);
test_Un_PSTH.PSTHMean_Off_norm = mean(test_Un.PSTHMean_Off./norm_Un, 1);
test_Un_PSTH.PSTHSEM_Off_norm = std(test_Un.PSTHMean_Off./norm_Un, 1)/sqrt(nClust_Un);
test_Un_PSTH.PSTHMean_On_norm = mean(test_Un.PSTHMean_On./norm_Un, 1);
test_Un_PSTH.PSTHSEM_On_norm = std(test_Un.PSTHMean_On./norm_Un, 1)/sqrt(nClust_Un);
test_Un_PSTH.PSTHBinCenters = test_Un.PSTHBinCenters(1,:);

% Normalized PSTH Values Activated
nClust_Act = height(test_Act);
norm_Act = max(test_Act.PSTHMean_Off,[],2);
test_Act_PSTH.PSTHMean_Off_norm = mean(test_Act.PSTHMean_Off./norm_Act, 1);
test_Act_PSTH.PSTHSEM_Off_norm = std(test_Act.PSTHMean_Off./norm_Act, 1)/sqrt(nClust_Act);
test_Act_PSTH.PSTHMean_On_norm = mean(test_Act.PSTHMean_On./norm_Act, 1);
test_Act_PSTH.PSTHSEM_On_norm = std(test_Act.PSTHMean_On./norm_Act, 1)/sqrt(nClust_Act);
test_Act_PSTH.PSTHBinCenters = test_Act.PSTHBinCenters(1,:);

% Normalized PSTH Values Inhibited
nClust_Inh = height(test_Inh);
norm_Inh = max(test_Inh.PSTHMean_Off,[],2);
test_Inh_PSTH.PSTHMean_Off_norm = mean(test_Inh.PSTHMean_Off./norm_Inh, 1);
test_Inh_PSTH.PSTHSEM_Off_norm = std(test_Inh.PSTHMean_Off./norm_Inh, 1)/sqrt(nClust_Inh);
test_Inh_PSTH.PSTHMean_On_norm = mean(test_Inh.PSTHMean_On./norm_Inh, 1);
test_Inh_PSTH.PSTHSEM_On_norm = std(test_Inh.PSTHMean_On./norm_Inh, 1)/sqrt(nClust_Inh);
test_Inh_PSTH.PSTHBinCenters = test_Inh.PSTHBinCenters(1,:);

% Normalized PSTH Values Other
nClust_Oth = height(test_Oth);
norm_Oth = max(test_Oth.PSTHMean_Off,[],2);
test_Oth_PSTH.PSTHMean_Off_norm = mean(test_Oth.PSTHMean_Off./norm_Oth, 1);
test_Oth_PSTH.PSTHSEM_Off_norm = std(test_Oth.PSTHMean_Off./norm_Oth, 1)/sqrt(nClust_Oth);
test_Oth_PSTH.PSTHMean_On_norm = mean(test_Oth.PSTHMean_On./norm_Oth, 1);
test_Oth_PSTH.PSTHSEM_On_norm = std(test_Oth.PSTHMean_On./norm_Oth, 1)/sqrt(nClust_Oth);
test_Oth_PSTH.PSTHBinCenters = test_Oth.PSTHBinCenters(1,:);