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


%% Match Visually Responsive Units (Full-Field Gratings) to Optotagged Units 

Subject = [];
RecDate = [];
ClusterN = [];
CellClass_OT = [];
CellClass_OG = [];
Area = [];
zeta_p_OG = [];
SpontRate = [];
ER_Off = [];
ER_On = [];
PSTHMean_Off = [];
PSTHSEM_Off = [];
PSTHMean_On = [];
PSTHSEM_On = [];
PSTHBinCenters = [];

for i = 1:height(DataOut_OG_GAD2.ClusterData)
    row = DataOut_OG_GAD2.ClusterData(i,:);
    TagCell = DataOut_OT_50ms.ClusterData(DataOut_OT_50ms.ClusterData.Subject == row.Subject ...
        & DataOut_OT_50ms.ClusterData.RecDate == row.RecDate & ...
        DataOut_OT_50ms.ClusterData.ClusterN == row.ClusterN,:);
    if ~isempty(TagCell)
        CellClass_OT = [CellClass_OT; TagCell.CellClass];
    else
        CellClass_OT = [CellClass_OT; "Untagged"];
    end
    
    Subject = [Subject; row.Subject];
    RecDate = [RecDate; row.RecDate];
    ClusterN = [ClusterN; row.ClusterN];
    CellClass_OG = [CellClass_OG; row.EffectOpto];
    Area = [Area; row.Area];
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

DataOut_Matched.All = table(Subject,RecDate,ClusterN,CellClass_OT,CellClass_OG,Area, zeta_p_OG,SpontRate,PSTHMean_Off, ...
    PSTHSEM_Off, PSTHMean_On, PSTHSEM_On, PSTHBinCenters);

DataOut_Matched.Untagged = DataOut_Matched.All(DataOut_Matched.All.CellClass_OT == "Untagged",:);
DataOut_Matched.GAD2 = DataOut_Matched.All(DataOut_Matched.All.CellClass_OT == "GAD2+",:);
DataOut_Matched.Act = DataOut_Matched.All(DataOut_Matched.All.CellClass_OT == "Activated",:);
DataOut_Matched.Inh = DataOut_Matched.All(DataOut_Matched.All.CellClass_OT == "Inhibited",:);
DataOut_Matched.Oth = DataOut_Matched.All(DataOut_Matched.All.CellClass_OT == "Other",:);

%% --- Z-Scored PSTH Values ---
% Untagged
Un_PSTH_off_z = (test_Un.PSTHMean_Off - mean(test_Un.SpontRate))./std(test_Un.SpontRate);
Un_PSTH.PSTHMean_Off_z = mean(Un_PSTH_off_z, 1);
Un_PSTH.PSTHSEM_Off_z = std(Un_PSTH_off_z)./sqrt(height(Un_PSTH_off_z));

Un_PSTH_on_z = (test_Un.PSTHMean_On - mean(test_Un.SpontRate))./std(test_Un.SpontRate);
Un_PSTH.PSTHMean_On_z = mean(Un_PSTH_on_z, 1);
Un_PSTH.PSTHSEM_On_z = std(Un_PSTH_on_z)./sqrt(height(Un_PSTH_on_z));

Un_PSTH.PSTHBinCenters = test_Un.PSTHBinCenters(1,:);

% Activated
Act_PSTH_off_z = (test_Act.PSTHMean_Off - mean(test_Act.SpontRate))./std(test_Act.SpontRate);
Act_PSTH.PSTHMean_Off_z = mean(Act_PSTH_off_z, 1);
Act_PSTH.PSTHSEM_Off_z = std(Act_PSTH_off_z)./sqrt(height(Act_PSTH_off_z));

Act_PSTH_on_z = (test_Act.PSTHMean_On - mean(test_Act.SpontRate))./std(test_Act.SpontRate);
Act_PSTH.PSTHMean_On_z = mean(Act_PSTH_on_z, 1);
Act_PSTH.PSTHSEM_On_z = std(Act_PSTH_on_z)./sqrt(height(Act_PSTH_on_z));

Act_PSTH.PSTHBinCenters = test_Act.PSTHBinCenters(1,:);

% Inhibited
Inh_PSTH_off_z = (test_Inh.PSTHMean_Off - mean(test_Inh.SpontRate))./std(test_Inh.SpontRate);
Inh_PSTH.PSTHMean_Off_z = mean(Inh_PSTH_off_z, 1);
Inh_PSTH.PSTHSEM_Off_z = std(Inh_PSTH_off_z)./sqrt(height(Inh_PSTH_off_z));

Inh_PSTH_on_z = (test_Inh.PSTHMean_On - mean(test_Inh.SpontRate))./std(test_Inh.SpontRate);
Inh_PSTH.PSTHMean_On_z = mean(Inh_PSTH_on_z, 1);
Inh_PSTH.PSTHSEM_On_z = std(Inh_PSTH_on_z)./sqrt(height(Inh_PSTH_on_z));

Inh_PSTH.PSTHBinCenters = test_Inh.PSTHBinCenters(1,:);

% Other
Oth_PSTH_off_z = (test_Oth.PSTHMean_Off - mean(test_Oth.SpontRate))./std(test_Oth.SpontRate);
Oth_PSTH.PSTHMean_Off_z = mean(Oth_PSTH_off_z, 1);
Oth_PSTH.PSTHSEM_Off_z = std(Oth_PSTH_off_z)./sqrt(height(Oth_PSTH_off_z));

Oth_PSTH_on_z = (test_Oth.PSTHMean_On - mean(test_Oth.SpontRate))./std(test_Oth.SpontRate);
Oth_PSTH.PSTHMean_On_z = mean(Oth_PSTH_on_z, 1);
Oth_PSTH.PSTHSEM_On_z = std(Oth_PSTH_on_z)./sqrt(height(Oth_PSTH_on_z));

Oth_PSTH.PSTHBinCenters = test_Oth.PSTHBinCenters(1,:);

% --- Normalized PSTH Values (Baseline-Subtracted & ...) ---
nClust = height(test);
norm = max(test.PSTHMean_Off,[],2);
PSTHMean_Off_norm = (test.PSTHMean_Off - test.SpontRate)./norm;

