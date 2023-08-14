%% Figure 4: Further Analysis Optotagged Units

sAP = load("79155_20230512_AP.mat"); sAP = sAP.sAP;
GAD2_cells = [353; 366];

%% Zeta Scores of GAD2+ Units to Visual Stimulus & Laser (OG)

cols.GAD2 = [0 0.4470 0.7410];
cols.xline = [143 143 143]/255;

ZetaP_OG = [];
ZetaP_LaserOG = [];
s_count_OG = [];
test = [];

stim = sAP.cellBlock{1,1}.vecStimOnTime;
laser_trials = sAP.cellBlock{1,1}.vecOptoOn;
laser_stim = sAP.cellBlock{1,1}.vecLaserOnTime;
maxDur = 0.15;

for a = 1:2
    i = GAD2_cells(a);
    st = sAP.sCluster(i).SpikeTimes;

    % Draw Plot
    figure; hold on;
    % plotRasterSplit(st,stim-0.2,maxDur,~laser_trials,[cols.on_GAD2; cols.off]);
    plotRaster2(st,laser_stim-0.05,maxDur, [], cols.GAD2);
    xticks([0 0.05 0.1 0.15]);
    xticklabels(["-0.05" "0" "0.05" "0.1"]);
    xlabel('Time from laser onset (s)');
    xline(0.05, '--', 'Color', cols.xline);
    title(num2str(i));
    fixfig; hold off;

    ZetaP_OG_Cl = zetatest(st,stim(~laser_trials),0.9); % Compute zetatest for Stimuli w/o Opto
    if ~isempty(ZetaP_OG_Cl)
        ZetaP_OG = [ZetaP_OG; ZetaP_OG_Cl];
    else
        ZetaP_OG = [ZetaP_OG; 0];
    end

    ZetaP_LaserOG_Cl = zetatest(st,laser_stim,0.03);
    ZetaP_LaserOG = [ZetaP_LaserOG, ZetaP_LaserOG_Cl];

    s_count_OG_Cl = sum(st > stim(1) & st < stim(end));
    s_count_OG = [s_count_OG; s_count_OG_Cl];
    
    test = [test; i];
end

%% Start Figure & Set Colors

figure;

cols.GAD2 = [0 0.4470 0.7410];
cols.Act = [0.8500 0.3250 0.0980];
cols.Inh = [0.9290 0.6940 0.1250];
cols.xline = [143 143 143]/255;
cols.xreg = [200 200 200]/255;
cols.xreg2 = [184 202 214]/255;
cols.err_bar = [26 35 126]/255;
cols.opto = [0.3010 0.7450 0.9330];

%% Receptive Field Mapping of GAD2+ Units

% Prep Output
RF_centers = [];
RF_diam = [];

% Prep Stim Data
structEP = sAP.cellBlock{1,2};  
vecStimOnSecs = structEP.vecStimOnTime;
vecStimOffSecs = structEP.vecStimOffTime;

% Get grid data
vecUniqueRects = unique(structEP.vecDstRect','rows'); %unique dst rects
vecUniqueStims = 1:length(vecUniqueRects);
vecStimIdx = zeros(size(structEP.vecDstRect,2),1);
for intStim = 1:length(vecUniqueRects)
    vecStimIdx(ismember(structEP.vecDstRect',vecUniqueRects(intStim,:),'rows')) = vecUniqueStims(intStim);
end

% Stim Center Locations
vecX_pix = unique(vecUniqueRects(:,1))+(vecUniqueRects(1,3)-unique(vecUniqueRects(1,1)))/2;
vecY_pix = unique(vecUniqueRects(:,2))+(vecUniqueRects(1,4)-unique(vecUniqueRects(1,2)))/2;

% Prep
matAvgRespAll = NaN(numel(vecY_pix),numel(vecX_pix), 2);

for a = 1:2
    i = GAD2_cells(a);
    st = sAP.sCluster(i).SpikeTimes; % Load in spike times
    dblZetaP = zetatest(st,vecStimOnSecs,0.9); % Compute zetatest for Stimuli w/o Opto -> Visually responsive neurons
    if dblZetaP < 0.01
        vecRate = zeros(1,structEP.intTrialNum); % Initialize horizontal vector, each element contains total # of spikes during stimulus presentation on that trial
        for intTrial = 1:structEP.intTrialNum
            vecSpikeT = st(st>vecStimOnSecs(intTrial)&st<vecStimOffSecs(intTrial)); % Grabs spike times that lie within stim presentation time for this trial
            vecRate(intTrial) = numel(vecSpikeT)/(vecStimOffSecs(intTrial)-vecStimOnSecs(intTrial)); % Counts number of spikes + divides by trial duration to get spks/s
        end
        matAvgResp = NaN(numel(vecY_pix),numel(vecX_pix)); % NaN Array where each element corresponds to a unique stim center; Will contain average number of spikes per stim location
    
        for intLoc = vecUniqueStims
            matAvgResp(intLoc) = mean(vecRate(vecStimIdx==intLoc));
        end
        sParams.dblSecsFromPrevStimOff = 0.1; %s, for computing unit's baseline rate
        dblRateSpontaneous = computeRateSpontaneous(vecSpikeT,vecStimOnSecs,vecStimOffSecs,sParams);
        matAvgRespAll(:,:,a) = matAvgResp-dblRateSpontaneous; % Mean response minus baseline firing rate for that neuron/cluster!
    end
    
    % Interpolate Into 1° Resolution
    Xq = repmat(linspace(1,17,102), [54 1]);
    Yq = repmat(linspace(1,9,54)', [1 102]);

    matAvgRespInterp = interp2(matAvgRespAll(:,:,a),Xq,Yq);

    % Get Coordinates
    X_coord = 51:-1:-51;
    X_coord(X_coord == 0) = [];
    Y_coord = [27:-1:-27]';
    Y_coord(Y_coord == 0) = [];

    % RF Center
    [center_y, center_x] = find(matAvgRespInterp == max(matAvgRespInterp,[],"all"));
    RF_centers = [RF_centers; [Y_coord(center_y), X_coord(center_x)]];

    % RF Size
    [rf_y, rf_x] = find(matAvgRespInterp >= max(matAvgRespInterp,[],"all")/2);
    diam_y = max(rf_y) - min(rf_y);
    diam_x = max(rf_x) - min(rf_x);
    RF_diam = [RF_diam; diam_x diam_y];
    % RF_size = [RF_size; sum(matAvgRespInterp >= max(matAvgRespInterp,[],"all")/2)];

    %Plot RF
    subplot(2,3,1 + 3*(a-1)); hold on;
    imagesc(matAvgRespInterp); colormap('bone');
    xlabel("Azimuth (Vis. Deg.)"); ylabel("Elevation (Vis. Deg.)");
    xticks(linspace(2,101,5)); xticklabels(["-50" "-25" "0" "25" "50"]);
    yticks(linspace(3,52,5)); yticklabels(["25" "12.5" "0" "-12.5" "-25"]);
    % outline_rf = boundary(rf_x, rf_y); plot(outline_rf, 'r--');
    visboundaries(matAvgRespInterp >= max(matAvgRespInterp,[],"all")/2, "Color", cols.Inh, "EnhanceVisibility", 0, "LineStyle","--");
    viscircles([find(X_coord == -15) find(Y_coord == 1)], 0.25, "Color" ,cols.GAD2);
    viscircles([find(X_coord == -15) find(Y_coord == 1)], 16, "Color" ,cols.GAD2, "LineStyle", ":", "EnhanceVisibility", 0);
    axis image;
    % colorbar;
    fixfig; hold off;
end

%% Size Tuning

structEP = sAP.cellBlock{1,4};  
vecStimOnTime = structEP.vecStimOnTime;
vecStimOffTime = structEP.vecStimOffTime;
vecTrialStimSize = [structEP.sStimObject(structEP.vecTrialStimTypes).StimulusSize_deg];

vecStimSizes = structEP.sStimParams.vecStimulusSize_deg';
resp_all = NaN(length(vecStimSizes),1,2);

BinEdge = [-0.5 0 0.2];
binDur = [BinEdge(2) - BinEdge(1), BinEdge(3) - BinEdge(2)];

colors = [cols.GAD2 cols.opto];

for a = 1:2
    i = GAD2_cells(a);
    st = sAP.sCluster(i).SpikeTimes; % Load in spike times
    ZetaP_ST_Cl = zetatest(st,vecStimOnTime,0.9); 
    if ZetaP_ST_Cl < 0.01
        for b = 1:length(vecStimSizes)
            trial_sel = vecStimOnTime(vecTrialStimSize == vecStimSizes(b));
            sCounts = zeros(numel(trial_sel),2);
            for intTrial=1:numel(trial_sel)
	            vecTheseEdges = BinEdge + trial_sel(intTrial); % Add stim onset time for this trial to rel. bin edges to get absolute bin edges
	            [vecCounts,edges] = histcounts(st,vecTheseEdges);
	            sCounts(intTrial,1) = vecCounts(1); % Counts Spontaneous Rate
                sCounts_Opto(intTrial,2) = vecCounts(2); % Count for Visual Response
            end
            
            sRates = [sCounts_Opto(:,1)/binDur(1) sCounts_Opto(:,2)/binDur(2)];
            SpontRate = mean(sRates(:,1));
            EvokedResp = mean(sRates(:,2)) - SpontRate;

            resp_all(b,1,a) = EvokedResp;
        end
        subplot(2,3,2 + 3*(a-1)); hold on;
        bar(resp_all(:,1,a), "FaceColor", colors([1:3] + 3*(a-1))); xticks(1:7); xticklabels(num2str(vecStimSizes) + "°");
        xlabel("Stimulus Size"); ylabel("Evoked Response (spks/s)"); fixfig;
    end
end



% vecStimPosX_deg: -14.7009
% vecStimPosY_deg: 0.0597

%% Violations1ms: Check Whether GAD2+ Cell and Activated Cells are the Same

ActCells_Rec = ActTab(ActTab.Subject == '79155' & ActTab.RecDate == '20230512',:);
cells = [GadTab.ClusterN; ActCells_Rec.ClusterN];

combos = nchoosek(cells,2);
% combos = combos(2:end,:);

sAP = load("79155_20230512_AP.mat"); sAP = sAP.sAP;

vecViolations1ms = zeros(length(combos),1);
for i = 1:length(combos)
    st = sortrows([sAP.sCluster(combos(i,1)).SpikeTimes; sAP.sCluster(combos(i,2)).SpikeTimes]);
    sOut = getClusterQuality(st,0);
    vecViolations1ms(i) = sOut.dblViolIdx1ms;
end

combos_viol = [combos vecViolations1ms];


%% Presence

Presence = [];

for a = 1:2
    i = GAD2_cells(a);
    st = sAP.sCluster(i).SpikeTimes;
    BinEdges = [1:ceil(max(sAP.sCluster(353).SpikeTimes))];
    BinCenters = BinEdges(1:end-1) + (diff(BinEdges)/2);
    s_counts = histcounts(st, BinEdges);
    Presence_Cl = sum(s_counts > 0)/(numel(BinEdges) - 1);
    Presence = [Presence; Presence_Cl];

    subplot(2,3,3 + 3*(a-1)); 
    % histogram(st, BinEdges);
    plot(BinCenters, s_counts, "Color", colors([1:3] + 3*(a-1)));
    % title(string(i) + " Presence: " + string(Presence_Cl));
    xregion(sAP.cellBlock{1,1}.vecStimOnTime(1),sAP.cellBlock{1,1}.vecStimOffTime(end), "FaceColor", cols.xreg2); % OptoGrating Stim Times
    text(310,50,"OG", "FontSize", 15);
    xregion(sAP.cellBlock{1,2}.vecStimOnTime(1),sAP.cellBlock{1,2}.vecStimOffTime(end), "FaceColor", cols.xreg); % Optotagging Stim Times
    text(1760,50,"RF", "FontSize", 15);
    xregion(sAP.cellBlock{1,3}.vecStimOnTime(1),sAP.cellBlock{1,3}.vecStimOffTime(end), "FaceColor", cols.xreg2); % RF mapper Stim Times
    text(2650,50,"OT", "FontSize", 15)
    xregion(sAP.cellBlock{1,4}.vecStimOnTime(1),sAP.cellBlock{1,4}.vecStimOffTime(end), "FaceColor", cols.xreg); % Size Tuning Stim Times
    text(3570,50,"ST", "FontSize", 15)
    xlim([0 BinEdges(end)])
    xlabel('Time from start rec. (s)'); ylabel("Firing Rate (spks/s)"); fixfig
end

%% Save

saveas(gcf, 'D:\NIN\Analysis-GAD2\Plots\Figure4.png');
savefig(gcf, 'D:\NIN\Analysis-GAD2\Plots\Figure4.fig');



