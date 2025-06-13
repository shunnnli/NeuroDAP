% Shun_figure_patch.m

%% Define data path

clear; close all;
addpath(genpath(osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Analysis/NeuroDAP/Methods')));

rootPath = osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Project valence/Patch/Combined');
[~,~,~,~,blueWhiteRed,~,bluePurpleRed] = loadColors;
today = char(datetime('today','Format','yyyyMMdd'));

% Select cell and DA trend data
expPath = uipickfiles('FilterSpec',rootPath,'Prompt','Select combined mat files');
DAtrend_path = uipickfiles('FilterSpec',osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Project valence/Results/Summary-LastSessions'),...
                           'Prompt','Select DAtrend.mat')';

% Load combined_cells and combined_epochs
for i = 1:length(length(expPath))
    dirsplit = split(expPath{i},filesep); 
    filename = dirsplit{end};
    disp(['Loading: ', filename]);
    load(fullfile(expPath{i}));
    disp(['Finished: ',filename, ' loaded']);
end
resultPath = fileparts(expPath{1});
disp(strcat('resultPath: ',resultPath));

% Load DAtrend struct
load(DAtrend_path{1});
disp('Finished: DAtrend loaded');

%% (Optional) Extract response trace

close all;
animalRange = 'SL206';
[~,~] = getResponseTraces(combined_cells,animalRange=animalRange,plot=true);

%% (Analysis) Calculate animal avg EI charge & amplitude index

% EPSC_peaks = cellfun(@(x) x.summary.EPSC.peakAvg, combined_cells.Stats);
% IPSC_peaks = cellfun(@(x) x.summary.IPSC.peakAvg, combined_cells.Stats);
% EPSC_aucs = cellfun(@(x) x.summary.EPSC.aucAvg, combined_cells.Stats);
% IPSC_aucs = cellfun(@(x) x.summary.IPSC.aucAvg, combined_cells.Stats);

EPSC_stats = getCellStats(combined_cells,'exci',type='min',average=true);
IPSC_stats = getCellStats(combined_cells,'inhi',type='max',average=true);
EPSC_peaks = vertcat(EPSC_stats.response);
IPSC_peaks = vertcat(IPSC_stats.response);

EPSC_stats = getCellStats(combined_cells,'exci',type='auc',average=true);
IPSC_stats = getCellStats(combined_cells,'inhi',type='auc',average=true);
EPSC_aucs = vertcat(EPSC_stats.response);
IPSC_aucs = vertcat(IPSC_stats.response);

% EI statistics
% 1. EPSC + IPSC
EIsum_peaks = IPSC_peaks + EPSC_peaks;
EIsum_aucs = IPSC_aucs + EPSC_aucs;

% 2. EPSC/IPSC
% EIratio_peaks = abs(EPSC_peaks) ./ abs(IPSC_peaks);
% EIratio_aucs = abs(EPSC_aucs) ./ abs(IPSC_aucs);

% 3. EI index (-1: all excitatory & 1: all inhibitory)
EIindex_peaks = (abs(IPSC_peaks)-abs(EPSC_peaks)) ./ (abs(IPSC_peaks)+abs(EPSC_peaks));
EIindex_aucs = (abs(IPSC_aucs)-abs(EPSC_aucs)) ./ (abs(IPSC_aucs)+abs(EPSC_aucs));

% 4. z-scored EI index
bootstrapping = @(B,nBoot) reshape(B(randi(numel(B), numel(B), nBoot)), [], 1);
centering = @(X,B) (X - mean(B)) ./ std(B);
% Step 1: calculate baseline distribution
randomIdx = strcmpi('Random',combined_cells.Task) & combined_cells.Analyze;
baselineEPSC_peaks = bootstrapping(EPSC_peaks(randomIdx), 5000); %EPSC_peaks(randomIdx);
baselineIPSC_peaks = bootstrapping(IPSC_peaks(randomIdx), 5000);%IPSC_peaks(randomIdx);
baselineEPSC_aucs = bootstrapping(EPSC_aucs(randomIdx), 5000);%EPSC_aucs(randomIdx);
baselineIPSC_aucs = bootstrapping(IPSC_aucs(randomIdx), 5000);%IPSC_aucs(randomIdx);
% Step 2: calculate zscore of EPSC and IPSC stats based on baseline
EPSC_peaks_centered = centering(EPSC_peaks,baselineEPSC_peaks);
IPSC_peaks_centered = centering(IPSC_peaks,baselineIPSC_peaks);
EPSC_aucs_centered = centering(EPSC_aucs,baselineEPSC_aucs);
IPSC_aucs_centered = centering(IPSC_aucs,baselineIPSC_aucs);
% Step 3: calculate EI index
EIindex_peaks_centered = (abs(IPSC_peaks_centered)-abs(EPSC_peaks_centered)) ./ (abs(IPSC_peaks_centered)+abs(EPSC_peaks_centered));
EIindex_aucs_centered = (abs(IPSC_aucs_centered)-abs(EPSC_aucs_centered)) ./ (abs(IPSC_aucs_centered)+abs(EPSC_aucs_centered));


% Calculate animal average
animalList = unique({DAtrend.animal})';
animalEIindex_peaks = cellfun(@(a) mean(EIindex_peaks(strcmp(combined_cells.Animal, a))), animalList);
animalEIindex_aucs = cellfun(@(a) mean(EIindex_aucs(strcmp(combined_cells.Animal, a))), animalList);

% Calculate cell count per animal
cellsWithDAIdx = ismember(combined_cells.Animal, animalList);
cellsWithDA = combined_cells(cellsWithDAIdx,:);
nCells = groupcounts(cellsWithDA.Animal);

%% Figure 3





%% Fig 3: select sample cells for each condition

randomSampleAnimal = {'SL255','SL285','SL333'}; randomSampleCell = [7,2,10];
rewardSampleAnimal = {'SL344','SL343','SL354'}; rewardSampleCell = [5,1,4];
punishSampleAnimal = {'SL207','SL207','SL213'}; punishSampleCell = [3,6,10];

randomCellIdx = cellfun(@(anm,cnum) ...
                find(strcmp(combined_cells.Animal, anm) & combined_cells.Cell == cnum), ...
                randomSampleAnimal, num2cell(randomSampleCell));
rewardCellIdx = cellfun(@(anm,cnum) ...
                find(strcmp(combined_cells.Animal, anm) & combined_cells.Cell == cnum), ...
                rewardSampleAnimal, num2cell(rewardSampleCell));
punishCellIdx = cellfun(@(anm,cnum) ...
                find(strcmp(combined_cells.Animal, anm) & combined_cells.Cell == cnum), ...
                punishSampleAnimal, num2cell(punishSampleCell));

randomColor = [.7 .7 .7];
rewardColor = [0 158 115]./255;
punishColor = [135 104 247]./255;

close all; initializeFig(0.5,1); t = tiledlayout('flow');
t.TileSpacing = 'none';
t.Padding = 'none';
[~,~] = getResponseTraces(combined_cells,cellList=randomCellIdx,color=randomColor,...
                          plot=true,ylim=[-1550,1000],plotIndividual=true,newFig=false);
[~,~] = getResponseTraces(combined_cells,cellList=rewardCellIdx,color=rewardColor,...
                          plot=true,ylim=[-1500,1550],plotIndividual=true,newFig=false);
[~,~] = getResponseTraces(combined_cells,cellList=punishCellIdx,color=punishColor,...
                          plot=true,ylim=[-1100,500],plotIndividual=true,newFig=false);

% saveFigures(gcf,'Fig3-SampleCells',resultPath,savePDF=true,savePNG=true);


%% (Figure 3) Plot reversal cells

rewardAnimals = {'SL206','SL343','SL344','SL345','SL353','SL354','SL355'};
punishAnimals = {'SL207','SL213','SL253','SL342','SL044'};

randomIdx = find(strcmpi('Random',combined_cells.Task) & combined_cells.Analyze);
rewardIdx = find(ismember(combined_cells.Animal, rewardAnimals) & ...
               contains(combined_cells.Task, 'pairing'));
punishIdx = find(ismember(combined_cells.Animal, punishAnimals) & ...
               contains(combined_cells.Task, 'pairing'));
rewardCtrlIdx = [find(strcmpi('Reward control',combined_cells.Task) & combined_cells.Analyze);...
                 find(strcmpi('Reward pairing',combined_cells.Task) & ~combined_cells.Analyze)];
punishCtrlIdx = [find(strcmpi('Punish control',combined_cells.Task) & combined_cells.Analyze);...
                 find(strcmpi('Punish pairing',combined_cells.Task) & ~combined_cells.Analyze)];
baselineIdx = [randomIdx; rewardCtrlIdx; punishCtrlIdx];


% Plot cell EI based on tasks

% Select data to plot
figureName = 'CellEI-Reversal';
groupIdx = {randomIdx,rewardIdx,punishIdx,rewardCtrlIdx,punishCtrlIdx}; 
plotGroup = [1,1,1,0,0];
groupNames = {'Random','Reward','Punish','Reward Ctrl','Punish Ctrl'};

% Set color
opacity = 0.5;
rewardColor = [0 158 115]./255; rewardCtrlColor = 1 - opacity*(1-rewardColor);
punishColor = [135 104 247]./255; punishCtrlColor = 1 - opacity*(1-punishColor);
groupColors = {[.7 .7 .7],rewardColor,punishColor,...
                rewardCtrlColor,punishCtrlColor};

% Plot figure
close all;
plotCellEI(combined_cells,groupIdx,nboot=10000,centerEI=false,...
           plotGroup=plotGroup,groupColors=groupColors,groupNames=groupNames,...
           save=false,figureName=figureName,resultPath=resultPath,print=true);

%% Figure 4





%% (Fig 4) Plot sample DA stat vs trials and animal EPSC IPSC trace

% In paper: SL066, SL068, SL316
rowIdx = [8 9 6];
dotColor = [243 186 0]/255;
earlyColor = [102 204 204]/255;
endingColor = [0 153 153]/255;

close all; initializeFig(0.4,0.7);
tiledlayout(length(rowIdx),3);

for i = 1:length(rowIdx)
    r = rowIdx(i);
    trace = DAtrend(r).amp.raw;
    trace_smoothed = DAtrend(r).amp.smoothed;
    
    nexttile(i*3-2,[1 2]);
    scatter(1:length(trace),trace,150,dotColor,'filled',MarkerFaceAlpha=0.8); hold on
    plot(trace_smoothed,color=[dotColor,0.3],LineWidth=4);
    yline(0,'--',LineWidth=2,color=[.6 .6 .6]);
    xlim([1,length(trace)]);
    xticks([1, 20, length(trace)-19, length(trace)]);
    xticklabels({'1','20','n-19','n'}); yticks([]);
    xlabel('trials'); ylabel('DA amp');
    
    % Add shading
    yl = ylim; % get current y-limits
    % shade first 20 points
    rectangle('Position',[1,yl(1),19,yl(2)-yl(1)], ...
        'FaceColor', earlyColor, 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    % shade last 20 points
    rectangle('Position',[(length(trace)-19),yl(1),20,yl(2)-yl(1)], ...
        'FaceColor', endingColor,  'FaceAlpha', 0.1, 'EdgeColor', 'none');
    % draw vertical boundary lines
    xline(1, 'LineWidth',3, 'Color', earlyColor);
    xline(20, 'LineWidth',3, 'Color', earlyColor);
    xline(length(trace)-19, 'LineWidth',3, 'Color', endingColor);
    xline(length(trace), 'LineWidth',3, 'Color', endingColor);
    % push the rectangles to the bottom so they don't cover your lines/dots
    uistack(findobj(gca,'Type','rectangle'),'bottom');
    
    % Fit and plot best fit line for the first 20 dots
    x1 = 2:19;%1:20;
    y1 = trace_smoothed(x1);
    p1 = polyfit(x1, y1, 1);
    yfit1 = polyval(p1, x1);
    plot(x1, yfit1, color=[earlyColor,0.8], LineWidth=8);
    
    % Fit and plot best fit line for the last 20 dots
    x2 = (length(trace)-18):length(trace)-1;%(length(trace)-19):length(trace);
    y2 = trace_smoothed(x2);
    p2 = polyfit(x2, y2, 1);
    yfit2 = polyval(p2, x2);
    plot(x2, yfit2, color=[endingColor,0.8], LineWidth=8);

    % Plot EPSC and IPSC of each cell
    nexttile(i*3,[1 1]);
    [~,~] = getResponseTraces(combined_cells,animalRange=animalList(r),...
                              plot=true,overlay=true,newFig=false,...
                              plotLegend=false,plotNormalized=false,plotIndividual=false);
end
% saveFigures(gcf,'Fig4-schematic',resultPath,saveFIG=true,savePDF=true,savePNG=true);


%% (Fig 4) Calculate slopes of last n trials

trialWindow = 15:40;

recentSlope_max = zeros(length(DAtrend),1);
recentSlope_min = zeros(length(DAtrend),1);
recentSlope_avg = zeros(length(DAtrend),1);
recentSlope_amp = zeros(length(DAtrend),1);

distantSlope_max = zeros(length(DAtrend),1);
distantSlope_min = zeros(length(DAtrend),1);
distantSlope_avg = zeros(length(DAtrend),1);
distantSlope_amp = zeros(length(DAtrend),1);

sessionValue_max = zeros(length(DAtrend),1);
sessionValue_min = zeros(length(DAtrend),1);
sessionValue_avg = zeros(length(DAtrend),1);
sessionValue_amp = zeros(length(DAtrend),1);

% Calculate average slopes in trialWindow
for a = 1:length(DAtrend)
    % Get last n trial slope
    startTrial = max(1,DAtrend(a).nTrials-trialWindow);
    recentSlope_max(a) = mean(DAtrend(a).max.slopeMap.smoothed.map(startTrial,end));
    recentSlope_min(a) = mean(DAtrend(a).min.slopeMap.smoothed.map(startTrial,end));
    recentSlope_avg(a) = mean(DAtrend(a).avg.slopeMap.smoothed.map(startTrial,end));
    recentSlope_amp(a) = mean(DAtrend(a).amp.slopeMap.smoothed.map(startTrial,end));

    % Get first n trial slope
    endingTrial = min(DAtrend(a).nTrials,trialWindow);
    distantSlope_max(a) = mean(DAtrend(a).max.slopeMap.smoothed.map(1,endingTrial));
    distantSlope_min(a) = mean(DAtrend(a).min.slopeMap.smoothed.map(1,endingTrial));
    distantSlope_avg(a) = mean(DAtrend(a).avg.slopeMap.smoothed.map(1,endingTrial));
    distantSlope_amp(a) = mean(DAtrend(a).amp.slopeMap.smoothed.map(1,endingTrial));

    sessionValue_max(a) = mean(DAtrend(a).max.raw,'omitnan');
    sessionValue_min(a) = mean(DAtrend(a).min.raw,'omitnan');
    sessionValue_avg(a) = mean(DAtrend(a).avg.raw,'omitnan');
    sessionValue_amp(a) = mean(DAtrend(a).amp.raw,'omitnan');
end

% Define the eight slope vectors in a cell array and label them
slopeVectors = {recentSlope_max, recentSlope_min, recentSlope_avg, recentSlope_amp};
slopeNames = {'CueMax Slope (Short)', 'CueMin Slope (Short)', 'CueAvg Slope (Short)', 'CueAmp Slope (Short)'};

% Set color
upColor = [0, 158, 115]/255; %[122, 201, 67]/255;
downColor = [135, 104, 247]/255; %[84, 137, 45]/255; 
stableColor = [198, 156, 109]/255;

% Set group names
% numGroups = 2;
% groupNames = {'Down DA','Up DA'};
% groupColors = {downColor,upColor};

numGroups = 3;
groupNames = {'Down DA','Stable DA','Up DA'};
groupColors = {downColor,stableColor,upColor};


% (Analysis) kmeans clustering to define up/down/stable animals
% Preallocate cell arrays to save grouping results for each slope vector
downAnimals_kmeans = cell(numel(slopeVectors),1);
stableAnimals_kmeans = cell(numel(slopeVectors),1);
upAnimals_kmeans = cell(numel(slopeVectors),1);

initializeFig(0.5,1); tiledlayout(2,2);

for i = 1:numel(slopeVectors)
    nexttile;
    data = slopeVectors{i};
    
    % Perform k-means clustering with replicates for stability
    [groupIdx, groupCenters] = kmeans(data, numGroups, Replicates=100);
    
    % Sort group centers so that:
    %   Group 1: obvious negative (lowest center)
    %   Group 2: ambiguous (middle center)
    %   Group 3: obvious positive (highest center)
    [~, sortOrder] = sort(groupCenters);
    mapping = zeros(numGroups,1);
    mapping(sortOrder) = 1:numGroups;
    groupLabels = mapping(groupIdx);
    
    % Save grouping results
    if numGroups == 2
        downAnimals_kmeans{i} = find(groupLabels == 1);
        upAnimals_kmeans{i} = find(groupLabels == 2);
    elseif numGroups == 3
        downAnimals_kmeans{i} = find(groupLabels == 1);
        stableAnimals_kmeans{i} = find(groupLabels == 2);
        upAnimals_kmeans{i} = find(groupLabels == 3);
    end
    
    hold on;
    for j = 1:numGroups
        idx = groupLabels == j;
        scatter(find(idx), data(idx), 200, groupColors{j}, 'filled', DisplayName=groupNames{j});
    end
    hold off;
    
    title(slopeNames{i});
    xlabel('Animal');
    ylabel('Slope');
    legend('Location','southeast');
end
sgtitle('K-means classification');
% saveFigures(gcf,'K-means',resultPath,saveFIG=true,savePDF=true,savePNG=true);

% (Analysis) Plot cell EI clustered by kmeans
% Set up getCellIndices function
getCellIndices = @(indices) cell2mat(arrayfun(@(idx) find(strcmpi(combined_cells.Animal, animalList{idx})), indices, 'UniformOutput', false));

for i = 4:4%numel(slopeNames)

    % Get the animal indices from quant
    animalIdxDown = downAnimals_kmeans{i}; 
    cellIdxDown   = getCellIndices(animalIdxDown);
    animalIdxUp = upAnimals_kmeans{i};
    cellIdxUp     = getCellIndices(animalIdxUp);
    if numGroups == 3
        animalIdxStable = stableAnimals_kmeans{i}; 
        cellIdxStable = getCellIndices(animalIdxStable);
    end
    
    if numGroups == 2
        groupIdx = {cellIdxDown,cellIdxUp};
        plotGroup = [1, 1];
    elseif numGroups == 3
        groupIdx = {cellIdxStable,cellIdxDown,cellIdxUp};
        plotGroup = [1, 1, 1];
        groupNames = {'Stable DA','Down DA','Up DA'};
        groupColors = {stableColor,downColor,upColor};
    end
    
    % Create a figure name based on the current criterion (e.g., "CellEI-CueMaxSlopeShort")
    % Removing spaces and parentheses for a clean filename.
    figureName = ['CellEI-' regexprep(slopeNames{i}, '[ ()]', ''),'-kmeans'];
    close all;
    plotCellEI(combined_cells,groupIdx,centerEI=false,...
           plotGroup=plotGroup,groupColors=groupColors,groupNames=groupNames,...
           save=false,figureName=figureName,resultPath=resultPath,print=true);
end
















%% (Fig 3 Supp) QC summary

EPSC_qc = getCellStats(combined_cells,'exci',type='qc',average=true);
IPSC_qc = getCellStats(combined_cells,'inhi',type='qc',average=true);

% Plot Rs, Cm, Rm for each cell
EPSC_Rs = cellfun(@(x) x.Rs, {EPSC_qc.QC})'; 
EPSC_Rs_header = cellfun(@(x) x.Rs_headerString, {EPSC_qc.QC})'; 
% EPSC_Rs(90) = EPSC_Rs_header(90); EPSC_Rs_header(60) = EPSC_Rs(60);
EPSC_Rs_final = min(EPSC_Rs,EPSC_Rs_header);
IPSC_Rs = cellfun(@(x) x.Rs, {IPSC_qc.QC})';
IPSC_Rs_header = cellfun(@(x) x.Rs_headerString, {IPSC_qc.QC})';
IPSC_Rs_final = min(IPSC_Rs,IPSC_Rs_header);
Rs = min(EPSC_Rs_final,IPSC_Rs_final); 

% Plot QC verus response
EPSCcolor = [255 157 33]/255; %[122, 201, 67]/255;
IPSCcolor = [71 144 253]/255; %[84, 137, 45]/255; 
ctrlColor = [.7 .7 .7];
dotSize = 200;

initializeFig(0.4,0.3); tiledlayout('flow');

nexttile;
lmE   = fitlm(Rs, abs(EPSC_peaks));
slopeE = lmE.Coefficients.Estimate(2);
pvalE  = lmE.Coefficients.pValue(2);
lmI   = fitlm(Rs, abs(IPSC_peaks));
slopeI = lmI.Coefficients.Estimate(2);
pvalI  = lmI.Coefficients.pValue(2);
scatter(Rs,abs(EPSC_peaks),dotSize,EPSCcolor,"filled",MarkerFaceAlpha=0.5); hold on;
scatter(Rs,abs(IPSC_peaks),dotSize,IPSCcolor,"filled",MarkerFaceAlpha=0.5); hold on;
xFit = linspace(min(Rs), max(Rs), 100)';
yFitE = slopeE*xFit + lmE.Coefficients.Estimate(1);
yFitI = slopeI*xFit + lmI.Coefficients.Estimate(1);
plot(xFit, yFitE, 'Color', EPSCcolor, 'LineWidth', 5);
plot(xFit, yFitI, 'Color', IPSCcolor, 'LineWidth', 5);
xlabel('Rs'), ylabel('Amplitude'); %xlim([0 35]);
title('Rs vs amplitude');

nexttile;
lmE   = fitlm(Rs, abs(EPSC_aucs));
slopeE = lmE.Coefficients.Estimate(2);
pvalE  = lmE.Coefficients.pValue(2);
lmI   = fitlm(Rs, abs(IPSC_aucs));
slopeI = lmI.Coefficients.Estimate(2);
pvalI  = lmI.Coefficients.pValue(2);
scatter(Rs,abs(EPSC_aucs),dotSize,EPSCcolor,"filled",MarkerFaceAlpha=0.5); hold on;
scatter(Rs,abs(IPSC_aucs),dotSize,IPSCcolor,"filled",MarkerFaceAlpha=0.5); hold on;
xFit = linspace(min(Rs), max(Rs), 100)';
yFitE = slopeE*xFit + lmE.Coefficients.Estimate(1);
yFitI = slopeI*xFit + lmI.Coefficients.Estimate(1);
plot(xFit, yFitE, 'Color', EPSCcolor, 'LineWidth', 5);
plot(xFit, yFitI, 'Color', IPSCcolor, 'LineWidth', 5);
xlabel('Rs'), ylabel('Charge'); %xlim([0 40]);
title('Rs vs charge');

nexttile;
peakColor = [157 106 110]/255;
aucColor = [224 186 125]/255;
lmE   = fitlm(Rs, abs(EIindex_peaks));
slopeE = lmE.Coefficients.Estimate(2);
pvalE  = lmE.Coefficients.pValue(2);
lmI   = fitlm(Rs, abs(EIindex_aucs));
slopeI = lmI.Coefficients.Estimate(2);
pvalI  = lmI.Coefficients.pValue(2);
scatter(Rs,abs(EIindex_peaks),dotSize,peakColor,"filled",MarkerFaceAlpha=0.5); hold on;
scatter(Rs,abs(EIindex_aucs),dotSize,aucColor,"filled",MarkerFaceAlpha=0.5); hold on;
xFit = linspace(min(Rs), max(Rs), 100)';
yFitE = slopeE*xFit + lmE.Coefficients.Estimate(1);
yFitI = slopeI*xFit + lmI.Coefficients.Estimate(1);
plot(xFit, yFitE, 'Color', peakColor, 'LineWidth', 5);
plot(xFit, yFitI, 'Color', aucColor, 'LineWidth', 5);
xlabel('Rs'), ylabel('EI index'); %xlim([0 35]);
title('Rs vs EI index');


%% (Fig 3 Supp) Day of patching summary

rewardColor = [0 158 115]./255;
punishColor = [135 104 247]./255;

% days = {'1st pairing, D1','1st pairing, D2','1st pairing, D3+',...
%         '2nd pairing, D1','2nd pairing, D2','2nd pairing, D3+',...
%         '3rd pairing, D1','3rd pairing, D2','3rd pairing, D3+'}; 
% y = [0 3; 1 0; 5 2;...
%      7 1; 0 0; 8 5;...
%      1 1; 0 0; 0 0];

days = {'Pairing D1','Pairing D2-3','Pairing D3+'}; 
y = [8 5; 7 5; 7 2];

b = bar(y, EdgeColor='none');
b(1).FaceColor = rewardColor;
b(2).FaceColor = punishColor;
xticklabels(days);
box off

%% (Figure 3 Supp) Plot all cells


randomIdx = find(strcmpi('Random',combined_cells.Task) & combined_cells.Analyze);
rewardIdx = find(contains(combined_cells.Task, 'reward pairing',IgnoreCase=true));
punishIdx = find(contains(combined_cells.Task, 'punish pairing',IgnoreCase=true));
rewardCtrlIdx = [find(strcmpi('Reward control',combined_cells.Task) & combined_cells.Analyze);...
                 find(strcmpi('Reward pairing',combined_cells.Task) & ~combined_cells.Analyze)];
punishCtrlIdx = [find(strcmpi('Punish control',combined_cells.Task) & combined_cells.Analyze);...
                 find(strcmpi('Punish pairing',combined_cells.Task) & ~combined_cells.Analyze)];
baselineIdx = [randomIdx; rewardCtrlIdx; punishCtrlIdx];


% Plot cell EI based on tasks

% Select data to plot
figureName = 'CellEI-AllCells';
groupIdx = {randomIdx,rewardIdx,punishIdx,rewardCtrlIdx,punishCtrlIdx}; 
plotGroup = [1,1,1,0,0];
groupNames = {'Random','Reward','Punish','Reward Ctrl','Punish Ctrl'};

% Set color
opacity = 0.5;
rewardColor = [0 158 115]./255; rewardCtrlColor = 1 - opacity*(1-rewardColor);
punishColor = [135 104 247]./255; punishCtrlColor = 1 - opacity*(1-punishColor);
groupColors = {[.7 .7 .7],rewardColor,punishColor,...
                rewardCtrlColor,punishCtrlColor};

% Plot figure
close all;
plotCellEI(combined_cells,groupIdx,nboot=10000,centerEI=true,...
           plotGroup=plotGroup,groupColors=groupColors,groupNames=groupNames,...
           save=false,figureName=figureName,resultPath=resultPath,print=true);


%% (Figure 3 Supp) Plot reward or punish control cells


randomIdx = find(strcmpi('Random',combined_cells.Task));
rewardIdx = find(contains(combined_cells.Task, 'reward pairing',IgnoreCase=true));
punishIdx = find(contains(combined_cells.Task, 'punish pairing',IgnoreCase=true));
rewardCtrlIdx = [find(strcmpi('Reward control',combined_cells.Task))];
punishCtrlIdx = [find(strcmpi('Punish control',combined_cells.Task))];
baselineIdx = [randomIdx; rewardCtrlIdx; punishCtrlIdx];


% Plot cell EI based on tasks

% Select data to plot
figureName = 'CellEI-CtrlCells';
groupIdx = {randomIdx,rewardIdx,punishIdx,rewardCtrlIdx,punishCtrlIdx}; 
plotGroup = [1,0,0,1,1];
groupNames = {'Random','Reward','Punish','Reward Ctrl','Punish Ctrl'};

% Set color
opacity = 0.5;
rewardColor = [0 158 115]./255; rewardCtrlColor = 1 - opacity*(1-rewardColor);
punishColor = [135 104 247]./255; punishCtrlColor = 1 - opacity*(1-punishColor);
groupColors = {[.7 .7 .7],rewardColor,punishColor,...
                rewardCtrlColor,punishCtrlColor};

% Plot figure
close all;
plotCellEI(combined_cells,groupIdx,nboot=10000,centerEI=false,...
           plotGroup=plotGroup,groupColors=groupColors,groupNames=groupNames,...
           save=false,figureName=figureName,resultPath=resultPath,print=true);


%% (Fig 4 supp) plot cell EI removing reversal cells

% Find reversal cells
rewardAnimals = {'SL206','SL343','SL344','SL345','SL354','SL355'};
punishAnimals = {'SL207','SL213','SL253','SL342','SL044'};

% Exclude reversal cells
rewardIdx = find(strcmpi('Reward pairing',combined_cells.Task) & combined_cells.Analyze);
punishIdx = find(strcmpi('Punish pairing',combined_cells.Task) & combined_cells.Analyze);
rewardReversalIdx = find(ismember(combined_cells.Animal, rewardAnimals) & ...
               contains(combined_cells.Task, 'pairing'));
punishReversalIdx = find(ismember(combined_cells.Animal, punishAnimals) & ...
               contains(combined_cells.Task, 'pairing'));
allReversalIdx = [rewardReversalIdx;punishReversalIdx];
rewardIdx_excluded = setdiff(rewardIdx,allReversalIdx);
punishIdx_excluded = setdiff(punishIdx,allReversalIdx);
cellIdxStable_excluded = setdiff(cellIdxStable,allReversalIdx);
cellIdxUp_excluded = setdiff(cellIdxUp,allReversalIdx);
cellIdxDown_excluded = setdiff(cellIdxDown,allReversalIdx);


% Other categories
randomIdx = find(strcmpi('Random',combined_cells.Task) & combined_cells.Analyze);
rewardCtrlIdx = [find(strcmpi('Reward control',combined_cells.Task) & combined_cells.Analyze);...
                 find(strcmpi('Reward pairing',combined_cells.Task) & ~combined_cells.Analyze)];
punishCtrlIdx = [find(strcmpi('Punish control',combined_cells.Task) & combined_cells.Analyze);...
                 find(strcmpi('Punish pairing',combined_cells.Task) & ~combined_cells.Analyze)];
baselineIdx = [randomIdx; rewardCtrlIdx; punishCtrlIdx];


% Set color
opacity = 0.5;
rewardColor = [0 158 115]./255; rewardCtrlColor = 1 - opacity*(1-rewardColor);
punishColor = [135 104 247]./255; punishCtrlColor = 1 - opacity*(1-punishColor);
groupColors = {[.7 .7 .7],rewardColor,punishColor,...
                rewardCtrlColor,punishCtrlColor};


% % Plot cell EI based on task
% groupIdx = {randomIdx,rewardIdx_excluded,punishIdx_excluded,rewardCtrlIdx,punishCtrlIdx}; 
% plotGroup = [0,1,1,0,0];
% groupNames = {'Random','Reward','Punish','Reward Ctrl','Punish Ctrl'};
% % Plot figure
% close all;
% figureName = 'CellEI-Task-ReversalExcluded';
% plotCellEI(combined_cells,groupIdx,nboot=10000,centerEI=false,...
%            plotGroup=plotGroup,groupColors=groupColors,groupNames=groupNames,...
%            save=true,figureName=figureName,resultPath=resultPath,print=true);


% Plot cell EI based on DA slope
groupIdx = {cellIdxStable_excluded,cellIdxDown_excluded,cellIdxUp_excluded};
plotGroup = [1, 1, 1];
groupNames = {'Stable DA','Down DA','Up DA'};
groupColors = {stableColor,downColor,upColor};
close all;
figureName = 'CellEI-DAslope-ReversalExcluded';
plotCellEI(combined_cells,groupIdx,centerEI=false,...
       plotGroup=plotGroup,groupColors=groupColors,groupNames=groupNames,...
       save=true,figureName=figureName,resultPath=resultPath,print=true);

%% Fig 4 supp: compare last 30 trial, first 30 trial, value

trialWindow = 20;
sweepWindow = 20;
endingSlope = zeros(length(DAtrend),1); 
startingSlope = zeros(length(DAtrend),1);
sessionValue = zeros(length(DAtrend),1);

for a = 1:length(DAtrend)
    % Get last n trial slope
    startTrial = max(1,DAtrend(a).nTrials-trialWindow);
    endingSlope(a) = DAtrend(a).amp.slopeMap.smoothed.map(startTrial,end);

    % Get first n trial slope
    endingTrial = min(DAtrend(a).nTrials,trialWindow);
    startingSlope(a) = DAtrend(a).amp.slopeMap.smoothed.map(1,endingTrial);

    % Get session average DA value
    sessionValue(a) = mean(DAtrend(a).amp.raw,'omitnan');
end


% Plot figure
close all;
initializeFig(0.8,1); 
master = tiledlayout('flow');
endingColor = [0, 153, 153] / 255; %[128, 179, 255]./255;
startingColor = [102, 204, 204] / 255; %[152, 201, 163]./255;
valueColor = [.7 .7 .7];
sweepColor = [184 140 1]/255;
LineWidth = 5;
nboot = 1000;

% Plot first trial slope vs animal EI
X = animalEIindex_peaks; Y = startingSlope; color = startingColor;
[model,p_value,fit_boot] = fitScatter(X,Y,type='weighted',weights=nCells); slope = model(1);
[r_spearman_starting, p_spearman_starting] = corr(X(:), Y(:), 'Type', 'Spearman');

nexttile(master,1); 
children = tiledlayout(master,4,1);
children.Layout.Tile = 1; children.Layout.TileSpan = [1 1];
children.TileSpacing = 'compact'; children.Padding = 'tight'; axis off;

nexttile(children,1,[1 1]);
histogram(fit_boot(:,1),200,EdgeColor=addOpacity(color,0.5),FaceColor=addOpacity(color,0.5)); hold on; 
xline(slope,'-',['p=',num2str(p_value)],Color=color,LineWidth=3);
box off; c = gca; c.YAxis(1).Visible = 'off';

nexttile(children,2,[3 1]);
scatter(X,Y,200,color,"filled"); hold on;
xFit = linspace(min(X), max(X), 100);
yFit = polyval(model, xFit);
plot(xFit, yFit, color=color, LineWidth=LineWidth);
xlabel('Animal EI index');
ylabel('Distant DA slope');



% Plot last trial slope vs animal EI
X = animalEIindex_peaks; Y = endingSlope; color = endingColor;
[model,p_value,fit_boot] = fitScatter(X,Y,type='weighted',weights=nCells); slope = model(1);
[r_spearman_ending, p_spearman_ending] = corr(X(:), Y(:), 'Type', 'Spearman');

nexttile(master,2); 
children = tiledlayout(master,4,1);
children.Layout.Tile = 2; children.Layout.TileSpan = [1 1];
children.TileSpacing = 'compact'; children.Padding = 'tight'; axis off;

nexttile(children,1,[1 1]);
histogram(fit_boot(:,1),200,EdgeColor=addOpacity(color,0.5),FaceColor=addOpacity(color,0.5)); hold on; 
xline(slope,'-',['p=',num2str(p_value)],Color=color,LineWidth=3);
box off; c = gca; c.YAxis(1).Visible = 'off';

nexttile(children,2,[3 1]);
scatter(X,Y,200,color,"filled"); hold on;
xFit = linspace(min(X), max(X), 100);
yFit = polyval(model, xFit);
plot(xFit, yFit, color=color, LineWidth=LineWidth);
xlabel('Animal EI index');
ylabel('Recent DA slope');



% Plot average trial value vs animal EI
X = animalEIindex_peaks; Y = sessionValue; color = valueColor;
[model,p_value,fit_boot] = fitScatter(X,Y,type='weighted',weights=nCells); slope = model(1);
[r_spearman_value, p_spearman_value] = corr(X(:), Y(:), 'Type', 'Spearman');

nexttile(master,3); 
children = tiledlayout(master,4,1);
children.Layout.Tile = 3; children.Layout.TileSpan = [1 1];
children.TileSpacing = 'compact'; children.Padding = 'tight'; axis off;

nexttile(children,1,[1 1]);
histogram(fit_boot(:,1),200,EdgeColor=addOpacity(color,0.5),FaceColor=addOpacity(color,0.5)); hold on; 
xline(slope,'-',['p=',num2str(p_value)],Color=color,LineWidth=3);
box off; c = gca; c.YAxis(1).Visible = 'off';

nexttile(children,2,[3 1]);
scatter(X,Y,200,color,"filled"); hold on;
xFit = linspace(min(X), max(X), 100);
yFit = polyval(model, xFit);
plot(xFit, yFit, color=color, LineWidth=LineWidth);
xlabel('Animal EI index');
ylabel('Average DA amplitude');


nexttile(master,4);
plotScatterBar(1,r_spearman_starting,style='bar',color=startingColor,LineWidth=LineWidth);
text(2,r_spearman_starting,['p = ',num2str(p_spearman_starting)]);
plotScatterBar(2,r_spearman_ending,style='bar',color=endingColor,LineWidth=LineWidth);
text(1,r_spearman_ending,['p = ',num2str(p_spearman_ending)]);
plotScatterBar(3,r_spearman_value,style='bar',color=valueColor,LineWidth=LineWidth);
text(3,r_spearman_value,['p = ',num2str(p_spearman_value)]);
xticks([1 2 3]);
xticklabels({'Distant DA slope','Recent DA slope','Average DA amplitude'});
ylabel('Spearman correlation coefficient');


% See trial window effects
trialWindowLength = 10:80;
ending_slopes = zeros(length(trialWindowLength),1);
ending_pvals = zeros(length(trialWindowLength),1);
ending_spearman = zeros(length(trialWindowLength),1);
ending_spearman_p = zeros(length(trialWindowLength),1);
starting_slopes = zeros(length(trialWindowLength),1);
starting_pvals = zeros(length(trialWindowLength),1);
starting_spearman = zeros(length(trialWindowLength),1);
starting_spearman_p = zeros(length(trialWindowLength),1);
ending_vals = zeros(length(trialWindowLength),1);
ending_vals_pvals = zeros(length(trialWindowLength),1);
ending_spearman_val = zeros(length(trialWindowLength),1);
ending_spearman_val_p = zeros(length(trialWindowLength),1);
starting_vals = zeros(length(trialWindowLength),1);
starting_vals_pvals = zeros(length(trialWindowLength),1);
starting_spearman_val = zeros(length(trialWindowLength),1);
starting_spearman_val_p = zeros(length(trialWindowLength),1);

for i = 1:length(trialWindowLength)
    % Get slopes
    trialWindow = trialWindowLength(i);
    endingSlope = zeros(length(DAtrend),1); 
    startingSlope = zeros(length(DAtrend),1);
    sessionValue = zeros(length(DAtrend),1);
    for a = 1:length(DAtrend)
        % Get last n trial slope
        startTrial = max(1,DAtrend(a).nTrials-trialWindow);
        endingSlope(a) = DAtrend(a).amp.slopeMap.raw.map(startTrial,end);
        % Get first n trial slope
        endingTrial = min(DAtrend(a).nTrials,trialWindow);
        startingSlope(a) = DAtrend(a).amp.slopeMap.raw.map(1,endingTrial);

        % Get last n trial slope
        startTrial = max(1,DAtrend(a).nTrials-trialWindow);
        endingVal(a) = mean(DAtrend(a).amp.raw(startTrial:end));
        % Get first n trial slope
        endingTrial = min(DAtrend(a).nTrials,trialWindow);
        startingVal(a) = mean(DAtrend(a).amp.raw(1:endingTrial));
    end

    % Calculate last n trial fitted slopes
    Y = endingSlope;
    [model,ending_pvals(i),~] = fitScatter(X,Y,weights=nCells);
    ending_slopes(i) = model(1);
    [ending_spearman(i),ending_spearman_p(i)] = corr(X(:), Y(:), 'Type', 'Spearman');

    % Calculate first n trial fitted slopes
    Y = startingSlope;
    [model,starting_pvals(i),~] = fitScatter(X,Y,weights=nCells);
    starting_slopes(i) = model(1);
    [starting_spearman(i),starting_spearman_p(i)] = corr(X(:), Y(:), 'Type', 'Spearman');

    % Calculate last n trial fitted slopes
    Y = endingVal;
    [model,ending_vals_pvals(i),~] = fitScatter(X,Y,weights=nCells);
    ending_vals(i) = model(1);
    [ending_spearman_val(i),ending_spearman_val_p(i)] = corr(X(:), Y(:), 'Type', 'Spearman');

    % Calculate first n trial fitted slopes
    Y = startingVal;
    [model,starting_vals_pvals(i),~] = fitScatter(X,Y,weights=nCells);
    starting_vals(i) = model(1);
    [starting_spearman_val(i),starting_spearman_val_p(i)] = corr(X(:), Y(:), 'Type', 'Spearman');
end


% Plot correlation vs window
nexttile(master,7);
plot(trialWindowLength,ending_spearman,LineWidth=LineWidth,color=endingColor); hold on;
plot(trialWindowLength,starting_spearman,LineWidth=LineWidth,color=startingColor); hold on;
plot(trialWindowLength,ending_spearman_val,LineWidth=LineWidth,color=[.5 .5 .5]); hold on;
plot(trialWindowLength,starting_spearman_val,LineWidth=LineWidth,color=valueColor); hold on;
xlabel('trial window length');
ylabel('Spearman correlation coefficient');
box off

nexttile(master,8);
plot(trialWindowLength,ending_spearman_p,LineWidth=LineWidth,color=endingColor); hold on;
plot(trialWindowLength,starting_spearman_p,LineWidth=LineWidth,color=startingColor); hold on;
plot(trialWindowLength,ending_spearman_val_p,LineWidth=LineWidth,color=[.5 .5 .5]); hold on;
plot(trialWindowLength,starting_spearman_val_p,LineWidth=LineWidth,color=valueColor); hold on;
yline(0.05,'--',color=[0.7 0.7 0.7],LineWidth=LineWidth);
xlabel('trial window length');
ylabel('p value');
box off

% Plot moving 20 trial window (aligned to start)
nWindows  = 80-1;                                % # of windows you want
nAnimals  = numel(DAtrend);
movSlopes = nan(nAnimals, nWindows-sweepWindow+1);
movAvgs   = nan(nAnimals, nWindows-sweepWindow+1);
% precompute your slope‐kernel once
xc = (1:sweepWindow)' - mean(1:sweepWindow);
b  = flipud(xc)/sum(xc.^2);
for a = 1:nAnimals
    raw = DAtrend(a).amp.smoothed(:);
    % truncate or pad to exactly nWindow samples
    data = raw(1 : min(end,nWindows));
    data(end+1:nWindows,1) = nan;
    
    % compute on the flipped data
    revSlopes = conv(data, b, 'valid'); 
    csum      = cumsum([0; data]);
    revMeans  = (csum(sweepWindow+1:end) - csum(1:end-sweepWindow))/sweepWindow;
    
    % flip the results back
    movSlopes(a,:) = -revSlopes;
    movAvgs(a,:)   = revMeans;
end

nWindows  = nWindows-sweepWindow+1;
sweeping_slopes = zeros(length(nWindows),1);
sweeping_pvals = zeros(length(nWindows),1);
sweeping_spearman = zeros(length(nWindows),1);
sweeping_spearman_p = zeros(length(nWindows),1);
sweeping_vals = zeros(length(nWindows),1);
sweeping_vals_p = zeros(length(nWindows),1);
sweeping_vals_spearman = zeros(length(nWindows),1);
sweeping_vals_spearman_p = zeros(length(nWindows),1);

for w = 1:nWindows
    % Calculate fitted slopes with animal EI
    Y = movSlopes(:,w);
    [model,sweeping_pvals(w),~] = fitScatter(X,Y,weights=nCells);
    sweeping_slopes(w) = model(1);
    [X_clean,Y_clean] = removeNaNs(X,Y);
    [sweeping_spearman(w),sweeping_spearman_p(w)] = corr(X_clean(:), Y_clean(:), 'Type', 'Spearman');
    % Calculate average val with animal EI
    Y = movAvgs(:,w);
    [model,sweeping_vals_p(w),~] = fitScatter(X,Y,weights=nCells);
    sweeping_vals(w) = model(1);
    [X_clean,Y_clean] = removeNaNs(X,Y);
    [sweeping_vals_spearman(w),sweeping_vals_spearman_p(w)] = corr(X_clean(:), Y_clean(:), 'Type', 'Spearman');
end


nexttile(master,11);
plot(1:nWindows,sweeping_spearman,LineWidth=LineWidth,color=sweepColor); hold on;
plot(1:nWindows,sweeping_vals_spearman,LineWidth=LineWidth,color=valueColor); hold on;
xlabel('trials');
ylabel('correlation coefficient');
box off

nexttile(master,12);
plot(1:nWindows,sweeping_spearman_p,LineWidth=LineWidth,color=sweepColor); hold on;
plot(1:nWindows,sweeping_vals_spearman_p,LineWidth=LineWidth,color=valueColor); hold on;
yline(0.05,'--',color=[0.7 0.7 0.7],LineWidth=LineWidth);
xlabel('trials');
ylabel('p value');
box off


% Plot moving 20 trial window (aligned to end)
nWindows  = 80-1;                                % # of windows you want
nAnimals  = numel(DAtrend);
movSlopes = nan(nAnimals, nWindows-sweepWindow+1);
movAvgs   = nan(nAnimals, nWindows-sweepWindow+1);
% precompute your slope‐kernel once
xc = (1:sweepWindow)' - mean(1:sweepWindow);
b  = flipud(xc)/sum(xc.^2);
for a = 1:nAnimals
    raw = DAtrend(a).amp.smoothed(:);
    L   = min(numel(raw), nWindows);
    data = [nan(nWindows-L,1);raw(1:L)];
    
    % flip the vector end-for-start
    rev = data(end:-1:1);
    
    % compute on the flipped data
    revSlopes = conv(rev, b, 'valid'); 
    csum      = cumsum([0; rev]);
    revMeans  = (csum(sweepWindow+1:end) - csum(1:end-sweepWindow))/sweepWindow;
    
    % flip the results back
    movSlopes(a,:) = -revSlopes;
    movAvgs(a,:)   = revMeans;
end

nWindows  = nWindows-sweepWindow+1;
sweeping_slopes = zeros(length(nWindows),1);
sweeping_pvals = zeros(length(nWindows),1);
sweeping_spearman = zeros(length(nWindows),1);
sweeping_spearman_p = zeros(length(nWindows),1);
sweeping_vals = zeros(length(nWindows),1);
sweeping_vals_p = zeros(length(nWindows),1);
sweeping_vals_spearman = zeros(length(nWindows),1);
sweeping_vals_spearman_p = zeros(length(nWindows),1);

for w = 1:nWindows
    % Calculate fitted slopes with animal EI
    Y = movSlopes(:,w);
    [model,sweeping_pvals(w),~] = fitScatter(X,Y,weights=nCells);
    sweeping_slopes(w) = model(1);
    [X_clean,Y_clean] = removeNaNs(X,Y);
    [sweeping_spearman(w),sweeping_spearman_p(w)] = corr(X_clean(:), Y_clean(:), 'Type', 'Spearman');
    % Calculate average val with animal EI
    Y = movAvgs(:,w);
    [model,sweeping_vals_p(w),~] = fitScatter(X,Y,weights=nCells);
    sweeping_vals(w) = model(1);
    [X_clean,Y_clean] = removeNaNs(X,Y);
    [sweeping_vals_spearman(w),sweeping_vals_spearman_p(w)] = corr(X_clean(:), Y_clean(:), 'Type', 'Spearman');
end


nexttile(master,15);
plot(1:nWindows,sweeping_spearman,LineWidth=LineWidth,color=sweepColor); hold on;
plot(1:nWindows,sweeping_vals_spearman,LineWidth=LineWidth,color=valueColor); hold on;
xlabel('trials to end');
ylabel('correlation coefficient');
box off

nexttile(master,16);
plot(1:nWindows,sweeping_spearman_p,LineWidth=LineWidth,color=sweepColor); hold on;
plot(1:nWindows,sweeping_vals_spearman_p,LineWidth=LineWidth,color=valueColor); hold on;
yline(0.05,'--',color=[0.7 0.7 0.7],LineWidth=LineWidth);
xlabel('trials to end');
ylabel('p value');
box off
% saveFigures(gcf,'Fig4-all',resultPath,saveFIG=true,savePDF=true,savePNG=true);


%% (Fig 1 supp) For TRN patching

% Plot example cell
close all;
initializeFig(0.5,1); tiledlayout('flow');
[~,~] = getResponseTraces(combined_cells,cellList=18,...
                          plot=true,overlay=false,newFig=false,ylim=[-20,20],...
                          plotLegend=false,plotNormalized=false,plotIndividual=false);

% Plot currents
% Step 1: get all currents
epscRespAmp    = []; epscRespCharge = [];
epscBaseAmp    = []; epscBaseCharge = [];
ipscRespAmp    = []; ipscRespCharge = [];
ipscBaseAmp    = []; ipscBaseCharge = [];

for i = 1:height(combined_cells)
    % pull out the Vhold vector for this row:
    vhold = combined_cells.Vhold{i};    % e.g. [-70; 0] or [0; -70], etc.
    S     = combined_cells.Stats{i};    % struct with .response and .baseline
    
    for j = 1:numel(vhold)
        if vhold(j) == -70
            % EPSC
            epscRespAmp   (end+1,1) = mean(S.response.min{j});
            epscRespCharge(end+1,1) = mean(S.response.auc{j});    % or .auc
            epscBaseAmp   (end+1,1) = mean(S.baseline.min{j});
            epscBaseCharge(end+1,1) = mean(S.baseline.auc{j});
            
        elseif vhold(j) == 0
            % IPSC
            ipscRespAmp   (end+1,1) = mean(S.response.max{j});
            ipscRespCharge(end+1,1) = mean(S.response.auc{j});
            ipscBaseAmp   (end+1,1) = mean(S.baseline.max{j});
            ipscBaseCharge(end+1,1) = mean(S.baseline.auc{j});
        end
    end
end

exciColor = [255 157 33]/255; %[122, 201, 67]/255;
inhiColor = [71 144 253]/255; %[84, 137, 45]/255; 
ctrlColor = [.7 .7 .7];

% nexttile;
% plotScatterBar(1,epscRespAmp,color=exciColor);
% plotScatterBar(2,epscBaseAmp,color=ctrlColor);
% xticks([1,2]); xticklabels({'EPSC amp','Baseline amp'});
% ylim([-50,50]); ylabel('Current (pA)');
% 
% nexttile;
% plotScatterBar(1,ipscRespAmp,color=inhiColor);
% plotScatterBar(2,ipscBaseAmp,color=ctrlColor);
% xticks([1,2]); xticklabels({'IPSC amp','Baseline amp'});
% ylim([-50,50]); ylabel('Current (pA)');

nexttile;
plotScatterBar(1,epscRespCharge,color=exciColor);
plotScatterBar(2,epscBaseCharge,color=ctrlColor);
xticks([1,2]); xticklabels({'EPSC charge','Baseline charge'});
ylim([-0.5,0.5]); ylabel('Charge (pC)');

nexttile;
plotScatterBar(1,ipscRespCharge,color=inhiColor);
plotScatterBar(2,ipscBaseCharge,color=ctrlColor);
xticks([1,2]); xticklabels({'IPSC charge','Baseline charge'});
ylim([-0.5,0.5]); ylabel('Charge (pC)');

nexttile;
histogram(epscBaseCharge,20,EdgeColor=ctrlColor,FaceColor=ctrlColor); hold on
histogram(epscRespCharge,20,EdgeColor=exciColor,FaceColor=exciColor); hold on
xlabel('Charge (pC)'); ylabel('Count');

nexttile;
histogram(ipscBaseCharge,20,EdgeColor=ctrlColor,FaceColor=ctrlColor); hold on
histogram(ipscRespCharge,20,EdgeColor=inhiColor,FaceColor=inhiColor); hold on
xlabel('Charge (pC)'); ylabel('Count');