disp(['*-*-*-* Running: scriptPlotAnalysisCellCC *-*-*-*'])

mouseAnalysisTablePath = [experimentDirectory filesep 'mouseAnalysis' filesep 'MouseAnalysisTable.mat'];
mouseParametersTablePath = [experimentDirectory filesep 'mouseAnalysis' filesep 'MouseParametersTable.mat'];

mouseAnalysisTable = load(mouseAnalysisTablePath);
mouseAnalysisTable = mouseAnalysisTable.mouseAnalysisTable;
mouseParametersTable = load(mouseParametersTablePath);
runGroupsParametersTable = mouseParametersTable.runGroupsParametersTable;

postTextSuptitle = '';
analyzedCells = unique(cell2mat(mouseAnalysisTable.cellName));

groupByCell = 0;

for iCell = 3%1:numel(analyzedCells)

iAnalyzedCell = analyzedCells(iCell);

cellName = ['cell' num2str(iAnalyzedCell)];

initializeFig(0.9,0.9);

nTileRows = 1;
nTileColumns = 1;
nTiles = nTileRows*nTileColumns;
masterLayout = tiledlayout(nTileRows,nTileColumns);
masterLayout.TileSpacing = 'compact';
masterLayout.Padding = 'compact';

colorBlue = [85, 161, 254]./255;
colorRed = [255, 50, 58]./255;
colorGrey = [192, 192, 192]./255;

colorBlueOverlay = [85, 161, 254]./255 * 0.6;
colorRedOverlay = [255, 50, 58]./255 * 0.6;
colorGreyOverlay = [192, 192, 192]./255 * 0.6;
    
%% ** Plot blue traces **

%% Collect data

cat1_queryRunFeature = 'areaVtPulse';
cat2_queryRunFeature = 'areaVtPulse';
cat3_queryRunFeature = 'areaVtPulse';
cat4_queryRunFeature = 'areaVtPulse';
cat5_queryRunFeature = 'areaVtPulse';
cat6_queryRunFeature = 'areaVtPulse';
cat7_queryRunFeature = 'areaVtPulse';
cat8_queryRunFeature = 'areaVtPulse';
cat9_queryRunFeature = 'areaVtPulse';

cat1_selectedCellsIDs = 5;
cat2_selectedCellsIDs = 5;
cat3_selectedCellsIDs = 5;
cat4_selectedCellsIDs = 5;
cat5_selectedCellsIDs = 5;
cat6_selectedCellsIDs = 5;
cat7_selectedCellsIDs = 5;
cat8_selectedCellsIDs = 5;
cat9_selectedCellsIDs = 5;

cat1_desiredEpochs = [];
cat2_desiredEpochs = [];
cat3_desiredEpochs = [];
cat4_desiredEpochs = [];
cat5_desiredEpochs = [];
cat6_desiredEpochs = [];
cat7_desiredEpochs = [];
cat8_desiredEpochs = [];
cat9_desiredEpochs = [];

cat1_selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0') & cell2mat(runGroupsParametersTable.amplitudeCurrentPulse) == -100 & cell2mat(runGroupsParametersTable.widthCurrentPulse) == 6000 & cell2mat(runGroupsParametersTable.widthLightPulse) == 1000 & cell2mat(runGroupsParametersTable.amplitudeLightPulse) == 0 & cell2mat(runGroupsParametersTable.delayLightPulse) == 25500);
cat2_selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0') & cell2mat(runGroupsParametersTable.amplitudeCurrentPulse) == -50 & cell2mat(runGroupsParametersTable.widthCurrentPulse) == 6000 & cell2mat(runGroupsParametersTable.widthLightPulse) == 1000 & cell2mat(runGroupsParametersTable.amplitudeLightPulse) == 0 & cell2mat(runGroupsParametersTable.delayLightPulse) == 25500);
cat3_selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0') & cell2mat(runGroupsParametersTable.amplitudeCurrentPulse) == 0 & cell2mat(runGroupsParametersTable.widthCurrentPulse) == 6000 & cell2mat(runGroupsParametersTable.widthLightPulse) == 1000 & cell2mat(runGroupsParametersTable.amplitudeLightPulse) == 0 & cell2mat(runGroupsParametersTable.delayLightPulse) == 25500);
cat4_selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0') & cell2mat(runGroupsParametersTable.amplitudeCurrentPulse) == 50 & cell2mat(runGroupsParametersTable.widthCurrentPulse) == 6000 & cell2mat(runGroupsParametersTable.widthLightPulse) == 1000 & cell2mat(runGroupsParametersTable.amplitudeLightPulse) == 0 & cell2mat(runGroupsParametersTable.delayLightPulse) == 25500);
cat5_selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0') & cell2mat(runGroupsParametersTable.amplitudeCurrentPulse) == 100 & cell2mat(runGroupsParametersTable.widthCurrentPulse) == 6000 & cell2mat(runGroupsParametersTable.widthLightPulse) == 1000 & cell2mat(runGroupsParametersTable.amplitudeLightPulse) == 0 & cell2mat(runGroupsParametersTable.delayLightPulse) == 25500);
cat6_selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0') & cell2mat(runGroupsParametersTable.amplitudeCurrentPulse) == 250 & cell2mat(runGroupsParametersTable.widthCurrentPulse) == 6000 & cell2mat(runGroupsParametersTable.widthLightPulse) == 1000 & cell2mat(runGroupsParametersTable.amplitudeLightPulse) == 0 & cell2mat(runGroupsParametersTable.delayLightPulse) == 25500);
cat7_selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0') & cell2mat(runGroupsParametersTable.amplitudeCurrentPulse) == 500 & cell2mat(runGroupsParametersTable.widthCurrentPulse) == 6000 & cell2mat(runGroupsParametersTable.widthLightPulse) == 1000 & cell2mat(runGroupsParametersTable.amplitudeLightPulse) == 0 & cell2mat(runGroupsParametersTable.delayLightPulse) == 25500);
cat8_selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0') & cell2mat(runGroupsParametersTable.amplitudeCurrentPulse) == 750 & cell2mat(runGroupsParametersTable.widthCurrentPulse) == 6000 & cell2mat(runGroupsParametersTable.widthLightPulse) == 1000 & cell2mat(runGroupsParametersTable.amplitudeLightPulse) == 0 & cell2mat(runGroupsParametersTable.delayLightPulse) == 25500);
cat9_selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0') & cell2mat(runGroupsParametersTable.amplitudeCurrentPulse) == 1000 & cell2mat(runGroupsParametersTable.widthCurrentPulse) == 6000 & cell2mat(runGroupsParametersTable.widthLightPulse) == 1000 & cell2mat(runGroupsParametersTable.amplitudeLightPulse) == 0 & cell2mat(runGroupsParametersTable.delayLightPulse) == 25500);

% No need to modify
cat1_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat1_selectedCellsIDs);
cat1_selectedRunGroupsAnalysisTable = cat1_selectedMembersTable(ismember(cat1_selectedMembersTable.runGroup,cat1_selectedRunGroups),:);
cat1_selectedRunGroupsAnalysisTable = extractEpochSubset(cat1_selectedRunGroupsAnalysisTable, cat1_desiredEpochs);

cat2_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat2_selectedCellsIDs);
cat2_selectedRunGroupsAnalysisTable = cat2_selectedMembersTable(ismember(cat2_selectedMembersTable.runGroup,cat2_selectedRunGroups),:);
cat2_selectedRunGroupsAnalysisTable = extractEpochSubset(cat2_selectedRunGroupsAnalysisTable, cat2_desiredEpochs);

cat3_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat3_selectedCellsIDs);
cat3_selectedRunGroupsAnalysisTable = cat3_selectedMembersTable(ismember(cat3_selectedMembersTable.runGroup,cat3_selectedRunGroups),:);
cat3_selectedRunGroupsAnalysisTable = extractEpochSubset(cat3_selectedRunGroupsAnalysisTable, cat3_desiredEpochs);

cat4_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat4_selectedCellsIDs);
cat4_selectedRunGroupsAnalysisTable = cat4_selectedMembersTable(ismember(cat4_selectedMembersTable.runGroup,cat4_selectedRunGroups),:);
cat4_selectedRunGroupsAnalysisTable = extractEpochSubset(cat4_selectedRunGroupsAnalysisTable, cat4_desiredEpochs);

cat5_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat5_selectedCellsIDs);
cat5_selectedRunGroupsAnalysisTable = cat5_selectedMembersTable(ismember(cat5_selectedMembersTable.runGroup,cat5_selectedRunGroups),:);
cat5_selectedRunGroupsAnalysisTable = extractEpochSubset(cat5_selectedRunGroupsAnalysisTable, cat5_desiredEpochs);

cat6_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat6_selectedCellsIDs);
cat6_selectedRunGroupsAnalysisTable = cat6_selectedMembersTable(ismember(cat6_selectedMembersTable.runGroup,cat6_selectedRunGroups),:);
cat6_selectedRunGroupsAnalysisTable = extractEpochSubset(cat6_selectedRunGroupsAnalysisTable, cat6_desiredEpochs);

cat7_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat7_selectedCellsIDs);
cat7_selectedRunGroupsAnalysisTable = cat7_selectedMembersTable(ismember(cat7_selectedMembersTable.runGroup,cat7_selectedRunGroups),:);
cat7_selectedRunGroupsAnalysisTable = extractEpochSubset(cat7_selectedRunGroupsAnalysisTable, cat7_desiredEpochs);

cat8_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat8_selectedCellsIDs);
cat8_selectedRunGroupsAnalysisTable = cat8_selectedMembersTable(ismember(cat8_selectedMembersTable.runGroup,cat8_selectedRunGroups),:);
cat8_selectedRunGroupsAnalysisTable = extractEpochSubset(cat8_selectedRunGroupsAnalysisTable, cat8_desiredEpochs);

cat9_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat9_selectedCellsIDs);
cat9_selectedRunGroupsAnalysisTable = cat9_selectedMembersTable(ismember(cat9_selectedMembersTable.runGroup,cat9_selectedRunGroups),:);
cat9_selectedRunGroupsAnalysisTable = extractEpochSubset(cat9_selectedRunGroupsAnalysisTable, cat9_desiredEpochs);

plotData = generatePlotDataStruct;
plotData.queryRunFeature.cat1 = cat1_queryRunFeature;
plotData.queryRunFeature.cat2 = cat2_queryRunFeature;
plotData.queryRunFeature.cat3 = cat3_queryRunFeature;
plotData.queryRunFeature.cat4 = cat4_queryRunFeature;
plotData.queryRunFeature.cat5 = cat5_queryRunFeature;
plotData.queryRunFeature.cat6 = cat6_queryRunFeature;
plotData.queryRunFeature.cat7 = cat7_queryRunFeature;
plotData.queryRunFeature.cat8 = cat8_queryRunFeature;
plotData.queryRunFeature.cat9 = cat9_queryRunFeature;
plotData.selectedRunGroupsAnalysisTable.cat1 = cat1_selectedRunGroupsAnalysisTable;
plotData.selectedRunGroupsAnalysisTable.cat2 = cat2_selectedRunGroupsAnalysisTable;
plotData.selectedRunGroupsAnalysisTable.cat3 = cat3_selectedRunGroupsAnalysisTable;
plotData.selectedRunGroupsAnalysisTable.cat4 = cat4_selectedRunGroupsAnalysisTable;
plotData.selectedRunGroupsAnalysisTable.cat5 = cat5_selectedRunGroupsAnalysisTable;
plotData.selectedRunGroupsAnalysisTable.cat6 = cat6_selectedRunGroupsAnalysisTable;
plotData.selectedRunGroupsAnalysisTable.cat7 = cat7_selectedRunGroupsAnalysisTable;
plotData.selectedRunGroupsAnalysisTable.cat8 = cat8_selectedRunGroupsAnalysisTable;
plotData.selectedRunGroupsAnalysisTable.cat9 = cat9_selectedRunGroupsAnalysisTable;

%% Plot
nexttile(masterLayout,1,[1 1]); axis off;

options = generateOptionsStruct(groupByCell);
options.plotColor.cat1 = colorGrey;
options.plotColor.cat2 = colorGrey;
options.plotColor.cat3 = colorGrey;
options.plotColor.cat4 = colorGrey;
options.plotColor.cat5 = colorGrey;
options.plotColor.cat6 = colorGrey;
options.plotColor.cat7 = colorGrey;
options.plotColor.cat8 = colorGrey;
options.plotColor.cat9 = colorGrey;
options.plotColorOverlay.cat1 = colorGreyOverlay;
options.plotColorOverlay.cat2 = colorGreyOverlay;
options.plotColorOverlay.cat3 = colorGreyOverlay;
options.plotColorOverlay.cat4 = colorGreyOverlay;
options.plotColorOverlay.cat5 = colorGreyOverlay;
options.plotColorOverlay.cat6 = colorGreyOverlay;
options.plotColorOverlay.cat7 = colorGreyOverlay;
options.plotColorOverlay.cat8 = colorGreyOverlay;
options.plotColorOverlay.cat9 = colorGreyOverlay;
options.plotTitle = ['Input-output curve'];
options.xLabelText.cat1 = ['-100 pA'];
options.xLabelText.cat2 = ['-50 pA'];
options.xLabelText.cat3 = ['0 pA'];
options.xLabelText.cat4 = ['50 pA'];
options.xLabelText.cat5 = ['100 pA'];
options.xLabelText.cat6 = ['250 pA'];
options.xLabelText.cat7 = ['500 pA'];
options.xLabelText.cat8 = ['750 pA'];
options.xLabelText.cat9 = ['1000 pA'];
options.yLabelText = ['AUC [mVs]'];

generateLinePlot(plotData,options)

%% Collect data

cat1_queryRunFeature = 'areaVtPulse';
cat2_queryRunFeature = 'areaVtPulse';
cat3_queryRunFeature = 'areaVtPulse';
cat4_queryRunFeature = 'areaVtPulse';
cat5_queryRunFeature = 'areaVtPulse';
cat6_queryRunFeature = 'areaVtPulse';
cat7_queryRunFeature = 'areaVtPulse';
cat8_queryRunFeature = 'areaVtPulse';
cat9_queryRunFeature = 'areaVtPulse';

cat1_selectedCellsIDs = 5;
cat2_selectedCellsIDs = 5;
cat3_selectedCellsIDs = 5;
cat4_selectedCellsIDs = 5;
cat5_selectedCellsIDs = 5;
cat6_selectedCellsIDs = 5;
cat7_selectedCellsIDs = 5;
cat8_selectedCellsIDs = 5;
cat9_selectedCellsIDs = 5;

cat1_desiredEpochs = [];
cat2_desiredEpochs = [];
cat3_desiredEpochs = [];
cat4_desiredEpochs = [];
cat5_desiredEpochs = [];
cat6_desiredEpochs = [];
cat7_desiredEpochs = [];
cat8_desiredEpochs = [];
cat9_desiredEpochs = [];

cat1_selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.amplitudeCurrentPulse) == -100 & cell2mat(runGroupsParametersTable.widthCurrentPulse) == 6000 & cell2mat(runGroupsParametersTable.widthLightPulse) == 12000 & cell2mat(runGroupsParametersTable.amplitudeLightPulse) == 5 & cell2mat(runGroupsParametersTable.delayLightPulse) == 19000);
cat2_selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.amplitudeCurrentPulse) == -50 & cell2mat(runGroupsParametersTable.widthCurrentPulse) == 6000 & cell2mat(runGroupsParametersTable.widthLightPulse) == 12000 & cell2mat(runGroupsParametersTable.amplitudeLightPulse) == 5 & cell2mat(runGroupsParametersTable.delayLightPulse) == 19000);
cat3_selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.amplitudeCurrentPulse) == 0 & cell2mat(runGroupsParametersTable.widthCurrentPulse) == 6000 & cell2mat(runGroupsParametersTable.widthLightPulse) == 12000 & cell2mat(runGroupsParametersTable.amplitudeLightPulse) == 5 & cell2mat(runGroupsParametersTable.delayLightPulse) == 19000);
cat4_selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.amplitudeCurrentPulse) == 50 & cell2mat(runGroupsParametersTable.widthCurrentPulse) == 6000 & cell2mat(runGroupsParametersTable.widthLightPulse) == 12000 & cell2mat(runGroupsParametersTable.amplitudeLightPulse) == 5 & cell2mat(runGroupsParametersTable.delayLightPulse) == 19000);
cat5_selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.amplitudeCurrentPulse) == 100 & cell2mat(runGroupsParametersTable.widthCurrentPulse) == 6000 & cell2mat(runGroupsParametersTable.widthLightPulse) == 12000 & cell2mat(runGroupsParametersTable.amplitudeLightPulse) == 5 & cell2mat(runGroupsParametersTable.delayLightPulse) == 19000);
cat6_selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.amplitudeCurrentPulse) == 250 & cell2mat(runGroupsParametersTable.widthCurrentPulse) == 6000 & cell2mat(runGroupsParametersTable.widthLightPulse) == 12000 & cell2mat(runGroupsParametersTable.amplitudeLightPulse) == 5 & cell2mat(runGroupsParametersTable.delayLightPulse) == 19000);
cat7_selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.amplitudeCurrentPulse) == 500 & cell2mat(runGroupsParametersTable.widthCurrentPulse) == 6000 & cell2mat(runGroupsParametersTable.widthLightPulse) == 12000 & cell2mat(runGroupsParametersTable.amplitudeLightPulse) == 5 & cell2mat(runGroupsParametersTable.delayLightPulse) == 19000);
cat8_selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.amplitudeCurrentPulse) == 750 & cell2mat(runGroupsParametersTable.widthCurrentPulse) == 6000 & cell2mat(runGroupsParametersTable.widthLightPulse) == 12000 & cell2mat(runGroupsParametersTable.amplitudeLightPulse) == 5 & cell2mat(runGroupsParametersTable.delayLightPulse) == 19000);
cat9_selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.amplitudeCurrentPulse) == 1000 & cell2mat(runGroupsParametersTable.widthCurrentPulse) == 6000 & cell2mat(runGroupsParametersTable.widthLightPulse) == 12000 & cell2mat(runGroupsParametersTable.amplitudeLightPulse) == 5 & cell2mat(runGroupsParametersTable.delayLightPulse) == 19000);

% No need to modify
cat1_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat1_selectedCellsIDs);
cat1_selectedRunGroupsAnalysisTable = cat1_selectedMembersTable(ismember(cat1_selectedMembersTable.runGroup,cat1_selectedRunGroups),:);
cat1_selectedRunGroupsAnalysisTable = extractEpochSubset(cat1_selectedRunGroupsAnalysisTable, cat1_desiredEpochs);

cat2_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat2_selectedCellsIDs);
cat2_selectedRunGroupsAnalysisTable = cat2_selectedMembersTable(ismember(cat2_selectedMembersTable.runGroup,cat2_selectedRunGroups),:);
cat2_selectedRunGroupsAnalysisTable = extractEpochSubset(cat2_selectedRunGroupsAnalysisTable, cat2_desiredEpochs);

cat3_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat3_selectedCellsIDs);
cat3_selectedRunGroupsAnalysisTable = cat3_selectedMembersTable(ismember(cat3_selectedMembersTable.runGroup,cat3_selectedRunGroups),:);
cat3_selectedRunGroupsAnalysisTable = extractEpochSubset(cat3_selectedRunGroupsAnalysisTable, cat3_desiredEpochs);

cat4_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat4_selectedCellsIDs);
cat4_selectedRunGroupsAnalysisTable = cat4_selectedMembersTable(ismember(cat4_selectedMembersTable.runGroup,cat4_selectedRunGroups),:);
cat4_selectedRunGroupsAnalysisTable = extractEpochSubset(cat4_selectedRunGroupsAnalysisTable, cat4_desiredEpochs);

cat5_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat5_selectedCellsIDs);
cat5_selectedRunGroupsAnalysisTable = cat5_selectedMembersTable(ismember(cat5_selectedMembersTable.runGroup,cat5_selectedRunGroups),:);
cat5_selectedRunGroupsAnalysisTable = extractEpochSubset(cat5_selectedRunGroupsAnalysisTable, cat5_desiredEpochs);

cat6_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat6_selectedCellsIDs);
cat6_selectedRunGroupsAnalysisTable = cat6_selectedMembersTable(ismember(cat6_selectedMembersTable.runGroup,cat6_selectedRunGroups),:);
cat6_selectedRunGroupsAnalysisTable = extractEpochSubset(cat6_selectedRunGroupsAnalysisTable, cat6_desiredEpochs);

cat7_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat7_selectedCellsIDs);
cat7_selectedRunGroupsAnalysisTable = cat7_selectedMembersTable(ismember(cat7_selectedMembersTable.runGroup,cat7_selectedRunGroups),:);
cat7_selectedRunGroupsAnalysisTable = extractEpochSubset(cat7_selectedRunGroupsAnalysisTable, cat7_desiredEpochs);

cat8_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat8_selectedCellsIDs);
cat8_selectedRunGroupsAnalysisTable = cat8_selectedMembersTable(ismember(cat8_selectedMembersTable.runGroup,cat8_selectedRunGroups),:);
cat8_selectedRunGroupsAnalysisTable = extractEpochSubset(cat8_selectedRunGroupsAnalysisTable, cat8_desiredEpochs);

cat9_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat9_selectedCellsIDs);
cat9_selectedRunGroupsAnalysisTable = cat9_selectedMembersTable(ismember(cat9_selectedMembersTable.runGroup,cat9_selectedRunGroups),:);
cat9_selectedRunGroupsAnalysisTable = extractEpochSubset(cat9_selectedRunGroupsAnalysisTable, cat9_desiredEpochs);

plotData = generatePlotDataStruct;
plotData.queryRunFeature.cat1 = cat1_queryRunFeature;
plotData.queryRunFeature.cat2 = cat2_queryRunFeature;
plotData.queryRunFeature.cat3 = cat3_queryRunFeature;
plotData.queryRunFeature.cat4 = cat4_queryRunFeature;
plotData.queryRunFeature.cat5 = cat5_queryRunFeature;
plotData.queryRunFeature.cat6 = cat6_queryRunFeature;
plotData.queryRunFeature.cat7 = cat7_queryRunFeature;
plotData.queryRunFeature.cat8 = cat8_queryRunFeature;
plotData.queryRunFeature.cat9 = cat9_queryRunFeature;
plotData.selectedRunGroupsAnalysisTable.cat1 = cat1_selectedRunGroupsAnalysisTable;
plotData.selectedRunGroupsAnalysisTable.cat2 = cat2_selectedRunGroupsAnalysisTable;
plotData.selectedRunGroupsAnalysisTable.cat3 = cat3_selectedRunGroupsAnalysisTable;
plotData.selectedRunGroupsAnalysisTable.cat4 = cat4_selectedRunGroupsAnalysisTable;
plotData.selectedRunGroupsAnalysisTable.cat5 = cat5_selectedRunGroupsAnalysisTable;
plotData.selectedRunGroupsAnalysisTable.cat6 = cat6_selectedRunGroupsAnalysisTable;
plotData.selectedRunGroupsAnalysisTable.cat7 = cat7_selectedRunGroupsAnalysisTable;
plotData.selectedRunGroupsAnalysisTable.cat8 = cat8_selectedRunGroupsAnalysisTable;
plotData.selectedRunGroupsAnalysisTable.cat9 = cat9_selectedRunGroupsAnalysisTable;

%% Plot
hold on;

options = generateOptionsStruct(groupByCell);
options.plotColor.cat1 = colorBlue;
options.plotColor.cat2 = colorBlue;
options.plotColor.cat3 = colorBlue;
options.plotColor.cat4 = colorBlue;
options.plotColor.cat5 = colorBlue;
options.plotColor.cat6 = colorBlue;
options.plotColor.cat7 = colorBlue;
options.plotColor.cat8 = colorBlue;
options.plotColor.cat9 = colorBlue;
options.plotColorOverlay.cat1 = colorBlueOverlay;
options.plotColorOverlay.cat2 = colorBlueOverlay;
options.plotColorOverlay.cat3 = colorBlueOverlay;
options.plotColorOverlay.cat4 = colorBlueOverlay;
options.plotColorOverlay.cat5 = colorBlueOverlay;
options.plotColorOverlay.cat6 = colorBlueOverlay;
options.plotColorOverlay.cat7 = colorBlueOverlay;
options.plotColorOverlay.cat8 = colorBlueOverlay;
options.plotColorOverlay.cat9 = colorBlueOverlay;
options.plotTitle = ['Input-output curve'];
options.xLabelText.cat1 = ['-100 pA'];
options.xLabelText.cat2 = ['-50 pA'];
options.xLabelText.cat3 = ['0 pA'];
options.xLabelText.cat4 = ['50 pA'];
options.xLabelText.cat5 = ['100 pA'];
options.xLabelText.cat6 = ['250 pA'];
options.xLabelText.cat7 = ['500 pA'];
options.xLabelText.cat8 = ['750 pA'];
options.xLabelText.cat9 = ['1000 pA'];
options.yLabelText = ['AUC [mVs]'];
options.xLimits = [0.5,9.5];

generateLinePlot(plotData,options)

end