disp(['*-*-*-* Running: scriptPlotAnalysisCellInhibitoryOpsinTesting *-*-*-*'])

mouseAnalysisTablePath = [experimentDirectory filesep 'mouseAnalysis' filesep 'MouseAnalysisTable.mat'];
mouseParametersTablePath = [experimentDirectory filesep 'mouseAnalysis' filesep 'MouseParametersTable.mat'];

mouseAnalysisTable = load(mouseAnalysisTablePath);
mouseAnalysisTable = mouseAnalysisTable.mouseAnalysisTable;
mouseParametersTable = load(mouseParametersTablePath);
runGroupsParametersTable = mouseParametersTable.runGroupsParametersTable;

postTextSuptitle = '';
analyzedCells = unique(cell2mat(mouseAnalysisTable.cellName));

groupByCell = 0;

for iCell = 1:numel(analyzedCells)
    
iAnalyzedCell = analyzedCells(iCell);

cellName = ['cell' num2str(iAnalyzedCell)];

initializeFig(0.9,0.9);

nTileRows = 2;
nTileColumns = 4;
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

%% ** Plot blue quality RS0 **

%% Collect data

cat1_queryRunFeature = 'qualityRS0';
cat2_queryRunFeature = 'qualityRS0';
cat3_queryRunFeature = 'qualityRS0';

cat1_selectedCellsIDs = iAnalyzedCell;
cat2_selectedCellsIDs = iAnalyzedCell;
cat3_selectedCellsIDs = iAnalyzedCell;

cat1_desiredEpochs = [];
cat2_desiredEpochs = [];
cat3_desiredEpochs = [];

cat1_queryParameters = [];
cat2_queryParameters = [];
cat3_queryParameters = [];

cat1_queryParameters.amplitudeBlue = 1;
cat1_queryParameters.pulseWidthBlue = 5000;
cat1_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'blue','inhibitory',cat1_queryParameters);

cat2_queryParameters.amplitudeBlue = 3;
cat2_queryParameters.pulseWidthBlue = 5000;
cat2_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'blue','inhibitory',cat2_queryParameters);

cat3_queryParameters.amplitudeBlue = 5;
cat3_queryParameters.pulseWidthBlue = 5000;
cat3_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'blue','inhibitory',cat3_queryParameters);


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

plotData = generatePlotDataStruct;
plotData.queryRunFeature.cat1 = cat1_queryRunFeature;
plotData.queryRunFeature.cat2 = cat2_queryRunFeature;
plotData.queryRunFeature.cat3 = cat3_queryRunFeature;
plotData.selectedRunGroupsAnalysisTable.cat1 = cat1_selectedRunGroupsAnalysisTable;
plotData.selectedRunGroupsAnalysisTable.cat2 = cat2_selectedRunGroupsAnalysisTable;
plotData.selectedRunGroupsAnalysisTable.cat3 = cat3_selectedRunGroupsAnalysisTable;

%% Plot
nexttile(masterLayout,1,[1 1]); axis off;

options = generateOptionsStruct(groupByCell);
options.plotColor.cat1 = colorBlue;
options.plotColor.cat2 = colorBlue*0.85;
options.plotColor.cat3 = colorBlue*0.7;
options.plotTitle = ['Recording quality'];
options.xLabelText.cat1 = ['Acquisition'];
options.xLabelText.cat2 = ['Acquisition'];
options.xLabelText.cat3 = ['Acquisition'];
options.yLabelText = ['Series resistance [M\Omega]'];
options.legendText.cat1 = ['1 mV'];
options.legendText.cat2 = ['3 mV'];
options.legendText.cat3 = ['5 mV'];

generateDotplot(plotData,options)
    
%% ** Plot blue traces **

%% Collect data

cat1_selectedCellsIDs = iAnalyzedCell;
cat2_selectedCellsIDs = iAnalyzedCell;
cat3_selectedCellsIDs = iAnalyzedCell;

cat1_desiredEpochs = [];
cat2_desiredEpochs = [];
cat3_desiredEpochs = [];

cat1_queryParameters = [];
cat2_queryParameters = [];
cat3_queryParameters = [];

cat1_queryParameters.amplitudeBlue = 1;
cat1_queryParameters.pulseWidthBlue = 5000;
cat1_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'blue','inhibitory',cat1_queryParameters);

cat2_queryParameters.amplitudeBlue = 3;
cat2_queryParameters.pulseWidthBlue = 5000;
cat2_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'blue','inhibitory',cat2_queryParameters);

cat3_queryParameters.amplitudeBlue = 5;
cat3_queryParameters.pulseWidthBlue = 5000;
cat3_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'blue','inhibitory',cat3_queryParameters);


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

plotData = generatePlotDataStruct;
plotData.queryRunFeature.cat1 = cat1_queryRunFeature;
plotData.queryRunFeature.cat2 = cat2_queryRunFeature;
plotData.queryRunFeature.cat3 = cat3_queryRunFeature;
plotData.selectedRunGroupsAnalysisTable.cat1 = cat1_selectedRunGroupsAnalysisTable;
plotData.selectedRunGroupsAnalysisTable.cat2 = cat2_selectedRunGroupsAnalysisTable;
plotData.selectedRunGroupsAnalysisTable.cat3 = cat3_selectedRunGroupsAnalysisTable;

%% Plot
nexttile(masterLayout,2,[1 1]); axis off;

options = generateOptionsStruct(groupByCell);
options.plotColor.cat1 = colorBlue;
options.plotColor.cat2 = colorBlue*0.85;
options.plotColor.cat3 = colorBlue*0.7;
options.plotTitle = ['Full-field responses'];
options.legendText.cat1 = ['1 mV'];
options.legendText.cat2 = ['3 mV'];
options.legendText.cat3 = ['5 mV'];
options.xLimits = [0,20000];

generatePlotTraces(plotData,options)

%% ** Barplot blue heightPulsePeak different power **

%% Collect data

cat1_queryRunFeature = 'heightPulsePeak';
cat2_queryRunFeature = 'heightPulsePeak';
cat3_queryRunFeature = 'heightPulsePeak';

cat1_selectedCellsIDs = iAnalyzedCell;
cat2_selectedCellsIDs = iAnalyzedCell;
cat3_selectedCellsIDs = iAnalyzedCell;

cat1_desiredEpochs = [];
cat2_desiredEpochs = [];
cat3_desiredEpochs = [];

cat1_queryParameters = [];
cat2_queryParameters = [];
cat3_queryParameters = [];

cat1_queryParameters.amplitudeBlue = 1;
cat1_queryParameters.pulseWidthBlue = 5000;
cat1_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'blue','inhibitory',cat1_queryParameters);

cat2_queryParameters.amplitudeBlue = 3;
cat2_queryParameters.pulseWidthBlue = 5000;
cat2_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'blue','inhibitory',cat2_queryParameters);

cat3_queryParameters.amplitudeBlue = 5;
cat3_queryParameters.pulseWidthBlue = 5000;
cat3_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'blue','inhibitory',cat3_queryParameters);


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

plotData = generatePlotDataStruct;
plotData.queryRunFeature.cat1 = cat1_queryRunFeature;
plotData.queryRunFeature.cat2 = cat2_queryRunFeature;
plotData.queryRunFeature.cat3 = cat3_queryRunFeature;
plotData.selectedRunGroupsAnalysisTable.cat1 = cat1_selectedRunGroupsAnalysisTable;
plotData.selectedRunGroupsAnalysisTable.cat2 = cat2_selectedRunGroupsAnalysisTable;
plotData.selectedRunGroupsAnalysisTable.cat3 = cat3_selectedRunGroupsAnalysisTable;

%% Plot
nexttile(masterLayout,3,[1 1]); axis off;

options = generateOptionsStruct(groupByCell);
options.plotColor.cat1 = colorBlue;
options.plotColor.cat2 = colorBlue*0.85;
options.plotColor.cat3 = colorBlue*0.7;
options.plotColorOverlay.cat1 = colorBlueOverlay;
options.plotColorOverlay.cat2 = colorBlueOverlay*0.85;
options.plotColorOverlay.cat3 = colorBlueOverlay*0.7;
options.plotTitle = ['Response amplitude'];
options.xLabelText.cat1 = ['1 mV'];
options.xLabelText.cat2 = ['3 mV'];
options.xLabelText.cat3 = ['5 mV'];
options.yLabelText = ['Current [pA]'];

generateBarplot(plotData,options)

%% ** Barplot blue areaPulse different power **

%% Collect data

cat1_queryRunFeature = 'areaPulse';
cat2_queryRunFeature = 'areaPulse';
cat3_queryRunFeature = 'areaPulse';

cat1_selectedCellsIDs = iAnalyzedCell;
cat2_selectedCellsIDs = iAnalyzedCell;
cat3_selectedCellsIDs = iAnalyzedCell;

cat1_desiredEpochs = [];
cat2_desiredEpochs = [];
cat3_desiredEpochs = [];

cat1_queryParameters = [];
cat2_queryParameters = [];
cat3_queryParameters = [];

cat1_queryParameters.amplitudeBlue = 1;
cat1_queryParameters.pulseWidthBlue = 5000;
cat1_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'blue','inhibitory',cat1_queryParameters);

cat2_queryParameters.amplitudeBlue = 3;
cat2_queryParameters.pulseWidthBlue = 5000;
cat2_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'blue','inhibitory',cat2_queryParameters);

cat3_queryParameters.amplitudeBlue = 5;
cat3_queryParameters.pulseWidthBlue = 5000;
cat3_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'blue','inhibitory',cat3_queryParameters);


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

plotData = generatePlotDataStruct;
plotData.queryRunFeature.cat1 = cat1_queryRunFeature;
plotData.queryRunFeature.cat2 = cat2_queryRunFeature;
plotData.queryRunFeature.cat3 = cat3_queryRunFeature;
plotData.selectedRunGroupsAnalysisTable.cat1 = cat1_selectedRunGroupsAnalysisTable;
plotData.selectedRunGroupsAnalysisTable.cat2 = cat2_selectedRunGroupsAnalysisTable;
plotData.selectedRunGroupsAnalysisTable.cat3 = cat3_selectedRunGroupsAnalysisTable;

%% Plot
nexttile(masterLayout,4,[1 1]); axis off;

options = generateOptionsStruct(groupByCell);
options.plotColor.cat1 = colorBlue;
options.plotColor.cat2 = colorBlue*0.85;
options.plotColor.cat3 = colorBlue*0.7;
options.plotColorOverlay.cat1 = colorBlueOverlay;
options.plotColorOverlay.cat2 = colorBlueOverlay*0.85;
options.plotColorOverlay.cat3 = colorBlueOverlay*0.7;
options.plotTitle = ['Response charge'];
options.xLabelText.cat1 = ['1 mV'];
options.xLabelText.cat2 = ['3 mV'];
options.xLabelText.cat3 = ['5 mV'];
options.yLabelText = ['Charge [pC]'];

generateBarplot(plotData,options)

%% ** Plot blue quality RS0 **

%% Collect data

cat1_queryRunFeature = 'qualityRS0';
cat2_queryRunFeature = 'qualityRS0';
cat3_queryRunFeature = 'qualityRS0';
cat4_queryRunFeature = 'qualityRS0';
cat5_queryRunFeature = 'qualityRS0';

cat1_selectedCellsIDs = iAnalyzedCell;
cat2_selectedCellsIDs = iAnalyzedCell;
cat3_selectedCellsIDs = iAnalyzedCell;
cat4_selectedCellsIDs = iAnalyzedCell;
cat5_selectedCellsIDs = iAnalyzedCell;

cat1_desiredEpochs = [];
cat2_desiredEpochs = [];
cat3_desiredEpochs = [];
cat4_desiredEpochs = [];
cat5_desiredEpochs = [];

cat1_queryParameters = [];
cat2_queryParameters = [];
cat3_queryParameters = [];
cat4_queryParameters = [];
cat5_queryParameters = [];

cat1_queryParameters.amplitudeBlue = 5;
cat1_queryParameters.pulseWidthBlue = 250;
cat1_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'blue','inhibitory',cat1_queryParameters);

cat2_queryParameters.amplitudeBlue = 5;
cat2_queryParameters.pulseWidthBlue = 500;
cat2_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'blue','inhibitory',cat2_queryParameters);

cat3_queryParameters.amplitudeBlue = 5;
cat3_queryParameters.pulseWidthBlue = 2000;
cat3_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'blue','inhibitory',cat3_queryParameters);

cat4_queryParameters.amplitudeBlue = 5;
cat4_queryParameters.pulseWidthBlue = 5000;
cat4_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'blue','inhibitory',cat4_queryParameters);

cat5_queryParameters.amplitudeBlue = 5;
cat5_queryParameters.pulseWidthBlue = 10000;
cat5_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'blue','inhibitory',cat5_queryParameters);


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

plotData = generatePlotDataStruct;
plotData.queryRunFeature.cat1 = cat1_queryRunFeature;
plotData.queryRunFeature.cat2 = cat2_queryRunFeature;
plotData.queryRunFeature.cat3 = cat3_queryRunFeature;
plotData.queryRunFeature.cat4 = cat4_queryRunFeature;
plotData.queryRunFeature.cat5 = cat5_queryRunFeature;
plotData.selectedRunGroupsAnalysisTable.cat1 = cat1_selectedRunGroupsAnalysisTable;
plotData.selectedRunGroupsAnalysisTable.cat2 = cat2_selectedRunGroupsAnalysisTable;
plotData.selectedRunGroupsAnalysisTable.cat3 = cat3_selectedRunGroupsAnalysisTable;
plotData.selectedRunGroupsAnalysisTable.cat4 = cat4_selectedRunGroupsAnalysisTable;
plotData.selectedRunGroupsAnalysisTable.cat5 = cat5_selectedRunGroupsAnalysisTable;

%% Plot
nexttile(masterLayout,5,[1 1]); axis off;

options = generateOptionsStruct(groupByCell);
options.plotColor.cat1 = colorBlue;
options.plotColor.cat2 = colorBlue*0.9;
options.plotColor.cat3 = colorBlue*0.8;
options.plotColor.cat4 = colorBlue*0.7;
options.plotColor.cat5 = colorBlue*0.6;
options.plotTitle = ['Recording quality'];
options.xLabelText.cat1 = ['Acquisition'];
options.xLabelText.cat2 = ['Acquisition'];
options.xLabelText.cat3 = ['Acquisition'];
options.xLabelText.cat4 = ['Acquisition'];
options.xLabelText.cat5 = ['Acquisition'];
options.yLabelText = ['Series resistance [M\Omega]'];
options.legendText.cat1 = ['25 ms'];
options.legendText.cat2 = ['50 ms'];
options.legendText.cat3 = ['200 ms'];
options.legendText.cat4 = ['500 ms'];
options.legendText.cat5 = ['1000 ms'];

generateDotplot(plotData,options)
    
%% ** Plot blue traces **

%% Collect data

cat1_selectedCellsIDs = iAnalyzedCell;
cat2_selectedCellsIDs = iAnalyzedCell;
cat3_selectedCellsIDs = iAnalyzedCell;
cat4_selectedCellsIDs = iAnalyzedCell;
cat5_selectedCellsIDs = iAnalyzedCell;

cat1_desiredEpochs = [];
cat2_desiredEpochs = [];
cat3_desiredEpochs = [];
cat4_desiredEpochs = [];
cat5_desiredEpochs = [];

cat1_queryParameters = [];
cat2_queryParameters = [];
cat3_queryParameters = [];
cat4_queryParameters = [];
cat5_queryParameters = [];

cat1_queryParameters.amplitudeBlue = 5;
cat1_queryParameters.pulseWidthBlue = 250;
cat1_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'blue','inhibitory',cat1_queryParameters);

cat2_queryParameters.amplitudeBlue = 5;
cat2_queryParameters.pulseWidthBlue = 500;
cat2_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'blue','inhibitory',cat2_queryParameters);

cat3_queryParameters.amplitudeBlue = 5;
cat3_queryParameters.pulseWidthBlue = 2000;
cat3_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'blue','inhibitory',cat3_queryParameters);

cat4_queryParameters.amplitudeBlue = 5;
cat4_queryParameters.pulseWidthBlue = 5000;
cat4_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'blue','inhibitory',cat4_queryParameters);

cat5_queryParameters.amplitudeBlue = 5;
cat5_queryParameters.pulseWidthBlue = 10000;
cat5_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'blue','inhibitory',cat5_queryParameters);


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

plotData = generatePlotDataStruct;
plotData.queryRunFeature.cat1 = cat1_queryRunFeature;
plotData.queryRunFeature.cat2 = cat2_queryRunFeature;
plotData.queryRunFeature.cat3 = cat3_queryRunFeature;
plotData.queryRunFeature.cat4 = cat4_queryRunFeature;
plotData.queryRunFeature.cat5 = cat5_queryRunFeature;
plotData.selectedRunGroupsAnalysisTable.cat1 = cat1_selectedRunGroupsAnalysisTable;
plotData.selectedRunGroupsAnalysisTable.cat2 = cat2_selectedRunGroupsAnalysisTable;
plotData.selectedRunGroupsAnalysisTable.cat3 = cat3_selectedRunGroupsAnalysisTable;
plotData.selectedRunGroupsAnalysisTable.cat4 = cat4_selectedRunGroupsAnalysisTable;
plotData.selectedRunGroupsAnalysisTable.cat5 = cat5_selectedRunGroupsAnalysisTable;

%% Plot
nexttile(masterLayout,6,[1 1]); axis off;

options = generateOptionsStruct(groupByCell);
options = generateOptionsStruct(groupByCell);
options.plotColor.cat1 = colorBlue;
options.plotColor.cat2 = colorBlue*0.9;
options.plotColor.cat3 = colorBlue*0.8;
options.plotColor.cat4 = colorBlue*0.7;
options.plotColor.cat5 = colorBlue*0.6;
options.plotTitle = ['Full-field responses'];
options.legendText.cat1 = ['25 ms'];
options.legendText.cat2 = ['50 ms'];
options.legendText.cat3 = ['200 ms'];
options.legendText.cat4 = ['500 ms'];
options.legendText.cat5 = ['1000 ms'];
options.xLimits = [0,20000];

generatePlotTraces(plotData,options)

%% ** Barplot blue heightPulsePeak different width **

%% Collect data

cat1_queryRunFeature = 'heightPulsePeak';
cat2_queryRunFeature = 'heightPulsePeak';
cat3_queryRunFeature = 'heightPulsePeak';
cat4_queryRunFeature = 'heightPulsePeak';
cat5_queryRunFeature = 'heightPulsePeak';

cat1_selectedCellsIDs = iAnalyzedCell;
cat2_selectedCellsIDs = iAnalyzedCell;
cat3_selectedCellsIDs = iAnalyzedCell;
cat4_selectedCellsIDs = iAnalyzedCell;
cat5_selectedCellsIDs = iAnalyzedCell;

cat1_desiredEpochs = [];
cat2_desiredEpochs = [];
cat3_desiredEpochs = [];
cat4_desiredEpochs = [];
cat5_desiredEpochs = [];

cat1_queryParameters = [];
cat2_queryParameters = [];
cat3_queryParameters = [];
cat4_queryParameters = [];
cat5_queryParameters = [];

cat1_queryParameters.amplitudeBlue = 5;
cat1_queryParameters.pulseWidthBlue = 250;
cat1_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'blue','inhibitory',cat1_queryParameters);

cat2_queryParameters.amplitudeBlue = 5;
cat2_queryParameters.pulseWidthBlue = 500;
cat2_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'blue','inhibitory',cat2_queryParameters);

cat3_queryParameters.amplitudeBlue = 5;
cat3_queryParameters.pulseWidthBlue = 2000;
cat3_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'blue','inhibitory',cat3_queryParameters);

cat4_queryParameters.amplitudeBlue = 5;
cat4_queryParameters.pulseWidthBlue = 5000;
cat4_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'blue','inhibitory',cat4_queryParameters);

cat5_queryParameters.amplitudeBlue = 5;
cat5_queryParameters.pulseWidthBlue = 10000;
cat5_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'blue','inhibitory',cat5_queryParameters);


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

plotData = generatePlotDataStruct;
plotData.queryRunFeature.cat1 = cat1_queryRunFeature;
plotData.queryRunFeature.cat2 = cat2_queryRunFeature;
plotData.queryRunFeature.cat3 = cat3_queryRunFeature;
plotData.queryRunFeature.cat4 = cat4_queryRunFeature;
plotData.queryRunFeature.cat5 = cat5_queryRunFeature;
plotData.selectedRunGroupsAnalysisTable.cat1 = cat1_selectedRunGroupsAnalysisTable;
plotData.selectedRunGroupsAnalysisTable.cat2 = cat2_selectedRunGroupsAnalysisTable;
plotData.selectedRunGroupsAnalysisTable.cat3 = cat3_selectedRunGroupsAnalysisTable;
plotData.selectedRunGroupsAnalysisTable.cat4 = cat4_selectedRunGroupsAnalysisTable;
plotData.selectedRunGroupsAnalysisTable.cat5 = cat5_selectedRunGroupsAnalysisTable;

%% Plot
nexttile(masterLayout,7,[1 1]); axis off;

options = generateOptionsStruct(groupByCell);
options.plotColor.cat1 = colorBlue;
options.plotColor.cat2 = colorBlue*0.9;
options.plotColor.cat3 = colorBlue*0.8;
options.plotColor.cat4 = colorBlue*0.7;
options.plotColor.cat5 = colorBlue*0.6;
options.plotColorOverlay.cat1 = colorBlueOverlay;
options.plotColorOverlay.cat2 = colorBlueOverlay*0.9;
options.plotColorOverlay.cat3 = colorBlueOverlay*0.8;
options.plotColorOverlay.cat4 = colorBlueOverlay*0.7;
options.plotColorOverlay.cat5 = colorBlueOverlay*0.6;
options.plotTitle = ['Response amplitude'];
options.xLabelText.cat1 = ['25 ms'];
options.xLabelText.cat2 = ['50 ms'];
options.xLabelText.cat3 = ['200 ms'];
options.xLabelText.cat4 = ['500 ms'];
options.xLabelText.cat5 = ['1000 ms'];
options.yLabelText = ['Current [pA]'];

generateBarplot(plotData,options)

%% ** Barplot blue areaPulse different widths **

%% Collect data

cat1_queryRunFeature = 'areaPulse';
cat2_queryRunFeature = 'areaPulse';
cat3_queryRunFeature = 'areaPulse';
cat4_queryRunFeature = 'areaPulse';
cat5_queryRunFeature = 'areaPulse';

cat1_selectedCellsIDs = iAnalyzedCell;
cat2_selectedCellsIDs = iAnalyzedCell;
cat3_selectedCellsIDs = iAnalyzedCell;
cat4_selectedCellsIDs = iAnalyzedCell;
cat5_selectedCellsIDs = iAnalyzedCell;

cat1_desiredEpochs = [];
cat2_desiredEpochs = [];
cat3_desiredEpochs = [];
cat4_desiredEpochs = [];
cat5_desiredEpochs = [];

cat1_queryParameters = [];
cat2_queryParameters = [];
cat3_queryParameters = [];
cat4_queryParameters = [];
cat5_queryParameters = [];

cat1_queryParameters.amplitudeBlue = 5;
cat1_queryParameters.pulseWidthBlue = 250;
cat1_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'blue','inhibitory',cat1_queryParameters);

cat2_queryParameters.amplitudeBlue = 5;
cat2_queryParameters.pulseWidthBlue = 500;
cat2_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'blue','inhibitory',cat2_queryParameters);

cat3_queryParameters.amplitudeBlue = 5;
cat3_queryParameters.pulseWidthBlue = 2000;
cat3_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'blue','inhibitory',cat3_queryParameters);

cat4_queryParameters.amplitudeBlue = 5;
cat4_queryParameters.pulseWidthBlue = 5000;
cat4_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'blue','inhibitory',cat4_queryParameters);

cat5_queryParameters.amplitudeBlue = 5;
cat5_queryParameters.pulseWidthBlue = 10000;
cat5_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'blue','inhibitory',cat5_queryParameters);


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

plotData = generatePlotDataStruct;
plotData.queryRunFeature.cat1 = cat1_queryRunFeature;
plotData.queryRunFeature.cat2 = cat2_queryRunFeature;
plotData.queryRunFeature.cat3 = cat3_queryRunFeature;
plotData.queryRunFeature.cat4 = cat4_queryRunFeature;
plotData.queryRunFeature.cat5 = cat5_queryRunFeature;
plotData.selectedRunGroupsAnalysisTable.cat1 = cat1_selectedRunGroupsAnalysisTable;
plotData.selectedRunGroupsAnalysisTable.cat2 = cat2_selectedRunGroupsAnalysisTable;
plotData.selectedRunGroupsAnalysisTable.cat3 = cat3_selectedRunGroupsAnalysisTable;
plotData.selectedRunGroupsAnalysisTable.cat4 = cat4_selectedRunGroupsAnalysisTable;
plotData.selectedRunGroupsAnalysisTable.cat5 = cat5_selectedRunGroupsAnalysisTable;

%% Plot
nexttile(masterLayout,8,[1 1]); axis off;

options = generateOptionsStruct(groupByCell);
options.plotColor.cat1 = colorBlue;
options.plotColor.cat2 = colorBlue*0.9;
options.plotColor.cat3 = colorBlue*0.8;
options.plotColor.cat4 = colorBlue*0.7;
options.plotColor.cat5 = colorBlue*0.6;
options.plotColorOverlay.cat1 = colorBlueOverlay;
options.plotColorOverlay.cat2 = colorBlueOverlay*0.9;
options.plotColorOverlay.cat3 = colorBlueOverlay*0.8;
options.plotColorOverlay.cat4 = colorBlueOverlay*0.7;
options.plotColorOverlay.cat5 = colorBlueOverlay*0.6;
options.plotTitle = ['Response charge'];
options.xLabelText.cat1 = ['25 ms'];
options.xLabelText.cat2 = ['50 ms'];
options.xLabelText.cat3 = ['200 ms'];
options.xLabelText.cat4 = ['500 ms'];
options.xLabelText.cat5 = ['1000 ms'];
options.yLabelText = ['Charge [pC]'];

generateBarplot(plotData,options)

%% Title and saving

plotMainTitle = ['Cell ', num2str(iAnalyzedCell),' - ',experimentName, postTextSuptitle];
sgtitle(plotMainTitle,'FontSize',14,'Interpreter', 'none')

plotName = ['cell', num2str(iAnalyzedCell), '_', 'Plot_' experimentName, postTextSuptitle];
savePlotPath = [experimentDirectory filesep 'mouseAnalysis'];
plotFullPath = fullfile(savePlotPath, plotName);

if savePlots == 1
    exportgraphics(gcf, fullfile(strcat(plotFullPath,'.pdf')), 'ContentType', 'vector');
    saveas(gcf, plotFullPath, 'fig')
end

end