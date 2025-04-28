disp(['*-*-*-* Running: scriptPlotAnalysisSummary *-*-*-*'])

mouseAnalysisTablePath = [experimentDirectory filesep 'mouseAnalysis' filesep 'MouseAnalysisTable.mat'];
mouseParametersTablePath = [experimentDirectory filesep 'mouseAnalysis' filesep 'MouseParametersTable.mat'];
cellsConnectionInfoPath = [experimentDirectory filesep 'mouseAnalysis' filesep 'CellsConnectionInfo.mat'];

mouseAnalysisTable = load(mouseAnalysisTablePath);
mouseAnalysisTable = mouseAnalysisTable.mouseAnalysisTable;
mouseParametersTable = load(mouseParametersTablePath);
runGroupsParametersTable = mouseParametersTable.runGroupsParametersTable;
cellsConnectionInfo = load(cellsConnectionInfoPath);
cellsConnectionInfo = cellsConnectionInfo.cellsConnectionInfo;

postTextSuptitle = '';
analyzedCells = unique(cell2mat(mouseAnalysisTable.cellName));

initializeFig(0.9,0.9);

masterLayout = tiledlayout(2,7);
masterLayout.TileSpacing = 'compact';
masterLayout.Padding = 'compact';

colorBlue = [85, 161, 254]./255;
colorRed = [255, 50, 58]./255;
colorGrey = [192, 192, 192]./255;

colorBlueOverlay = [85, 161, 254]./255 * 0.6;
colorRedOverlay = [255, 50, 58]./255 * 0.6;
colorGreyOverlay = [192, 192, 192]./255 * 0.6;

groupByCell = 1;

nonConnectedCellsBlue = cell2mat(cellsConnectionInfo.blue.nonConnectedCells);
connectedCellsBlue = cell2mat(cellsConnectionInfo.blue.connectedCells);
nonConnectedCellsRed = cell2mat(cellsConnectionInfo.red.nonConnectedCells);
connectedCellsRed = cell2mat(cellsConnectionInfo.red.connectedCells);
doubleConnectedCells = cellsConnectionInfo.doubleConnectedCells;

%% ** Plot blue traces **

%% Collect data

cat1_selectedCellsIDs = nonConnectedCellsBlue;
cat2_selectedCellsIDs = connectedCellsBlue;
cat1_desiredEpochs = [];
cat2_desiredEpochs = [];
cat1_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'blue',responseTypeBlue);
cat2_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'blue',responseTypeBlue);

% No need to modify
cat1_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat1_selectedCellsIDs);
cat1_selectedRunGroupsAnalysisTable = cat1_selectedMembersTable(ismember(cat1_selectedMembersTable.runGroup,cat1_selectedRunGroups),:);
cat1_selectedRunGroupsAnalysisTable = extractEpochSubset(cat1_selectedRunGroupsAnalysisTable, cat1_desiredEpochs);

cat2_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat2_selectedCellsIDs);
cat2_selectedRunGroupsAnalysisTable = cat2_selectedMembersTable(ismember(cat2_selectedMembersTable.runGroup,cat2_selectedRunGroups),:);
cat2_selectedRunGroupsAnalysisTable = extractEpochSubset(cat2_selectedRunGroupsAnalysisTable, cat2_desiredEpochs);

plotData = generatePlotDataStruct;
plotData.selectedRunGroupsAnalysisTable.cat1 = cat1_selectedRunGroupsAnalysisTable;
plotData.selectedRunGroupsAnalysisTable.cat2 = cat2_selectedRunGroupsAnalysisTable;

%% Plot
nexttile(masterLayout,1,[1 2]); axis off;

options = generateOptionsStruct(groupByCell);
options.plotColor.cat1 = colorGrey;
options.plotColor.cat2 = colorBlue;
options.plotTitle = ['Full-field responses'];
options.legendText.cat1 = ['Non-connected'];
options.legendText.cat2 = ['Connected'];

generatePlotTraces(plotData,options)

%% ** Barplot blue heightControlPeak vs heightPulsePeak **

%% Collect data

cat1_queryRunFeature = 'heightControlPeak';
cat2_queryRunFeature = 'heightPulsePeak';
cat1_selectedCellsIDs = connectedCellsBlue;
cat2_selectedCellsIDs = connectedCellsBlue;
cat1_desiredEpochs = [];
cat2_desiredEpochs = [];
cat1_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'blue',responseTypeBlue);
cat2_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'blue',responseTypeBlue);

% No need to modify
cat1_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat1_selectedCellsIDs);
cat1_selectedRunGroupsAnalysisTable = cat1_selectedMembersTable(ismember(cat1_selectedMembersTable.runGroup,cat1_selectedRunGroups),:);
cat1_selectedRunGroupsAnalysisTable = extractEpochSubset(cat1_selectedRunGroupsAnalysisTable, cat1_desiredEpochs);

cat2_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat2_selectedCellsIDs);
cat2_selectedRunGroupsAnalysisTable = cat2_selectedMembersTable(ismember(cat2_selectedMembersTable.runGroup,cat2_selectedRunGroups),:);
cat2_selectedRunGroupsAnalysisTable = extractEpochSubset(cat2_selectedRunGroupsAnalysisTable, cat2_desiredEpochs);

plotData = generatePlotDataStruct;
plotData.queryRunFeature.cat1 = cat1_queryRunFeature;
plotData.queryRunFeature.cat2 = cat2_queryRunFeature;
plotData.selectedRunGroupsAnalysisTable.cat1 = cat1_selectedRunGroupsAnalysisTable;
plotData.selectedRunGroupsAnalysisTable.cat2 = cat2_selectedRunGroupsAnalysisTable;

%% Plot
nexttile(masterLayout,3,[1 1]); axis off;

options = generateOptionsStruct(groupByCell);
options.plotColor.cat1 = colorGrey;
options.plotColor.cat2 = colorBlue;
options.plotColorOverlay.cat1 = colorGreyOverlay;
options.plotColorOverlay.cat2 = colorBlueOverlay;
options.plotTitle = ['Response amplitude' newline 'Connected cells'];
options.xLabelText.cat1 = ['Control'];
options.xLabelText.cat2 = ['Light'];
options.yLabelText = ['Current [pA]'];

generateBarplot(plotData,options)

%% ** Barplot blue areaControl vs areaPulse **

%% Collect data

cat1_queryRunFeature = 'areaControl';
cat2_queryRunFeature = 'areaPulse';
cat1_selectedCellsIDs = connectedCellsBlue;
cat2_selectedCellsIDs = connectedCellsBlue;
cat1_desiredEpochs = [];
cat2_desiredEpochs = [];
cat1_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'blue',responseTypeBlue);
cat2_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'blue',responseTypeBlue);

% No need to modify
cat1_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat1_selectedCellsIDs);
cat1_selectedRunGroupsAnalysisTable = cat1_selectedMembersTable(ismember(cat1_selectedMembersTable.runGroup,cat1_selectedRunGroups),:);
cat1_selectedRunGroupsAnalysisTable = extractEpochSubset(cat1_selectedRunGroupsAnalysisTable, cat1_desiredEpochs);

cat2_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat2_selectedCellsIDs);
cat2_selectedRunGroupsAnalysisTable = cat2_selectedMembersTable(ismember(cat2_selectedMembersTable.runGroup,cat2_selectedRunGroups),:);
cat2_selectedRunGroupsAnalysisTable = extractEpochSubset(cat2_selectedRunGroupsAnalysisTable, cat2_desiredEpochs);

plotData = generatePlotDataStruct;
plotData.queryRunFeature.cat1 = cat1_queryRunFeature;
plotData.queryRunFeature.cat2 = cat2_queryRunFeature;
plotData.selectedRunGroupsAnalysisTable.cat1 = cat1_selectedRunGroupsAnalysisTable;
plotData.selectedRunGroupsAnalysisTable.cat2 = cat2_selectedRunGroupsAnalysisTable;

%% Plot
nexttile(masterLayout,4,[1 1]); axis off;

options = generateOptionsStruct(groupByCell);
options.plotColor.cat1 = colorGrey;
options.plotColor.cat2 = colorBlue;
options.plotColorOverlay.cat1 = colorGreyOverlay;
options.plotColorOverlay.cat2 = colorBlueOverlay;
options.plotTitle = ['Response charge' newline 'Connected cells'];
options.xLabelText.cat1 = ['Control'];
options.xLabelText.cat2 = ['Light'];
options.yLabelText = ['Charge [pC]'];

generateBarplot(plotData,options)

%% ** Barplot blue non-connected areaControl vs areaPulse **

%% Collect data

cat1_queryRunFeature = 'areaControl';
cat2_queryRunFeature = 'areaPulse';
cat1_selectedCellsIDs = nonConnectedCellsBlue;
cat2_selectedCellsIDs = nonConnectedCellsBlue;
cat1_desiredEpochs = [];
cat2_desiredEpochs = [];
cat1_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'blue',responseTypeBlue);
cat2_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'blue',responseTypeBlue);

% No need to modify
cat1_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat1_selectedCellsIDs);
cat1_selectedRunGroupsAnalysisTable = cat1_selectedMembersTable(ismember(cat1_selectedMembersTable.runGroup,cat1_selectedRunGroups),:);
cat1_selectedRunGroupsAnalysisTable = extractEpochSubset(cat1_selectedRunGroupsAnalysisTable, cat1_desiredEpochs);

cat2_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat2_selectedCellsIDs);
cat2_selectedRunGroupsAnalysisTable = cat2_selectedMembersTable(ismember(cat2_selectedMembersTable.runGroup,cat2_selectedRunGroups),:);
cat2_selectedRunGroupsAnalysisTable = extractEpochSubset(cat2_selectedRunGroupsAnalysisTable, cat2_desiredEpochs);

plotData = generatePlotDataStruct;
plotData.queryRunFeature.cat1 = cat1_queryRunFeature;
plotData.queryRunFeature.cat2 = cat2_queryRunFeature;
plotData.selectedRunGroupsAnalysisTable.cat1 = cat1_selectedRunGroupsAnalysisTable;
plotData.selectedRunGroupsAnalysisTable.cat2 = cat2_selectedRunGroupsAnalysisTable;

%% Plot
nexttile(masterLayout,5,[1 1]); axis off;

options = generateOptionsStruct(groupByCell);
options.plotColor.cat1 = colorGrey;
options.plotColor.cat2 = colorBlue;
options.plotColorOverlay.cat1 = colorGreyOverlay;
options.plotColorOverlay.cat2 = colorBlueOverlay;
options.plotTitle = ['Response charge' newline 'Non-connected cells'];
options.xLabelText.cat1 = ['Control'];
options.xLabelText.cat2 = ['Light'];
options.yLabelText = ['Charge [pC]'];

generateBarplot(plotData,options)

%% ** Plot red traces **

%% Collect data

cat1_selectedCellsIDs = nonConnectedCellsRed;
cat2_selectedCellsIDs = connectedCellsRed;
cat1_desiredEpochs = [];
cat2_desiredEpochs = [];
cat1_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'red',responseTypeRed);
cat2_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'red',responseTypeRed);

% No need to modify
cat1_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat1_selectedCellsIDs);
cat1_selectedRunGroupsAnalysisTable = cat1_selectedMembersTable(ismember(cat1_selectedMembersTable.runGroup,cat1_selectedRunGroups),:);
cat1_selectedRunGroupsAnalysisTable = extractEpochSubset(cat1_selectedRunGroupsAnalysisTable, cat1_desiredEpochs);

cat2_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat2_selectedCellsIDs);
cat2_selectedRunGroupsAnalysisTable = cat2_selectedMembersTable(ismember(cat2_selectedMembersTable.runGroup,cat2_selectedRunGroups),:);
cat2_selectedRunGroupsAnalysisTable = extractEpochSubset(cat2_selectedRunGroupsAnalysisTable, cat2_desiredEpochs);

plotData = generatePlotDataStruct;
plotData.selectedRunGroupsAnalysisTable.cat1 = cat1_selectedRunGroupsAnalysisTable;
plotData.selectedRunGroupsAnalysisTable.cat2 = cat2_selectedRunGroupsAnalysisTable;

%% Plot
nexttile(masterLayout,8,[1 2]); axis off;

options = generateOptionsStruct(groupByCell);
options.plotColor.cat1 = colorGrey;
options.plotColor.cat2 = colorRed;
options.plotTitle = ['Full-field responses'];
options.legendText.cat1 = ['Non-connected'];
options.legendText.cat2 = ['Connected'];

generatePlotTraces(plotData,options)

%% ** Barplot red heightControlPeak vs heightPulsePeak **

%% Collect data

cat1_queryRunFeature = 'heightControlPeak';
cat2_queryRunFeature = 'heightPulsePeak';
cat1_selectedCellsIDs = connectedCellsRed;
cat2_selectedCellsIDs = connectedCellsRed;
cat1_desiredEpochs = [];
cat2_desiredEpochs = [];
cat1_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'red',responseTypeRed);
cat2_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'red',responseTypeRed);

% No need to modify
cat1_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat1_selectedCellsIDs);
cat1_selectedRunGroupsAnalysisTable = cat1_selectedMembersTable(ismember(cat1_selectedMembersTable.runGroup,cat1_selectedRunGroups),:);
cat1_selectedRunGroupsAnalysisTable = extractEpochSubset(cat1_selectedRunGroupsAnalysisTable, cat1_desiredEpochs);

cat2_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat2_selectedCellsIDs);
cat2_selectedRunGroupsAnalysisTable = cat2_selectedMembersTable(ismember(cat2_selectedMembersTable.runGroup,cat2_selectedRunGroups),:);
cat2_selectedRunGroupsAnalysisTable = extractEpochSubset(cat2_selectedRunGroupsAnalysisTable, cat2_desiredEpochs);

plotData = generatePlotDataStruct;
plotData.queryRunFeature.cat1 = cat1_queryRunFeature;
plotData.queryRunFeature.cat2 = cat2_queryRunFeature;
plotData.selectedRunGroupsAnalysisTable.cat1 = cat1_selectedRunGroupsAnalysisTable;
plotData.selectedRunGroupsAnalysisTable.cat2 = cat2_selectedRunGroupsAnalysisTable;

%% Plot
nexttile(masterLayout,10,[1 1]); axis off;

options = generateOptionsStruct(groupByCell);
options.plotColor.cat1 = colorGrey;
options.plotColor.cat2 = colorRed;
options.plotColorOverlay.cat1 = colorGreyOverlay;
options.plotColorOverlay.cat2 = colorRedOverlay;
options.plotTitle = ['Response amplitude' newline 'Connected cells'];
options.xLabelText.cat1 = ['Control'];
options.xLabelText.cat2 = ['Light'];
options.yLabelText = ['Current [pA]'];

generateBarplot(plotData,options)

%% ** Barplot red areaControl vs areaPulse **

%% Collect data

cat1_queryRunFeature = 'areaControl';
cat2_queryRunFeature = 'areaPulse';
cat1_selectedCellsIDs = connectedCellsRed;
cat2_selectedCellsIDs = connectedCellsRed;
cat1_desiredEpochs = [];
cat2_desiredEpochs = [];
cat1_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'red',responseTypeRed);
cat2_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'red',responseTypeRed);

% No need to modify
cat1_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat1_selectedCellsIDs);
cat1_selectedRunGroupsAnalysisTable = cat1_selectedMembersTable(ismember(cat1_selectedMembersTable.runGroup,cat1_selectedRunGroups),:);
cat1_selectedRunGroupsAnalysisTable = extractEpochSubset(cat1_selectedRunGroupsAnalysisTable, cat1_desiredEpochs);

cat2_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat2_selectedCellsIDs);
cat2_selectedRunGroupsAnalysisTable = cat2_selectedMembersTable(ismember(cat2_selectedMembersTable.runGroup,cat2_selectedRunGroups),:);
cat2_selectedRunGroupsAnalysisTable = extractEpochSubset(cat2_selectedRunGroupsAnalysisTable, cat2_desiredEpochs);

plotData = generatePlotDataStruct;
plotData.queryRunFeature.cat1 = cat1_queryRunFeature;
plotData.queryRunFeature.cat2 = cat2_queryRunFeature;
plotData.selectedRunGroupsAnalysisTable.cat1 = cat1_selectedRunGroupsAnalysisTable;
plotData.selectedRunGroupsAnalysisTable.cat2 = cat2_selectedRunGroupsAnalysisTable;

%% Plot
nexttile(masterLayout,11,[1 1]); axis off;

options = generateOptionsStruct(groupByCell);
options.plotColor.cat1 = colorGrey;
options.plotColor.cat2 = colorRed;
options.plotColorOverlay.cat1 = colorGreyOverlay;
options.plotColorOverlay.cat2 = colorRedOverlay;
options.plotTitle = ['Response charge' newline 'Connected cells'];
options.xLabelText.cat1 = ['Control'];
options.xLabelText.cat2 = ['Light'];
options.yLabelText = ['Charge [pC]'];

generateBarplot(plotData,options)

%% ** Barplot red non-connected areaControl vs areaPulse **

%% Collect data

cat1_queryRunFeature = 'areaControl';
cat2_queryRunFeature = 'areaPulse';
cat1_selectedCellsIDs = nonConnectedCellsRed;
cat2_selectedCellsIDs = nonConnectedCellsRed;
cat1_desiredEpochs = [];
cat2_desiredEpochs = [];
cat1_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'red',responseTypeRed);
cat2_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'red',responseTypeRed);

% No need to modify
cat1_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat1_selectedCellsIDs);
cat1_selectedRunGroupsAnalysisTable = cat1_selectedMembersTable(ismember(cat1_selectedMembersTable.runGroup,cat1_selectedRunGroups),:);
cat1_selectedRunGroupsAnalysisTable = extractEpochSubset(cat1_selectedRunGroupsAnalysisTable, cat1_desiredEpochs);

cat2_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat2_selectedCellsIDs);
cat2_selectedRunGroupsAnalysisTable = cat2_selectedMembersTable(ismember(cat2_selectedMembersTable.runGroup,cat2_selectedRunGroups),:);
cat2_selectedRunGroupsAnalysisTable = extractEpochSubset(cat2_selectedRunGroupsAnalysisTable, cat2_desiredEpochs);

plotData = generatePlotDataStruct;
plotData.queryRunFeature.cat1 = cat1_queryRunFeature;
plotData.queryRunFeature.cat2 = cat2_queryRunFeature;
plotData.selectedRunGroupsAnalysisTable.cat1 = cat1_selectedRunGroupsAnalysisTable;
plotData.selectedRunGroupsAnalysisTable.cat2 = cat2_selectedRunGroupsAnalysisTable;

%% Plot
nexttile(masterLayout,12,[1 1]); axis off;

options = generateOptionsStruct(groupByCell);
options.plotColor.cat1 = colorGrey;
options.plotColor.cat2 = colorRed;
options.plotColorOverlay.cat1 = colorGreyOverlay;
options.plotColorOverlay.cat2 = colorRedOverlay;
options.plotTitle = ['Response charge' newline 'Non-connected cells'];
options.xLabelText.cat1 = ['Control'];
options.xLabelText.cat2 = ['Light'];
options.yLabelText = ['Charge [pC]'];

generateBarplot(plotData,options)

%% Double-connected cells
if isempty(doubleConnectedCells) == 0
    
%% Full-field responses

%% Collect data

cat1_selectedCellsIDs = doubleConnectedCells;
cat2_selectedCellsIDs = doubleConnectedCells;
cat1_desiredEpochs = [];
cat2_desiredEpochs = [];
cat1_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'blue',responseTypeBlue);
cat2_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'red',responseTypeRed);

% No need to modify
cat1_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat1_selectedCellsIDs);
cat1_selectedRunGroupsAnalysisTable = cat1_selectedMembersTable(ismember(cat1_selectedMembersTable.runGroup,cat1_selectedRunGroups),:);
cat1_selectedRunGroupsAnalysisTable = extractEpochSubset(cat1_selectedRunGroupsAnalysisTable, cat1_desiredEpochs);

cat2_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat2_selectedCellsIDs);
cat2_selectedRunGroupsAnalysisTable = cat2_selectedMembersTable(ismember(cat2_selectedMembersTable.runGroup,cat2_selectedRunGroups),:);
cat2_selectedRunGroupsAnalysisTable = extractEpochSubset(cat2_selectedRunGroupsAnalysisTable, cat2_desiredEpochs);

plotData = generatePlotDataStruct;
plotData.selectedRunGroupsAnalysisTable.cat1 = cat1_selectedRunGroupsAnalysisTable;
plotData.selectedRunGroupsAnalysisTable.cat2 = cat2_selectedRunGroupsAnalysisTable;
    
%% Plot
nexttile(masterLayout,6,[1 2]); axis off;

options = generateOptionsStruct(groupByCell);
options.plotColor.cat1 = colorBlue;
options.plotColor.cat2 = colorRed;
options.plotTitle = ['Double-connected cells' newline 'Full-field responses'];
options.legendText.cat1 = ['Blue'];
options.legendText.cat2 = ['Red'];

generatePlotTraces(plotData,options)

%% Plot 1 double-connected

%% Collect data

cat1_queryRunFeature = 'areaPulse';
cat2_queryRunFeature = 'areaPulse';
cat1_selectedCellsIDs = doubleConnectedCells;
cat2_selectedCellsIDs = doubleConnectedCells;
cat1_desiredEpochs = [];
cat2_desiredEpochs = [];
cat1_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'blue',responseTypeBlue);
cat2_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'red',responseTypeRed);

% No need to modify
cat1_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat1_selectedCellsIDs);
cat1_selectedRunGroupsAnalysisTable = cat1_selectedMembersTable(ismember(cat1_selectedMembersTable.runGroup,cat1_selectedRunGroups),:);
cat1_selectedRunGroupsAnalysisTable = extractEpochSubset(cat1_selectedRunGroupsAnalysisTable, cat1_desiredEpochs);

cat2_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat2_selectedCellsIDs);
cat2_selectedRunGroupsAnalysisTable = cat2_selectedMembersTable(ismember(cat2_selectedMembersTable.runGroup,cat2_selectedRunGroups),:);
cat2_selectedRunGroupsAnalysisTable = extractEpochSubset(cat2_selectedRunGroupsAnalysisTable, cat2_desiredEpochs);

plotData = generatePlotDataStruct;
plotData.queryRunFeature.cat1 = cat1_queryRunFeature;
plotData.queryRunFeature.cat2 = cat2_queryRunFeature;
plotData.selectedRunGroupsAnalysisTable.cat1 = cat1_selectedRunGroupsAnalysisTable;
plotData.selectedRunGroupsAnalysisTable.cat2 = cat2_selectedRunGroupsAnalysisTable;

%% Plot
nexttile(masterLayout,13,[1 1]); axis off;

options = generateOptionsStruct(groupByCell);
options.plotColor.cat1 = colorBlue;
options.plotColor.cat2 = colorRed;
options.plotColorOverlay.cat1 = colorGreyOverlay;
options.plotColorOverlay.cat2 = colorRedOverlay;
options.plotTitle = ['Double-connected cells' newline 'Response charge'];
options.xLabelText.cat1 = ['Blue'];
options.xLabelText.cat2 = ['Red'];
options.yLabelText = ['Charge [pC]'];

generateBarplot(plotData,options)

%% Plot 2 double-connected

%% Collect data

cat1_queryRunFeature = 'timePulsePeak';
cat2_queryRunFeature = 'timePulsePeak';
cat1_selectedCellsIDs = doubleConnectedCells;
cat2_selectedCellsIDs = doubleConnectedCells;
cat1_desiredEpochs = [];
cat2_desiredEpochs = [];
cat1_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'blue',responseTypeBlue);
cat2_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'red',responseTypeRed);

% No need to modify
cat1_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat1_selectedCellsIDs);
cat1_selectedRunGroupsAnalysisTable = cat1_selectedMembersTable(ismember(cat1_selectedMembersTable.runGroup,cat1_selectedRunGroups),:);
cat1_selectedRunGroupsAnalysisTable = extractEpochSubset(cat1_selectedRunGroupsAnalysisTable, cat1_desiredEpochs);

cat2_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat2_selectedCellsIDs);
cat2_selectedRunGroupsAnalysisTable = cat2_selectedMembersTable(ismember(cat2_selectedMembersTable.runGroup,cat2_selectedRunGroups),:);
cat2_selectedRunGroupsAnalysisTable = extractEpochSubset(cat2_selectedRunGroupsAnalysisTable, cat2_desiredEpochs);

plotData = generatePlotDataStruct;
plotData.queryRunFeature.cat1 = cat1_queryRunFeature;
plotData.queryRunFeature.cat2 = cat2_queryRunFeature;
plotData.selectedRunGroupsAnalysisTable.cat1 = cat1_selectedRunGroupsAnalysisTable;
plotData.selectedRunGroupsAnalysisTable.cat2 = cat2_selectedRunGroupsAnalysisTable;

%% Plot
nexttile(masterLayout,14,[1 1]); axis off;

options = generateOptionsStruct(groupByCell);
options.plotColor.cat1 = colorBlue;
options.plotColor.cat2 = colorRed;
options.plotColorOverlay.cat1 = colorGreyOverlay;
options.plotColorOverlay.cat2 = colorRedOverlay;
options.plotTitle = ['Double-connected cells' newline 'Response delay'];
options.xLabelText.cat1 = ['Blue'];
options.xLabelText.cat2 = ['Red'];
options.yLabelText = ['Time [ms]'];

generateBarplot(plotData,options)

end
%% Title and saving

plotMainTitle = [experimentName, postTextSuptitle];
sgtitle(plotMainTitle,'FontSize',14,'Interpreter', 'none')

plotName = ['summaryPlot_' plotMainTitle];
savePlotPath = [experimentDirectory filesep 'mouseAnalysis'];
plotFullPath = fullfile(savePlotPath, plotName);

if savePlots == 1
    exportgraphics(gcf, fullfile(strcat(plotFullPath,'.pdf')), 'ContentType', 'vector');
    saveas(gcf, plotFullPath, 'fig')
end
