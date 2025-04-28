disp(['*-*-*-* Running: scriptPlotOverallExperimentTwoInputs *-*-*-*'])

overallExperimentDirectory = [matlabDirectory 'GroupedExperiments' filesep overallExperimentName];
experimentAnalysisTablePath = [overallExperimentDirectory filesep 'ExperimentAnalysisTable.mat'];
experimentParametersTablePath = [overallExperimentDirectory filesep 'ExperimentParametersTable.mat'];
experimentCellsConnectionInfoPath = [overallExperimentDirectory filesep 'ExperimentCellsConnectionInfo.mat'];

experimentAnalysisTable = load(experimentAnalysisTablePath);
experimentAnalysisTable = experimentAnalysisTable.experimentAnalysisTable;
experimentParametersTable = load(experimentParametersTablePath);
runGroupsParametersTable = experimentParametersTable.runGroupsParametersTable;
experimentCellsConnectionInfo = load(experimentCellsConnectionInfoPath);
experimentCellsConnectionInfo = experimentCellsConnectionInfo.experimentCellsConnectionInfo;

analyzedCells = unique(cell2mat(experimentAnalysisTable.overallCellName));

initializeFig(0.9,0.9);

masterLayout = tiledlayout(2,8);
masterLayout.TileSpacing = 'compact';
masterLayout.Padding = 'compact';

colorBlue = [85, 161, 254]./255;
colorRed = [255, 50, 58]./255;
colorGrey = [192, 192, 192]./255;
colorPurple = [232 22 224]./255;

colorBlueOverlay = [85, 161, 254]./255 * 0.6;
colorRedOverlay = [255, 50, 58]./255 * 0.6;
colorGreyOverlay = [192, 192, 192]./255 * 0.6;
colorPurpleOverlay = [232 22 224]./255 * 0.6;

groupByCell = 1;
isOverallExperimentPlot = 1;

nonConnectedCellsBlue = cell2mat(experimentCellsConnectionInfo.blue.nonConnectedCells);
connectedCellsBlue = cell2mat(experimentCellsConnectionInfo.blue.connectedCells);
nonConnectedCellsRed = cell2mat(experimentCellsConnectionInfo.red.nonConnectedCells);
connectedCellsRed = cell2mat(experimentCellsConnectionInfo.red.connectedCells);
doubleConnectedCells = cell2mat(experimentCellsConnectionInfo.doubleConnectedCells(~cellfun('isempty', experimentCellsConnectionInfo.doubleConnectedCells)));

experimentCellsConnectionInfo.doubleConnectedCells;

%% ** Plot blue traces **

%% Collect data

cat1_selectedCellsIDs = nonConnectedCellsBlue;
cat2_selectedCellsIDs = connectedCellsBlue;
cat1_desiredEpochs = [];
cat2_desiredEpochs = [];
cat1_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'blue',responseTypeBlue,queryParameters);
cat2_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'blue',responseTypeBlue,queryParameters);

% No need to modify
cat1_selectedMembersTable = extractSelectedMembersFromTable(experimentAnalysisTable, 'overallCellName', cat1_selectedCellsIDs);
cat1_selectedRunGroupsAnalysisTable = cat1_selectedMembersTable(ismember(cat1_selectedMembersTable.runGroup,cat1_selectedRunGroups),:);
cat1_selectedRunGroupsAnalysisTable = extractEpochSubset(cat1_selectedRunGroupsAnalysisTable, cat1_desiredEpochs);

cat2_selectedMembersTable = extractSelectedMembersFromTable(experimentAnalysisTable, 'overallCellName', cat2_selectedCellsIDs);
cat2_selectedRunGroupsAnalysisTable = cat2_selectedMembersTable(ismember(cat2_selectedMembersTable.runGroup,cat2_selectedRunGroups),:);
cat2_selectedRunGroupsAnalysisTable = extractEpochSubset(cat2_selectedRunGroupsAnalysisTable, cat2_desiredEpochs);

plotData = generatePlotDataStruct;
plotData.selectedRunGroupsAnalysisTable.cat1 = cat1_selectedRunGroupsAnalysisTable;
plotData.selectedRunGroupsAnalysisTable.cat2 = cat2_selectedRunGroupsAnalysisTable;

nonConnectedCellsBlueFiltered = unique(cell2mat(cat1_selectedRunGroupsAnalysisTable.overallCellName));
connectedCellsBlueFiltered = unique(cell2mat(cat2_selectedRunGroupsAnalysisTable.overallCellName));

%% Plot
nexttile(masterLayout,1,[1 2]); axis off;

options = generateOptionsStruct(groupByCell,isOverallExperimentPlot);
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
cat1_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'blue',responseTypeBlue,queryParameters);
cat2_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'blue',responseTypeBlue,queryParameters);

% No need to modify
cat1_selectedMembersTable = extractSelectedMembersFromTable(experimentAnalysisTable, 'overallCellName', cat1_selectedCellsIDs);
cat1_selectedRunGroupsAnalysisTable = cat1_selectedMembersTable(ismember(cat1_selectedMembersTable.runGroup,cat1_selectedRunGroups),:);
cat1_selectedRunGroupsAnalysisTable = extractEpochSubset(cat1_selectedRunGroupsAnalysisTable, cat1_desiredEpochs);

cat2_selectedMembersTable = extractSelectedMembersFromTable(experimentAnalysisTable, 'overallCellName', cat2_selectedCellsIDs);
cat2_selectedRunGroupsAnalysisTable = cat2_selectedMembersTable(ismember(cat2_selectedMembersTable.runGroup,cat2_selectedRunGroups),:);
cat2_selectedRunGroupsAnalysisTable = extractEpochSubset(cat2_selectedRunGroupsAnalysisTable, cat2_desiredEpochs);

plotData = generatePlotDataStruct;
plotData.queryRunFeature.cat1 = cat1_queryRunFeature;
plotData.queryRunFeature.cat2 = cat2_queryRunFeature;
plotData.selectedRunGroupsAnalysisTable.cat1 = cat1_selectedRunGroupsAnalysisTable;
plotData.selectedRunGroupsAnalysisTable.cat2 = cat2_selectedRunGroupsAnalysisTable;

%% Plot
nexttile(masterLayout,3,[1 1]); axis off;

options = generateOptionsStruct(groupByCell,isOverallExperimentPlot);
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
cat1_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'blue',responseTypeBlue,queryParameters);
cat2_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'blue',responseTypeBlue,queryParameters);

% No need to modify
cat1_selectedMembersTable = extractSelectedMembersFromTable(experimentAnalysisTable, 'overallCellName', cat1_selectedCellsIDs);
cat1_selectedRunGroupsAnalysisTable = cat1_selectedMembersTable(ismember(cat1_selectedMembersTable.runGroup,cat1_selectedRunGroups),:);
cat1_selectedRunGroupsAnalysisTable = extractEpochSubset(cat1_selectedRunGroupsAnalysisTable, cat1_desiredEpochs);

cat2_selectedMembersTable = extractSelectedMembersFromTable(experimentAnalysisTable, 'overallCellName', cat2_selectedCellsIDs);
cat2_selectedRunGroupsAnalysisTable = cat2_selectedMembersTable(ismember(cat2_selectedMembersTable.runGroup,cat2_selectedRunGroups),:);
cat2_selectedRunGroupsAnalysisTable = extractEpochSubset(cat2_selectedRunGroupsAnalysisTable, cat2_desiredEpochs);

plotData = generatePlotDataStruct;
plotData.queryRunFeature.cat1 = cat1_queryRunFeature;
plotData.queryRunFeature.cat2 = cat2_queryRunFeature;
plotData.selectedRunGroupsAnalysisTable.cat1 = cat1_selectedRunGroupsAnalysisTable;
plotData.selectedRunGroupsAnalysisTable.cat2 = cat2_selectedRunGroupsAnalysisTable;

%% Plot
nexttile(masterLayout,4,[1 1]); axis off;

options = generateOptionsStruct(groupByCell,isOverallExperimentPlot);
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
cat1_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'blue',responseTypeBlue,queryParameters);
cat2_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'blue',responseTypeBlue,queryParameters);

% No need to modify
cat1_selectedMembersTable = extractSelectedMembersFromTable(experimentAnalysisTable, 'overallCellName', cat1_selectedCellsIDs);
cat1_selectedRunGroupsAnalysisTable = cat1_selectedMembersTable(ismember(cat1_selectedMembersTable.runGroup,cat1_selectedRunGroups),:);
cat1_selectedRunGroupsAnalysisTable = extractEpochSubset(cat1_selectedRunGroupsAnalysisTable, cat1_desiredEpochs);

cat2_selectedMembersTable = extractSelectedMembersFromTable(experimentAnalysisTable, 'overallCellName', cat2_selectedCellsIDs);
cat2_selectedRunGroupsAnalysisTable = cat2_selectedMembersTable(ismember(cat2_selectedMembersTable.runGroup,cat2_selectedRunGroups),:);
cat2_selectedRunGroupsAnalysisTable = extractEpochSubset(cat2_selectedRunGroupsAnalysisTable, cat2_desiredEpochs);

plotData = generatePlotDataStruct;
plotData.queryRunFeature.cat1 = cat1_queryRunFeature;
plotData.queryRunFeature.cat2 = cat2_queryRunFeature;
plotData.selectedRunGroupsAnalysisTable.cat1 = cat1_selectedRunGroupsAnalysisTable;
plotData.selectedRunGroupsAnalysisTable.cat2 = cat2_selectedRunGroupsAnalysisTable;

%% Plot
nexttile(masterLayout,5,[1 1]); axis off;

options = generateOptionsStruct(groupByCell,isOverallExperimentPlot);
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
cat1_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'red',responseTypeRed,queryParameters);
cat2_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'red',responseTypeRed,queryParameters);

% No need to modify
cat1_selectedMembersTable = extractSelectedMembersFromTable(experimentAnalysisTable, 'overallCellName', cat1_selectedCellsIDs);
cat1_selectedRunGroupsAnalysisTable = cat1_selectedMembersTable(ismember(cat1_selectedMembersTable.runGroup,cat1_selectedRunGroups),:);
cat1_selectedRunGroupsAnalysisTable = extractEpochSubset(cat1_selectedRunGroupsAnalysisTable, cat1_desiredEpochs);

cat2_selectedMembersTable = extractSelectedMembersFromTable(experimentAnalysisTable, 'overallCellName', cat2_selectedCellsIDs);
cat2_selectedRunGroupsAnalysisTable = cat2_selectedMembersTable(ismember(cat2_selectedMembersTable.runGroup,cat2_selectedRunGroups),:);
cat2_selectedRunGroupsAnalysisTable = extractEpochSubset(cat2_selectedRunGroupsAnalysisTable, cat2_desiredEpochs);

plotData = generatePlotDataStruct;
plotData.selectedRunGroupsAnalysisTable.cat1 = cat1_selectedRunGroupsAnalysisTable;
plotData.selectedRunGroupsAnalysisTable.cat2 = cat2_selectedRunGroupsAnalysisTable;

nonConnectedCellsRedFiltered = unique(cell2mat(cat1_selectedRunGroupsAnalysisTable.overallCellName));
connectedCellsRedFiltered = unique(cell2mat(cat2_selectedRunGroupsAnalysisTable.overallCellName));

%% Plot
nexttile(masterLayout,9,[1 2]); axis off;

options = generateOptionsStruct(groupByCell,isOverallExperimentPlot);
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
cat1_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'red',responseTypeRed,queryParameters);
cat2_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'red',responseTypeRed,queryParameters);

% No need to modify
cat1_selectedMembersTable = extractSelectedMembersFromTable(experimentAnalysisTable, 'overallCellName', cat1_selectedCellsIDs);
cat1_selectedRunGroupsAnalysisTable = cat1_selectedMembersTable(ismember(cat1_selectedMembersTable.runGroup,cat1_selectedRunGroups),:);
cat1_selectedRunGroupsAnalysisTable = extractEpochSubset(cat1_selectedRunGroupsAnalysisTable, cat1_desiredEpochs);

cat2_selectedMembersTable = extractSelectedMembersFromTable(experimentAnalysisTable, 'overallCellName', cat2_selectedCellsIDs);
cat2_selectedRunGroupsAnalysisTable = cat2_selectedMembersTable(ismember(cat2_selectedMembersTable.runGroup,cat2_selectedRunGroups),:);
cat2_selectedRunGroupsAnalysisTable = extractEpochSubset(cat2_selectedRunGroupsAnalysisTable, cat2_desiredEpochs);

plotData = generatePlotDataStruct;
plotData.queryRunFeature.cat1 = cat1_queryRunFeature;
plotData.queryRunFeature.cat2 = cat2_queryRunFeature;
plotData.selectedRunGroupsAnalysisTable.cat1 = cat1_selectedRunGroupsAnalysisTable;
plotData.selectedRunGroupsAnalysisTable.cat2 = cat2_selectedRunGroupsAnalysisTable;

%% Plot
nexttile(masterLayout,11,[1 1]); axis off;

options = generateOptionsStruct(groupByCell,isOverallExperimentPlot);
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
cat1_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'red',responseTypeRed,queryParameters);
cat2_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'red',responseTypeRed,queryParameters);

% No need to modify
cat1_selectedMembersTable = extractSelectedMembersFromTable(experimentAnalysisTable, 'overallCellName', cat1_selectedCellsIDs);
cat1_selectedRunGroupsAnalysisTable = cat1_selectedMembersTable(ismember(cat1_selectedMembersTable.runGroup,cat1_selectedRunGroups),:);
cat1_selectedRunGroupsAnalysisTable = extractEpochSubset(cat1_selectedRunGroupsAnalysisTable, cat1_desiredEpochs);

cat2_selectedMembersTable = extractSelectedMembersFromTable(experimentAnalysisTable, 'overallCellName', cat2_selectedCellsIDs);
cat2_selectedRunGroupsAnalysisTable = cat2_selectedMembersTable(ismember(cat2_selectedMembersTable.runGroup,cat2_selectedRunGroups),:);
cat2_selectedRunGroupsAnalysisTable = extractEpochSubset(cat2_selectedRunGroupsAnalysisTable, cat2_desiredEpochs);

plotData = generatePlotDataStruct;
plotData.queryRunFeature.cat1 = cat1_queryRunFeature;
plotData.queryRunFeature.cat2 = cat2_queryRunFeature;
plotData.selectedRunGroupsAnalysisTable.cat1 = cat1_selectedRunGroupsAnalysisTable;
plotData.selectedRunGroupsAnalysisTable.cat2 = cat2_selectedRunGroupsAnalysisTable;

%% Plot
nexttile(masterLayout,12,[1 1]); axis off;

options = generateOptionsStruct(groupByCell,isOverallExperimentPlot);
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
cat1_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'red',responseTypeRed,queryParameters);
cat2_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'red',responseTypeRed,queryParameters);

% No need to modify
cat1_selectedMembersTable = extractSelectedMembersFromTable(experimentAnalysisTable, 'overallCellName', cat1_selectedCellsIDs);
cat1_selectedRunGroupsAnalysisTable = cat1_selectedMembersTable(ismember(cat1_selectedMembersTable.runGroup,cat1_selectedRunGroups),:);
cat1_selectedRunGroupsAnalysisTable = extractEpochSubset(cat1_selectedRunGroupsAnalysisTable, cat1_desiredEpochs);

cat2_selectedMembersTable = extractSelectedMembersFromTable(experimentAnalysisTable, 'overallCellName', cat2_selectedCellsIDs);
cat2_selectedRunGroupsAnalysisTable = cat2_selectedMembersTable(ismember(cat2_selectedMembersTable.runGroup,cat2_selectedRunGroups),:);
cat2_selectedRunGroupsAnalysisTable = extractEpochSubset(cat2_selectedRunGroupsAnalysisTable, cat2_desiredEpochs);

plotData = generatePlotDataStruct;
plotData.queryRunFeature.cat1 = cat1_queryRunFeature;
plotData.queryRunFeature.cat2 = cat2_queryRunFeature;
plotData.selectedRunGroupsAnalysisTable.cat1 = cat1_selectedRunGroupsAnalysisTable;
plotData.selectedRunGroupsAnalysisTable.cat2 = cat2_selectedRunGroupsAnalysisTable;

%% Plot
nexttile(masterLayout,13,[1 1]); axis off;

options = generateOptionsStruct(groupByCell,isOverallExperimentPlot);
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
cat1_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'blue',responseTypeBlue,queryParameters);
cat2_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'red',responseTypeRed,queryParameters);

% No need to modify
cat1_selectedMembersTable = extractSelectedMembersFromTable(experimentAnalysisTable, 'overallCellName', cat1_selectedCellsIDs);
cat1_selectedRunGroupsAnalysisTable = cat1_selectedMembersTable(ismember(cat1_selectedMembersTable.runGroup,cat1_selectedRunGroups),:);
cat1_selectedRunGroupsAnalysisTable = extractEpochSubset(cat1_selectedRunGroupsAnalysisTable, cat1_desiredEpochs);

cat2_selectedMembersTable = extractSelectedMembersFromTable(experimentAnalysisTable, 'overallCellName', cat2_selectedCellsIDs);
cat2_selectedRunGroupsAnalysisTable = cat2_selectedMembersTable(ismember(cat2_selectedMembersTable.runGroup,cat2_selectedRunGroups),:);
cat2_selectedRunGroupsAnalysisTable = extractEpochSubset(cat2_selectedRunGroupsAnalysisTable, cat2_desiredEpochs);

plotData = generatePlotDataStruct;
plotData.selectedRunGroupsAnalysisTable.cat1 = cat1_selectedRunGroupsAnalysisTable;
plotData.selectedRunGroupsAnalysisTable.cat2 = cat2_selectedRunGroupsAnalysisTable;
    
%% Plot
nexttile(masterLayout,14,[1 2]); axis off;

options = generateOptionsStruct(groupByCell,isOverallExperimentPlot);
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
cat1_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'blue',responseTypeBlue,queryParameters);
cat2_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'red',responseTypeRed,queryParameters);

% No need to modify
cat1_selectedMembersTable = extractSelectedMembersFromTable(experimentAnalysisTable, 'overallCellName', cat1_selectedCellsIDs);
cat1_selectedRunGroupsAnalysisTable = cat1_selectedMembersTable(ismember(cat1_selectedMembersTable.runGroup,cat1_selectedRunGroups),:);
cat1_selectedRunGroupsAnalysisTable = extractEpochSubset(cat1_selectedRunGroupsAnalysisTable, cat1_desiredEpochs);

cat2_selectedMembersTable = extractSelectedMembersFromTable(experimentAnalysisTable, 'overallCellName', cat2_selectedCellsIDs);
cat2_selectedRunGroupsAnalysisTable = cat2_selectedMembersTable(ismember(cat2_selectedMembersTable.runGroup,cat2_selectedRunGroups),:);
cat2_selectedRunGroupsAnalysisTable = extractEpochSubset(cat2_selectedRunGroupsAnalysisTable, cat2_desiredEpochs);

plotData = generatePlotDataStruct;
plotData.queryRunFeature.cat1 = cat1_queryRunFeature;
plotData.queryRunFeature.cat2 = cat2_queryRunFeature;
plotData.selectedRunGroupsAnalysisTable.cat1 = cat1_selectedRunGroupsAnalysisTable;
plotData.selectedRunGroupsAnalysisTable.cat2 = cat2_selectedRunGroupsAnalysisTable;

%% Plot
nexttile(masterLayout,16,[1 1]); axis off;

options = generateOptionsStruct(groupByCell,isOverallExperimentPlot);
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
cat1_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'blue',responseTypeBlue,queryParameters);
cat2_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'red',responseTypeRed,queryParameters);

% No need to modify
cat1_selectedMembersTable = extractSelectedMembersFromTable(experimentAnalysisTable, 'overallCellName', cat1_selectedCellsIDs);
cat1_selectedRunGroupsAnalysisTable = cat1_selectedMembersTable(ismember(cat1_selectedMembersTable.runGroup,cat1_selectedRunGroups),:);
cat1_selectedRunGroupsAnalysisTable = extractEpochSubset(cat1_selectedRunGroupsAnalysisTable, cat1_desiredEpochs);

cat2_selectedMembersTable = extractSelectedMembersFromTable(experimentAnalysisTable, 'overallCellName', cat2_selectedCellsIDs);
cat2_selectedRunGroupsAnalysisTable = cat2_selectedMembersTable(ismember(cat2_selectedMembersTable.runGroup,cat2_selectedRunGroups),:);
cat2_selectedRunGroupsAnalysisTable = extractEpochSubset(cat2_selectedRunGroupsAnalysisTable, cat2_desiredEpochs);

plotData = generatePlotDataStruct;
plotData.queryRunFeature.cat1 = cat1_queryRunFeature;
plotData.queryRunFeature.cat2 = cat2_queryRunFeature;
plotData.selectedRunGroupsAnalysisTable.cat1 = cat1_selectedRunGroupsAnalysisTable;
plotData.selectedRunGroupsAnalysisTable.cat2 = cat2_selectedRunGroupsAnalysisTable;

%% Plot
nexttile(masterLayout,8,[1 1]); axis off;

options = generateOptionsStruct(groupByCell,isOverallExperimentPlot);
options.plotColor.cat1 = colorBlue;
options.plotColor.cat2 = colorRed;
options.plotColorOverlay.cat1 = colorGreyOverlay;
options.plotColorOverlay.cat2 = colorRedOverlay;
options.plotTitle = ['Double-connected cells' newline 'Response delay'];
options.xLabelText.cat1 = ['Blue'];
options.xLabelText.cat2 = ['Red'];
options.yLabelText = ['Time [ms]'];

generateBarplot(plotData,options)

%% Pie chart
nexttile(masterLayout,6,[1 2]); axis off;

nConnectedCellsBlue = numel(connectedCellsBlueFiltered);
nNonConnectedCellsBlue = numel(nonConnectedCellsBlueFiltered);
nConnectedCellsRed = numel(connectedCellsRedFiltered);
nNonConnectedCellsRed = numel(nonConnectedCellsRedFiltered);
nDoubleConnectedCells = numel(intersect(connectedCellsBlueFiltered,connectedCellsRedFiltered));

nConnectedOnlyBlue = nConnectedCellsBlue-nDoubleConnectedCells;
nConnectedOnlyRed = nConnectedCellsRed-nDoubleConnectedCells;

nTotalCellsBlue = nConnectedCellsBlue + nNonConnectedCellsBlue;
nTotalCellsRed = nConnectedCellsRed + nNonConnectedCellsRed;

connectedCellsOnlyBlueLabel = ['Blue-only connected cells (' num2str(nConnectedOnlyBlue) '/' num2str(nTotalCellsBlue) ')'];
nonConnectedCellsBlueLabel = ['Blue non-connected cells (' num2str(nNonConnectedCellsBlue) '/' num2str(nTotalCellsBlue) ')'];
connectedCellsOnlyRedLabel = ['Red-only connected cells (' num2str(nConnectedOnlyRed) '/' num2str(nTotalCellsRed) ')'];
nonConnectedCellsRedLabel = ['Red non-connected cells (' num2str(nNonConnectedCellsRed) '/' num2str(nTotalCellsRed) ')'];
doubleConnectedCellsLabel = ['Double-connected cells (' num2str(nDoubleConnectedCells) '/' num2str(min(nTotalCellsBlue,nTotalCellsRed)) ')'];

chartCategories = {nonConnectedCellsBlueLabel,connectedCellsOnlyBlueLabel,doubleConnectedCellsLabel,connectedCellsOnlyRedLabel,nonConnectedCellsRedLabel};
chartNumbers = [nNonConnectedCellsBlue nConnectedOnlyBlue nDoubleConnectedCells nConnectedOnlyRed nNonConnectedCellsRed];
pieChart = pie(chartNumbers,chartCategories);

legend(chartCategories,'Location','northoutside','FontSize',12)

weightGrey = 0.75;
weightColor = 0.25;

pNonConnectedBlue = pieChart(1);
pNonConnectedBlue.FaceColor = (weightGrey * colorGrey) + (weightColor * colorBlue);
tNonConnectedBlue = pieChart(2);
tNonConnectedBlue.FontSize = 10;
pConnectedBlue = pieChart(3);
pConnectedBlue.FaceColor = colorBlue;
tConnectedBlue = pieChart(4);
tConnectedBlue.FontSize = 10;

pDoubleConnectedBlue = pieChart(5);
pDoubleConnectedBlue.FaceColor = colorPurple;
tDoubleConnectedBlue = pieChart(6);
tDoubleConnectedBlue.FontSize = 10;

pConnectedRed = pieChart(7);
pConnectedRed.FaceColor = colorRed;
tConnectedRed = pieChart(8);
tConnectedRed.FontSize = 10;
pNonConnectedRed = pieChart(9);
pNonConnectedRed.FaceColor = (weightGrey * colorGrey) + (weightColor * colorRed);
tNonConnectedRed = pieChart(10);
tNonConnectedRed.FontSize = 10;

delete(findobj(pieChart, 'Type', 'text'));

end
%% Title and saving

plotMainTitle = [overallExperimentName, postTextSuptitle];
sgtitle(plotMainTitle,'FontSize',14,'Interpreter', 'none')

plotName = ['summaryPlot_' plotMainTitle];
savePlotPath = overallExperimentDirectory;
plotFullPath = fullfile(savePlotPath, plotName);

if savePlots == 1
    exportgraphics(gcf, fullfile(strcat(plotFullPath,'.pdf')), 'ContentType', 'vector');
    saveas(gcf, plotFullPath, 'fig')
end
