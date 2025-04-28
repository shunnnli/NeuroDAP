disp(['*-*-*-* Running: scriptPlotAnalysisCell *-*-*-*'])

mouseAnalysisTablePath = [experimentDirectory filesep 'mouseAnalysis' filesep 'MouseAnalysisTable.mat'];
mouseParametersTablePath = [experimentDirectory filesep 'mouseAnalysis' filesep 'MouseParametersTable.mat'];

mouseAnalysisTable = load(mouseAnalysisTablePath);
mouseAnalysisTable = mouseAnalysisTable.mouseAnalysisTable;
mouseParametersTable = load(mouseParametersTablePath);
runGroupsParametersTable = mouseParametersTable.runGroupsParametersTable;

postTextSuptitle = '';
analyzedCells = unique(cell2mat(mouseAnalysisTable.cellName));

groupByCell = 0;

cellsConnectionInfo = [];
cellsConnectionInfo.blue = [];
cellsConnectionInfo.red = [];
cellsConnectionInfo.blue.nonConnectedCells = [];
cellsConnectionInfo.blue.connectedCells = [];
cellsConnectionInfo.red.nonConnectedCells = [];
cellsConnectionInfo.red.connectedCells = [];
cellsConnectionInfo.doubleConnectedCells = [];
cellsConnectionInfo.responseType = [];

for iCell = 1:numel(analyzedCells)
    
iAnalyzedCell = analyzedCells(iCell);

cellName = ['cell' num2str(iAnalyzedCell)];
cellsConnectionInfo.blue.(cellName) = [];
cellsConnectionInfo.red.(cellName) = [];

initializeFig(0.9,0.9);

nTileRows = 2;
nTileColumns = 7;
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
cat1_selectedCellsIDs = iAnalyzedCell;
cat1_desiredEpochs = [];
cat1_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'blue',responseTypeBlue);

% No need to modify
cat1_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat1_selectedCellsIDs);
cat1_selectedRunGroupsAnalysisTable = cat1_selectedMembersTable(ismember(cat1_selectedMembersTable.runGroup,cat1_selectedRunGroups),:);
cat1_selectedRunGroupsAnalysisTable = extractEpochSubset(cat1_selectedRunGroupsAnalysisTable, cat1_desiredEpochs);

plotData = generatePlotDataStruct;
plotData.queryRunFeature.cat1 = cat1_queryRunFeature;
plotData.selectedRunGroupsAnalysisTable.cat1 = cat1_selectedRunGroupsAnalysisTable;

%% Plot
nexttile(masterLayout,1,[1 1]); axis off;

options = generateOptionsStruct(groupByCell);
options.plotColor.cat1 = colorBlue;
options.plotTitle = ['Recording quality'];
options.xLabelText.cat1 = ['Acquisition'];
options.yLabelText = ['Series resistance [M\Omega]'];

generateDotplot(plotData,options)
    
%% ** Plot blue traces **

%% Collect data

cat1_selectedCellsIDs = iAnalyzedCell;
cat1_desiredEpochs = [];
cat1_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'blue',responseTypeBlue);

% No need to modify
cat1_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat1_selectedCellsIDs);
cat1_selectedRunGroupsAnalysisTable = cat1_selectedMembersTable(ismember(cat1_selectedMembersTable.runGroup,cat1_selectedRunGroups),:);
cat1_selectedRunGroupsAnalysisTable = extractEpochSubset(cat1_selectedRunGroupsAnalysisTable, cat1_desiredEpochs);

plotData = generatePlotDataStruct;
plotData.selectedRunGroupsAnalysisTable.cat1 = cat1_selectedRunGroupsAnalysisTable;

%% Plot
nexttile(masterLayout,2,[1 2]); axis off;

options = generateOptionsStruct(groupByCell);
options.plotColor.cat1 = colorBlue;
options.plotTitle = ['Full-field responses'];

generatePlotTraces(plotData,options)

%% ** Barplot blue heightControlPeak vs heightPulsePeak **

%% Collect data

cat1_queryRunFeature = 'heightControlPeak';
cat2_queryRunFeature = 'heightPulsePeak';
cat1_selectedCellsIDs = iAnalyzedCell;
cat2_selectedCellsIDs = iAnalyzedCell;
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
options.plotTitle = ['Response amplitude'];
options.xLabelText.cat1 = ['Control'];
options.xLabelText.cat2 = ['Light'];
options.yLabelText = ['Current [pA]'];

generateBarplot(plotData,options)
cellsConnectionInfo = determineIfConnected(cellsConnectionInfo,cellName,'blue',plotData,options);

%% ** Barplot blue areaControl vs areaPulse **

%% Collect data

cat1_queryRunFeature = 'areaControl';
cat2_queryRunFeature = 'areaPulse';
cat1_selectedCellsIDs = iAnalyzedCell;
cat2_selectedCellsIDs = iAnalyzedCell;
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
options.plotTitle = ['Response charge'];
options.xLabelText.cat1 = ['Control'];
options.xLabelText.cat2 = ['Light'];
options.yLabelText = ['Charge [pC]'];

generateBarplot(plotData,options)
cellsConnectionInfo = determineIfConnected(cellsConnectionInfo,cellName,'blue',plotData,options);

%% ** Barplot blue timePulsePeak **

%% Collect data

cat1_queryRunFeature = 'timePulsePeak';
cat1_selectedCellsIDs = iAnalyzedCell;
cat1_desiredEpochs = [];
cat1_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'blue',responseTypeBlue);

% No need to modify
cat1_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat1_selectedCellsIDs);
cat1_selectedRunGroupsAnalysisTable = cat1_selectedMembersTable(ismember(cat1_selectedMembersTable.runGroup,cat1_selectedRunGroups),:);
cat1_selectedRunGroupsAnalysisTable = extractEpochSubset(cat1_selectedRunGroupsAnalysisTable, cat1_desiredEpochs);

plotData = generatePlotDataStruct;
plotData.queryRunFeature.cat1 = cat1_queryRunFeature;
plotData.selectedRunGroupsAnalysisTable.cat1 = cat1_selectedRunGroupsAnalysisTable;

%% Plot
nexttile(masterLayout,6,[1 1]); axis off;

options = generateOptionsStruct(groupByCell);
options.plotColor.cat1 = colorBlue;
options.plotColorOverlay.cat1 = colorBlueOverlay;
options.plotTitle = ['Response delay'];
options.xLabelText.cat1 = ['Light'];
options.yLabelText = ['Time [ms]'];

generateBarplot(plotData,options)

%% ** Barplot blue different amplitudes **

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

queryParameters = []; queryParameters.amplitudeBlue = 1;
cat1_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'blue',responseTypeBlue,queryParameters);
queryParameters = []; queryParameters.amplitudeBlue = 3;
cat2_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'blue',responseTypeBlue,queryParameters);
queryParameters = []; queryParameters.amplitudeBlue = 5;
cat3_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'blue',responseTypeBlue,queryParameters);

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
nexttile(masterLayout,7,[1 1]); axis off;

options = generateOptionsStruct(groupByCell);
options.plotColor.cat1 = colorBlue;
options.plotColor.cat2 = colorBlue;
options.plotColor.cat3 = colorBlue;
options.plotColorOverlay.cat1 = colorBlueOverlay;
options.plotColorOverlay.cat2 = colorBlueOverlay;
options.plotColorOverlay.cat3 = colorBlueOverlay;
options.plotTitle = ['Different amplitudes'];
options.xLabelText.cat1 = ['1 mV'];
options.xLabelText.cat2 = ['3 mV'];
options.xLabelText.cat3 = ['5 mV'];
options.yLabelText = ['Current [pA]'];

generateBarplot(plotData,options)

%% ** Red quality RS0 **

%% Collect data

cat1_queryRunFeature = 'qualityRS0';
cat1_selectedCellsIDs = iAnalyzedCell;
cat1_desiredEpochs = [];
cat1_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'red',responseTypeRed);

% No need to modify
cat1_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat1_selectedCellsIDs);
cat1_selectedRunGroupsAnalysisTable = cat1_selectedMembersTable(ismember(cat1_selectedMembersTable.runGroup,cat1_selectedRunGroups),:);
cat1_selectedRunGroupsAnalysisTable = extractEpochSubset(cat1_selectedRunGroupsAnalysisTable, cat1_desiredEpochs);

plotData = generatePlotDataStruct;
plotData.queryRunFeature.cat1 = cat1_queryRunFeature;
plotData.selectedRunGroupsAnalysisTable.cat1 = cat1_selectedRunGroupsAnalysisTable;

%% Plot
nexttile(masterLayout,nTileColumns+1,[1 1]); axis off;

options = generateOptionsStruct(groupByCell);
options.plotColor.cat1 = colorRed;
options.plotTitle = ['Recording quality'];
options.xLabelText.cat1 = ['Acquisition'];
options.yLabelText = ['Series resistance [M\Omega]'];

generateDotplot(plotData,options)

%% ** Plot red traces **

%% Collect data

cat1_selectedCellsIDs = iAnalyzedCell;
cat1_desiredEpochs = [];
cat1_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'red',responseTypeRed);

% No need to modify
cat1_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat1_selectedCellsIDs);
cat1_selectedRunGroupsAnalysisTable = cat1_selectedMembersTable(ismember(cat1_selectedMembersTable.runGroup,cat1_selectedRunGroups),:);
cat1_selectedRunGroupsAnalysisTable = extractEpochSubset(cat1_selectedRunGroupsAnalysisTable, cat1_desiredEpochs);

plotData = generatePlotDataStruct;
plotData.selectedRunGroupsAnalysisTable.cat1 = cat1_selectedRunGroupsAnalysisTable;

%% Plot
nexttile(masterLayout,nTileColumns+2,[1 2]); axis off;

options = generateOptionsStruct(groupByCell);
options.plotColor.cat1 = colorRed;
options.plotTitle = ['Full-field responses'];

generatePlotTraces(plotData,options)

%% ** Barplot red heightControlPeak vs heightPulsePeak **

%% Collect data

cat1_queryRunFeature = 'heightControlPeak';
cat2_queryRunFeature = 'heightPulsePeak';
cat1_selectedCellsIDs = iAnalyzedCell;
cat2_selectedCellsIDs = iAnalyzedCell;
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
nexttile(masterLayout,nTileColumns+4,[1 1]); axis off;

options = generateOptionsStruct(groupByCell);
options.plotColor.cat1 = colorGrey;
options.plotColor.cat2 = colorRed;
options.plotColorOverlay.cat1 = colorGreyOverlay;
options.plotColorOverlay.cat2 = colorRedOverlay;
options.plotTitle = ['Response amplitude'];
options.xLabelText.cat1 = ['Control'];
options.xLabelText.cat2 = ['Light'];
options.yLabelText = ['Current [pA]'];

generateBarplot(plotData,options)
cellsConnectionInfo = determineIfConnected(cellsConnectionInfo,cellName,'red',plotData,options);

%% ** Barplot red areaControl vs areaPulse **

%% Collect data

cat1_queryRunFeature = 'areaControl';
cat2_queryRunFeature = 'areaPulse';
cat1_selectedCellsIDs = iAnalyzedCell;
cat2_selectedCellsIDs = iAnalyzedCell;
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
nexttile(masterLayout,nTileColumns+5,[1 1]); axis off;

options = generateOptionsStruct(groupByCell);
options.plotColor.cat1 = colorGrey;
options.plotColor.cat2 = colorRed;
options.plotColorOverlay.cat1 = colorGreyOverlay;
options.plotColorOverlay.cat2 = colorRedOverlay;
options.plotTitle = ['Response charge'];
options.xLabelText.cat1 = ['Control'];
options.xLabelText.cat2 = ['Light'];
options.yLabelText = ['Charge [pC]'];

generateBarplot(plotData,options)
cellsConnectionInfo = determineIfConnected(cellsConnectionInfo,cellName,'red',plotData,options);

%% ** Barplot red timePulsePeak **

%% Collect data

cat1_queryRunFeature = 'timePulsePeak';
cat1_selectedCellsIDs = iAnalyzedCell;
cat1_desiredEpochs = [];
cat1_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'red',responseTypeRed);

% No need to modify
cat1_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat1_selectedCellsIDs);
cat1_selectedRunGroupsAnalysisTable = cat1_selectedMembersTable(ismember(cat1_selectedMembersTable.runGroup,cat1_selectedRunGroups),:);
cat1_selectedRunGroupsAnalysisTable = extractEpochSubset(cat1_selectedRunGroupsAnalysisTable, cat1_desiredEpochs);

plotData = generatePlotDataStruct;
plotData.queryRunFeature.cat1 = cat1_queryRunFeature;
plotData.selectedRunGroupsAnalysisTable.cat1 = cat1_selectedRunGroupsAnalysisTable;

%% Plot
nexttile(masterLayout,nTileColumns+6,[1 1]); axis off;

options = generateOptionsStruct(groupByCell);
options.plotColor.cat1 = colorRed;
options.plotColorOverlay.cat1 = colorRedOverlay;
options.plotTitle = ['Response delay'];
options.xLabelText.cat1 = ['Light'];
options.yLabelText = ['Time [ms]'];

generateBarplot(plotData,options)

%% ** Barplot red different amplitudes **

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

queryParameters = []; queryParameters.amplitudeRed = 1;
cat1_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'red',responseTypeRed,queryParameters);
queryParameters = []; queryParameters.amplitudeRed = 3;
cat2_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'red',responseTypeRed,queryParameters);
queryParameters = []; queryParameters.amplitudeRed = 5;
cat3_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,'red',responseTypeRed,queryParameters);

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
nexttile(masterLayout,nTileColumns+7,[1 1]); axis off;

options = generateOptionsStruct(groupByCell);
options.plotColor.cat1 = colorRed;
options.plotColor.cat2 = colorRed;
options.plotColor.cat3 = colorRed;
options.plotColorOverlay.cat1 = colorRedOverlay;
options.plotColorOverlay.cat2 = colorRedOverlay;
options.plotColorOverlay.cat3 = colorRedOverlay;
options.plotTitle = ['Different amplitudes'];
options.xLabelText.cat1 = ['1 mV'];
options.xLabelText.cat2 = ['3 mV'];
options.xLabelText.cat3 = ['5 mV'];
options.yLabelText = ['Current [pA]'];

generateBarplot(plotData,options)

%% Cell connection

if any(cell2mat(cellsConnectionInfo.blue.(cellName))==0); cellsConnectionInfo.blue.nonConnectedCells{end+1} = iAnalyzedCell;
elseif all(cell2mat(cellsConnectionInfo.blue.(cellName))==1); cellsConnectionInfo.blue.connectedCells{end+1} = iAnalyzedCell; end

if any(cell2mat(cellsConnectionInfo.red.(cellName))==0); cellsConnectionInfo.red.nonConnectedCells{end+1} = iAnalyzedCell;
elseif all(cell2mat(cellsConnectionInfo.red.(cellName))==1); cellsConnectionInfo.red.connectedCells{end+1} = iAnalyzedCell; end

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

cellsConnectionInfo.doubleConnectedCells = intersect(cell2mat(cellsConnectionInfo.blue.connectedCells),cell2mat(cellsConnectionInfo.red.connectedCells));
cellsConnectionInfo.responseType.blue = responseTypeBlue;
cellsConnectionInfo.responseType.red = responseTypeRed;
disp('Blue connected cells: '); cell2mat(cellsConnectionInfo.blue.connectedCells)
disp('Red connected cells: '); cell2mat(cellsConnectionInfo.red.connectedCells)
disp('Double-connected cells: '); cellsConnectionInfo.doubleConnectedCells

saveCellsConnectionInfoName = 'CellsConnectionInfo.mat';
saveCellsConnectionInfoName = strrep(saveCellsConnectionInfoName, ' ', '');
saveCellsConnectionInfoPath = fullfile(savePlotPath, saveCellsConnectionInfoName);

if ~isfile(saveCellsConnectionInfoPath) || overwriteConnectionInfo == 1
     save(saveCellsConnectionInfoPath,'cellsConnectionInfo');
end
