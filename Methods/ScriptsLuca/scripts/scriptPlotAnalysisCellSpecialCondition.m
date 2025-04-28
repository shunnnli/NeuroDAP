disp(['*-*-*-* Running: scriptPlotAnalysisCellSpecialCondition *-*-*-*'])

mouseAnalysisTablePath = [experimentDirectory filesep 'mouseAnalysis' filesep 'MouseAnalysisTable.mat'];
mouseParametersTablePath = [experimentDirectory filesep 'mouseAnalysis' filesep 'MouseParametersTable.mat'];

mouseAnalysisTable = load(mouseAnalysisTablePath);
mouseAnalysisTable = mouseAnalysisTable.mouseAnalysisTable;
mouseParametersTable = load(mouseParametersTablePath);
runGroupsParametersTable = mouseParametersTable.runGroupsParametersTable;

postTextSuptitle = ['_','Temperature'];
analyzedCells = unique(cell2mat(mouseAnalysisTable.cellName));

groupByCell = 0;
channel = 'blue';   
iAnalyzedCell = ;

cellName = ['cell' num2str(iAnalyzedCell)];
cellsConnectionInfo.blue.(cellName) = [];
cellsConnectionInfo.red.(cellName) = [];

initializeFig(0.9,0.9);

nTileRows = 1;
nTileColumns = 3;
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
    
%% ** Plot traces **

%% Collect data

cat1_selectedCellsIDs = iAnalyzedCell;
cat2_selectedCellsIDs = iAnalyzedCell;

cat1_desiredEpochs = [];
cat2_desiredEpochs = [];

queryParameters = [];
queryParameters.cellTemperature = 24;
cat1_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,channel,responseTypeBlue,queryParameters);

queryParameters = [];
queryParameters.cellTemperature = 32;
cat2_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,channel,responseTypeBlue,queryParameters);

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
nexttile(masterLayout,1,[1 1]); axis off;

options = generateOptionsStruct(groupByCell);
options.plotColor.cat1 = colorGrey;
options.plotColor.cat2 = colorBlue;
options.plotTitle = ['Full-field responses'];
options.legendText.cat1 = ['24°C'];
options.legendText.cat1 = ['32°C'];

generatePlotTraces(plotData,options)

%% ** Barplot blue heightPulsePeak **

%% Collect data

cat1_queryRunFeature = 'heightPulsePeak';
cat2_queryRunFeature = 'heightPulsePeak';
cat1_selectedCellsIDs = iAnalyzedCell;
cat2_selectedCellsIDs = iAnalyzedCell;
cat1_desiredEpochs = [];
cat2_desiredEpochs = [];

queryParameters = [];
queryParameters.cellTemperature = 24;
cat1_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,channel,responseTypeBlue,queryParameters);

queryParameters = [];
queryParameters.cellTemperature = 32;
cat2_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,channel,responseTypeBlue,queryParameters);

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
nexttile(masterLayout,2,[1 1]); axis off;

options = generateOptionsStruct(groupByCell);
options.plotColor.cat1 = colorGrey;
options.plotColor.cat2 = colorBlue;
options.plotColorOverlay.cat1 = colorGreyOverlay;
options.plotColorOverlay.cat2 = colorBlueOverlay;
options.plotTitle = ['Response amplitude'];
options.xLabelText.cat1 = ['24°C'];
options.xLabelText.cat2 = ['32°C'];
options.yLabelText = ['Current [pA]'];

generateBarplot(plotData,options)

%% ** Barplot blue areaPulse **

%% Collect data

cat1_queryRunFeature = 'areaPulse';
cat2_queryRunFeature = 'areaPulse';
cat1_selectedCellsIDs = iAnalyzedCell;
cat2_selectedCellsIDs = iAnalyzedCell;
cat1_desiredEpochs = [];
cat2_desiredEpochs = [];

queryParameters = [];
queryParameters.cellTemperature = 24;
cat1_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,channel,responseTypeBlue,queryParameters);

queryParameters = [];
queryParameters.cellTemperature = 32;
cat2_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,channel,responseTypeBlue,queryParameters);

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
options.plotTitle = ['Response charge'];
options.xLabelText.cat1 = ['24°C'];
options.xLabelText.cat2 = ['32°C'];
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
