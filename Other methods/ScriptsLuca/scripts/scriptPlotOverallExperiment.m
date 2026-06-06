disp(['*-*-*-* Running: scriptPlotOverallExperiment *-*-*-*'])

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

masterLayout = tiledlayout(2,4);
masterLayout.TileSpacing = 'compact';
masterLayout.Padding = 'compact';

colorBlue = [85, 161, 254]./255;
colorRed = [255, 50, 58]./255;
colorGrey = [192, 192, 192]./255;

colorBlueOverlay = [85, 161, 254]./255 * 0.6;
colorRedOverlay = [255, 50, 58]./255 * 0.6;
colorGreyOverlay = [192, 192, 192]./255 * 0.6;

groupByCell = 1;
isOverallExperimentPlot = 1;

nonConnectedCellsBlue = cell2mat(experimentCellsConnectionInfo.blue.nonConnectedCells);
connectedCellsBlue = cell2mat(experimentCellsConnectionInfo.blue.connectedCells);
nonConnectedCellsRed = cell2mat(experimentCellsConnectionInfo.red.nonConnectedCells);
connectedCellsRed = cell2mat(experimentCellsConnectionInfo.red.connectedCells);

nonConnectedCellsExcitatory = cell2mat(experimentCellsConnectionInfo.excitatory.nonConnectedCells);
connectedCellsExcitatory = cell2mat(experimentCellsConnectionInfo.excitatory.connectedCells);
nonConnectedCellsInhibitory = cell2mat(experimentCellsConnectionInfo.inhibitory.nonConnectedCells);
connectedCellsInhibitory = cell2mat(experimentCellsConnectionInfo.inhibitory.connectedCells);

if strcmp(channel,'blue') 
    nonConnectedCells = nonConnectedCellsBlue; 
    connectedCells = connectedCellsBlue; 
    plotColorConnected = colorBlue; 
    plotColorOverlayConnected = colorBlueOverlay;
elseif strcmp(channel,'red') 
    nonConnectedCells = nonConnectedCellsRed; 
    connectedCells = connectedCellsRed; 
    plotColorConnected = colorRed; 
    plotColorOverlayConnected = colorRedOverlay; 
elseif strcmp(channel,'blue|red')
    if strcmp(responseType,'excitatory')
        nonConnectedCells = nonConnectedCellsExcitatory; 
        connectedCells = connectedCellsExcitatory; 
        plotColorConnected = colorBlue; 
        plotColorOverlayConnected = colorBlueOverlay;
    elseif strcmp(responseType,'inhibitory')
        nonConnectedCells = nonConnectedCellsInhibitory; 
        connectedCells = connectedCellsInhibitory; 
        plotColorConnected = colorRed; 
        plotColorOverlayConnected = colorRedOverlay;        
    end    
end
    
%% ** Plot traces **

%% Collect data

cat1_selectedCellsIDs = nonConnectedCells;
cat2_selectedCellsIDs = connectedCells;
cat1_desiredEpochs = [];
cat2_desiredEpochs = [];
cat1_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,channel,responseType,queryParameters);
cat2_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,channel,responseType,queryParameters);

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

nNonConnectedCells = numel(unique(cell2mat(cat1_selectedRunGroupsAnalysisTable.overallCellName)));
nConnectedCells = numel(unique(cell2mat(cat2_selectedRunGroupsAnalysisTable.overallCellName)));

%% Plot
nexttile(masterLayout,1,[1 2]); axis off;

options = generateOptionsStruct(groupByCell,isOverallExperimentPlot);
options.plotColor.cat1 = colorGrey;
options.plotColor.cat2 = plotColorConnected;
options.plotTitle = ['Full-field responses'];
options.legendText.cat1 = ['Non-connected'];
options.legendText.cat2 = ['Connected'];

generatePlotTraces(plotData,options)

%% ** Barplot heightControlPeak vs heightPulsePeak **

%% Collect data

cat1_queryRunFeature = 'heightControlPeak';
cat2_queryRunFeature = 'heightPulsePeak';
cat1_selectedCellsIDs = connectedCells;
cat2_selectedCellsIDs = connectedCells;
cat1_desiredEpochs = [];
cat2_desiredEpochs = [];
cat1_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,channel,responseType,queryParameters);
cat2_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,channel,responseType,queryParameters);

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
options.plotColor.cat2 = plotColorConnected;
options.plotColorOverlay.cat1 = colorGreyOverlay;
options.plotColorOverlay.cat2 = plotColorOverlayConnected;
options.plotTitle = ['Response amplitude' newline 'Connected cells'];
options.xLabelText.cat1 = ['Control'];
options.xLabelText.cat2 = ['Light'];
options.yLabelText = ['Current [pA]'];

generateBarplot(plotData,options)

%% ** Barplot areaControl vs areaPulse **

%% Collect data

cat1_queryRunFeature = 'areaControl';
cat2_queryRunFeature = 'areaPulse';
cat1_selectedCellsIDs = connectedCells;
cat2_selectedCellsIDs = connectedCells;
cat1_desiredEpochs = [];
cat2_desiredEpochs = [];
cat1_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,channel,responseType,queryParameters);
cat2_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,channel,responseType,queryParameters);

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
nexttile(masterLayout,6,[1 1]); axis off;

options = generateOptionsStruct(groupByCell,isOverallExperimentPlot);
options.plotColor.cat1 = colorGrey;
options.plotColor.cat2 = plotColorConnected;
options.plotColorOverlay.cat1 = colorGreyOverlay;
options.plotColorOverlay.cat2 = plotColorOverlayConnected;
options.plotTitle = ['Response charge' newline 'Connected cells'];
options.xLabelText.cat1 = ['Control'];
options.xLabelText.cat2 = ['Light'];
options.yLabelText = ['Charge [pC]'];

generateBarplot(plotData,options)

%% ** Barplot non-connected heightControlPeak vs heightPulsePeak **

% %% Collect data
% 
% cat1_queryRunFeature = 'heightControlPeak';
% cat2_queryRunFeature = 'heightPulsePeak';
% cat1_selectedCellsIDs = nonConnectedCells;
% cat2_selectedCellsIDs = nonConnectedCells;
% cat1_desiredEpochs = [];
% cat2_desiredEpochs = [];
% cat1_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,channel,responseType,queryParameters);
% cat2_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,channel,responseType,queryParameters);
% 
% % No need to modify
% cat1_selectedMembersTable = extractSelectedMembersFromTable(experimentAnalysisTable, 'overallCellName', cat1_selectedCellsIDs);
% cat1_selectedRunGroupsAnalysisTable = cat1_selectedMembersTable(ismember(cat1_selectedMembersTable.runGroup,cat1_selectedRunGroups),:);
% cat1_selectedRunGroupsAnalysisTable = extractEpochSubset(cat1_selectedRunGroupsAnalysisTable, cat1_desiredEpochs);
% 
% cat2_selectedMembersTable = extractSelectedMembersFromTable(experimentAnalysisTable, 'overallCellName', cat2_selectedCellsIDs);
% cat2_selectedRunGroupsAnalysisTable = cat2_selectedMembersTable(ismember(cat2_selectedMembersTable.runGroup,cat2_selectedRunGroups),:);
% cat2_selectedRunGroupsAnalysisTable = extractEpochSubset(cat2_selectedRunGroupsAnalysisTable, cat2_desiredEpochs);
% 
% plotData = generatePlotDataStruct;
% plotData.queryRunFeature.cat1 = cat1_queryRunFeature;
% plotData.queryRunFeature.cat2 = cat2_queryRunFeature;
% plotData.selectedRunGroupsAnalysisTable.cat1 = cat1_selectedRunGroupsAnalysisTable;
% plotData.selectedRunGroupsAnalysisTable.cat2 = cat2_selectedRunGroupsAnalysisTable;
% 
% %% Plot
% nexttile(masterLayout,7,[1 1]); axis off;
% 
% options = generateOptionsStruct(groupByCell,isOverallExperimentPlot);
% options.plotColor.cat1 = colorGrey;
% options.plotColor.cat2 = plotColorConnected;
% options.plotColorOverlay.cat1 = colorGreyOverlay;
% options.plotColorOverlay.cat2 = plotColorOverlayConnected;
% options.plotTitle = ['Response amplitude' newline 'Non-connected cells'];
% options.xLabelText.cat1 = ['Control'];
% options.xLabelText.cat2 = ['Light'];
% options.yLabelText = ['Current [pA]'];
% 
% generateBarplot(plotData,options)

% ** Barplot non-connected heightControlPeak vs heightPulsePeak **

%% Collect data

cat1_queryRunFeature = 'timePulsePeak';
cat1_selectedCellsIDs = connectedCells;
cat1_desiredEpochs = [];
cat1_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,channel,responseType,queryParameters);

% No need to modify
cat1_selectedMembersTable = extractSelectedMembersFromTable(experimentAnalysisTable, 'overallCellName', cat1_selectedCellsIDs);
cat1_selectedRunGroupsAnalysisTable = cat1_selectedMembersTable(ismember(cat1_selectedMembersTable.runGroup,cat1_selectedRunGroups),:);
cat1_selectedRunGroupsAnalysisTable = extractEpochSubset(cat1_selectedRunGroupsAnalysisTable, cat1_desiredEpochs);

plotData = generatePlotDataStruct;
plotData.queryRunFeature.cat1 = cat1_queryRunFeature;
plotData.selectedRunGroupsAnalysisTable.cat1 = cat1_selectedRunGroupsAnalysisTable;

%% Plot
nexttile(masterLayout,7,[1 1]); axis off;

options = generateOptionsStruct(groupByCell,isOverallExperimentPlot);
options.plotColor.cat1 = plotColorConnected;
options.plotColorOverlay.cat1 = plotColorOverlayConnected;
options.plotTitle = ['Response delay' newline 'Connected cells'];
options.xLabelText.cat1 = ['Light'];
options.yLabelText = ['Time [ms]'];

generateBarplot(plotData,options)

%% ** Barplot non-connected areaControl vs areaPulse **

%% Collect data

cat1_queryRunFeature = 'areaControl';
cat2_queryRunFeature = 'areaPulse';
cat1_selectedCellsIDs = nonConnectedCells;
cat2_selectedCellsIDs = nonConnectedCells;
cat1_desiredEpochs = [];
cat2_desiredEpochs = [];
cat1_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,channel,responseType,queryParameters);
cat2_selectedRunGroups = findSelectedRunGroups(runGroupsParametersTable,channel,responseType,queryParameters);

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
options.plotColor.cat1 = colorGrey;
options.plotColor.cat2 = plotColorConnected;
options.plotColorOverlay.cat1 = colorGreyOverlay;
options.plotColorOverlay.cat2 = plotColorOverlayConnected;
options.plotTitle = ['Response charge' newline 'Non-connected cells'];
options.xLabelText.cat1 = ['Control'];
options.xLabelText.cat2 = ['Light'];
options.yLabelText = ['Charge [pC]'];

generateBarplot(plotData,options)

%% ** Pie chart non-connected vs connected cell number **
nexttile(masterLayout,3,[1 2]); axis off;

% nConnectedCells = numel(connectedCells);
% nNonConnectedCells = numel(nonConnectedCells);
nTotalCells = nConnectedCells + nNonConnectedCells;

connectedCellsLabel = ['Connected cells (' num2str(nConnectedCells) '/' num2str(nTotalCells) ')'];
nonConnectedCellsLabel = ['Non-connected cells (' num2str(nNonConnectedCells) '/' num2str(nTotalCells) ')'];
chartCategories = {connectedCellsLabel,nonConnectedCellsLabel};
chartNumbers = [nConnectedCells nNonConnectedCells];
pieChart = pie(chartNumbers,chartCategories);

pConnected = pieChart(1);
pConnected.FaceColor = plotColorConnected;
tConnected = pieChart(2);
tConnected.FontSize = 14;
pNonConnected = pieChart(3);
pNonConnected.FaceColor = colorGrey;
tNonConnected = pieChart(4);
tNonConnected.FontSize = 14;

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

