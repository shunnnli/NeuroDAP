clear all;
close all;
clc;

% User environment (automatic)
user = getenv('USER'); 

% Indicate main directory and add it to path
matlabDirectory = ['/Users/' user '/Desktop/Thesis/Matlab/'];
mainDirectory = ['/Users/' user '/Desktop/Thesis/Matlab/Experiments/']; % including all Matlab toolboxes and necessary paths
addpath(genpath(matlabDirectory));

% Indicate experiment to analyze
experimentName = 'vGAT_JAWS_Ephys_7_240401';

% Set path
experimentDirectory = [mainDirectory experimentName];
analyzedMouseDataFolderName = 'mouseAnalysis';
analyzedMouseDataFolderPath = [experimentDirectory filesep analyzedMouseDataFolderName];
mouseTableName = 'MouseAnalysisTable.mat';
mouseTableName = strrep(mouseTableName, ' ', '');
mouseParametersTableName = 'MouseParametersTable.mat';
mouseParametersTableName = strrep(mouseParametersTableName, ' ', '');
mouseAnalysisTablePath = fullfile(analyzedMouseDataFolderPath, mouseTableName);
mouseParametersTablePath = fullfile(analyzedMouseDataFolderPath, mouseParametersTableName);

mouseAnalysisTable = load(mouseAnalysisTablePath);
mouseAnalysisTable = mouseAnalysisTable.mouseAnalysisTable;
mouseParametersTable = load(mouseParametersTablePath);
runGroupsParametersTable = mouseParametersTable.runGroupsParametersTable;

typeCC = 'comparisonWithWithoutLight';
savePlotPath = [experimentDirectory filesep 'mouseAnalysis'];

%%

queryOptoParameter =  'amplitudeCurrentPulse';
queryRunFeature = 'areaVtPulse';

selectedCellsIDs = [2,3];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0') & cell2mat(runGroupsParametersTable.widthCurrentPulse) == 6000);
selectedRunGroupsAnalysisTable = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroups),:);

queryOptoParameterValues = zeros(numel(selectedRunGroups),1);
queryOptoParameterArray = [];
queryRunFeatureData = {};
cellNameData= {};

groupsCounter = 1;

for iRunGroup = 1:numel(selectedRunGroups)
    
    rowsRunGroup = find(selectedRunGroupsAnalysisTable.runGroup(:) == selectedRunGroups(iRunGroup));
    queryOptoParameterValues(groupsCounter) = selectedRunGroupsAnalysisTable.optoParameters{rowsRunGroup(1)}.(queryOptoParameter);
    
    queryOptoParameterArray = [queryOptoParameterArray; repmat(queryOptoParameterValues(groupsCounter),numel(rowsRunGroup),1)];
    queryRunFeatureData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTable.(queryRunFeature)(rowsRunGroup));
    cellNameData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTable.cellName(rowsRunGroup));
    
    groupsCounter = groupsCounter + 1;
    
end

allData = vertcat(queryRunFeatureData{:});
cellNameDataArray = vertcat(cellNameData{:});

meansQueryRunFeatureDataAllCells = nan(numel(unique(cellNameDataArray)),numel(unique(queryOptoParameterValues)));


for iQueryOptoParameter = 1:length(queryRunFeatureData)

    cellNames = unique(cellNameData{iQueryOptoParameter});
    
    for iCell = 1:numel(cellNames)
        
        cellRows = find(cellNameData{iQueryOptoParameter} == cellNames(iCell)); 
        meanQueryRunFeatureDataCell = nanmean(queryRunFeatureData{iQueryOptoParameter}(cellRows));
        positionOptoParameterValue = find(sort(queryOptoParameterValues) == queryOptoParameterValues(iQueryOptoParameter));
        plotscatter = scatter(positionOptoParameterValue, meanQueryRunFeatureDataCell, 40, 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0.5 0.5 0.5]);
        plotscatter.Annotation.LegendInformation.IconDisplayStyle = 'off';
        hold on;
        
        meansQueryRunFeatureDataAllCells(iCell,positionOptoParameterValue) = meanQueryRunFeatureDataCell;
        
    end
    
end

meanQueryRunFeatureDataAllCells = nanmean(meansQueryRunFeatureDataAllCells,1);
stdQueryRunFeatureDataAllCells = nanstd(meansQueryRunFeatureDataAllCells,1);
xValues = linspace(1,numel(queryOptoParameterValues),numel(queryOptoParameterValues));
plot(xValues, meanQueryRunFeatureDataAllCells,'Color',[0 0 0],'LineWidth',2);
ploterrorbar = errorbar(xValues,meanQueryRunFeatureDataAllCells,stdQueryRunFeatureDataAllCells,'.','Color',[0 0 0],'LineWidth', 1.2, 'CapSize', 20);
ploterrorbar.Annotation.LegendInformation.IconDisplayStyle = 'off';
hold on;

hScatter = findobj(gca, 'Type', 'scatter');
uistack(hScatter, 'top');

xticks(xValues)
xticklabels([unique(queryOptoParameterValues)])
xlim([0.5, 9.5])
xlabel(queryOptoParameter);
ylabel(queryRunFeature);
box on;
title([queryOptoParameter,' vs ', queryRunFeature]);
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');

traceNoLight = selectedRunGroupsAnalysisTable.trace{51,1};

%% amplitudeCurrentPulse vs areaVtPulse -- withLightBeforeCurrent -- line plot

queryOptoParameter =  'amplitudeCurrentPulse';
queryRunFeature = 'areaVtPulse';

selectedCellsIDs = [2,3];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao2') & cell2mat(runGroupsParametersTable.widthCurrentPulse) == 6000 & cell2mat(runGroupsParametersTable.delayLightPulse) == 19000 & cell2mat(runGroupsParametersTable.amplitudeLightPulse) == 5);
selectedRunGroupsAnalysisTable = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroups),:);

queryOptoParameterValues = zeros(numel(selectedRunGroups),1);
queryOptoParameterArray = [];
queryRunFeatureData = {};
cellNameData= {};

groupsCounter = 1;

for iRunGroup = 1:numel(selectedRunGroups)
    
    rowsRunGroup = find(selectedRunGroupsAnalysisTable.runGroup(:) == selectedRunGroups(iRunGroup));
    queryOptoParameterValues(groupsCounter) = selectedRunGroupsAnalysisTable.optoParameters{rowsRunGroup(1)}.(queryOptoParameter);
    
    queryOptoParameterArray = [queryOptoParameterArray; repmat(queryOptoParameterValues(groupsCounter),numel(rowsRunGroup),1)];
    queryRunFeatureData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTable.(queryRunFeature)(rowsRunGroup));
    cellNameData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTable.cellName(rowsRunGroup));
    
    groupsCounter = groupsCounter + 1;
    
end

allData = vertcat(queryRunFeatureData{:});
cellNameDataArray = vertcat(cellNameData{:});

meansQueryRunFeatureDataAllCells = nan(numel(unique(cellNameDataArray)),numel(unique(queryOptoParameterValues)));

for iQueryOptoParameter = 1:length(queryRunFeatureData)

    cellNames = unique(cellNameData{iQueryOptoParameter});
    
    for iCell = 1:numel(cellNames)
        
        cellRows = find(cellNameData{iQueryOptoParameter} == cellNames(iCell)); 
        meanQueryRunFeatureDataCell = nanmean(queryRunFeatureData{iQueryOptoParameter}(cellRows));
        positionOptoParameterValue = find(sort(queryOptoParameterValues) == queryOptoParameterValues(iQueryOptoParameter));
        plotscatter = scatter(positionOptoParameterValue, meanQueryRunFeatureDataCell, 40, 'MarkerFaceColor', [0 0.4 0.7], 'MarkerEdgeColor', [0 0.4 0.7]);
        %plotscatter.Annotation.LegendInformation.IconDisplayStyle = 'off';
        hold on;
        
        meansQueryRunFeatureDataAllCells(iCell,positionOptoParameterValue) = meanQueryRunFeatureDataCell;
        
    end
    
end

meanQueryRunFeatureDataAllCells = nanmean(meansQueryRunFeatureDataAllCells,1);
stdQueryRunFeatureDataAllCells = nanstd(meansQueryRunFeatureDataAllCells,1);
xValues = linspace(1,numel(queryOptoParameterValues),numel(queryOptoParameterValues));
plot(xValues, meanQueryRunFeatureDataAllCells,'Color',[0 0.4470 0.7410],'LineWidth',2);
ploterrorbar = errorbar(xValues,meanQueryRunFeatureDataAllCells,stdQueryRunFeatureDataAllCells,'.','Color',[0 0 0],'LineWidth', 1.2, 'CapSize', 20);

hScatter = findobj(gca, 'Type', 'scatter');
uistack(hScatter, 'top');

xticks(xValues)
xticklabels([unique(queryOptoParameterValues)])
xlabel(queryOptoParameter);
ylabel(queryRunFeature);
box on;
title([queryOptoParameter,' vs ', queryRunFeature]);
legend('no Light','with Light','Location', 'northwest', 'FontSize',12)
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');

plotName = [experimentName '_' typeCC '_' queryOptoParameter,'_vs_', queryRunFeature];
plotFullPath = fullfile(savePlotPath, plotName);
% saveas(gcf, plotFullPath, 'pdf');
% saveas(gcf, plotFullPath, 'fig');

traceWithLight = selectedRunGroupsAnalysisTable.trace{42,1};%45

%%

figure;
hold on;
plot(traceNoLight,'Color',[0 0 0],'LineWidth',2)%[0.7 0.7 0.7]
hold on;
plot(traceWithLight,'Color',[0 0.4470 0.7410],'LineWidth',2)

xLimMin = 1.8;
xLimMax = 3.2;
xticks = linspace(xLimMin,xLimMax, 8);
xlim([xLimMin*10000, xLimMax*10000]);
set(gca, 'XTickLabel', num2str((xticks'-xLimMin)))
xlabel('Time [s]');
ylabel('Membrane potential [mV]');
box on;
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');
legend('no Light','with Light','Location', 'northwest', 'FontSize',12)

plotName = [experimentName '_' typeCC '_traces'];
plotFullPath = fullfile(savePlotPath, plotName);
saveas(gcf, plotFullPath, 'pdf');
saveas(gcf, plotFullPath, 'fig');

clear xticks