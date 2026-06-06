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
experimentName = 'Pitx_Jaws_Ephys_9_240326';

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

typeCC = 'noLight';
savePlotPath = [experimentDirectory filesep 'mouseAnalysis'];


%% amplitudeCurrentPulse vs areaVtPulse

queryOptoParameter =  'amplitudeCurrentPulse';
queryRunFeature = 'areaVtPulse';

selectedCellsIDs = [4];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0') & cell2mat(runGroupsParametersTable.widthCurrentPulse) == 6000);
selectedRunGroupsAnalysisTable = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroups),:);

queryOptoParameterValues = zeros(numel(selectedRunGroups),1);
queryOptoParameterArray = [];
queryRunFeatureData = {};
cellNameData = {};

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

figure;

[~,parametersPosition,~] = unique(queryOptoParameterValues);
queryRunFeatureDataArray = reshape(allData,numel(allData)/9,9);
scatter(parametersPosition, queryRunFeatureDataArray,'k','filled');
hold on;

meanQueryRunFeatureDataArray = nanmean(queryRunFeatureDataArray,1);
stdQueryRunFeatureDataArray = nanstd(queryRunFeatureDataArray,1);
xValues = linspace(1,numel(queryOptoParameterValues),numel(queryOptoParameterValues));
bplot = bar(xValues, meanQueryRunFeatureDataArray,'FaceColor',[0 0.4470 0.7410],'BarWidth', 0.5);
errorbar(xValues,meanQueryRunFeatureDataArray,stdQueryRunFeatureDataArray,'.','Color',[0.3 0.3 0.3],'LineWidth', 1.2, 'CapSize', 14)

bplot.ShowBaseLine = 'on';
hScatter = findobj(gca, 'Type', 'scatter');
uistack(hScatter, 'top');

xticks(xValues);
xticklabels([unique(queryOptoParameterValues)]);
xlabel(queryOptoParameter);
ylabel(queryRunFeature);
box on;
title([queryOptoParameter,' vs ', queryRunFeature]);
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');

plotName = [experimentName '_' typeCC '_' queryOptoParameter,'_vs_', queryRunFeature];
plotFullPath = fullfile(savePlotPath, plotName);
saveas(gcf, plotFullPath, 'pdf')
saveas(gcf, plotFullPath, 'fig')

%% amplitudeCurrentPulse vs nAPduringLight

queryOptoParameter =  'amplitudeCurrentPulse';
queryRunFeature = 'nAPtotal';

selectedCellsIDs = [4];
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

figure;

[~,parametersPosition,~] = unique(queryOptoParameterValues);
queryRunFeatureDataArray = reshape(allData,numel(allData)/9,9);
scatter(parametersPosition, queryRunFeatureDataArray,'k','filled');
hold on;

meanQueryRunFeatureDataArray = nanmean(queryRunFeatureDataArray,1);
stdQueryRunFeatureDataArray = nanstd(queryRunFeatureDataArray,1);
xValues = linspace(1,numel(queryOptoParameterValues),numel(queryOptoParameterValues));
bplot = bar(xValues, meanQueryRunFeatureDataArray,'FaceColor',[0 0.4470 0.7410],'BarWidth', 0.5);
errorbar(xValues,meanQueryRunFeatureDataArray,stdQueryRunFeatureDataArray,'.','Color',[0.3 0.3 0.3],'LineWidth', 1.2, 'CapSize', 14)

bplot.ShowBaseLine = 'on';
hScatter = findobj(gca, 'Type', 'scatter');
uistack(hScatter, 'top');

xticks(xValues)
xticklabels([unique(queryOptoParameterValues)])
xlabel(queryOptoParameter);
ylabel(queryRunFeature);
box on;
title([queryOptoParameter,' vs ', queryRunFeature]);
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');

plotName = [experimentName '_' typeCC '_' queryOptoParameter,'_vs_', queryRunFeature];
plotFullPath = fullfile(savePlotPath, plotName);
saveas(gcf, plotFullPath, 'pdf')
saveas(gcf, plotFullPath, 'fig')
