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

typeCC = 'withLight';
savePlotPath = [experimentDirectory filesep 'mouseAnalysis'];

%% widthLightPulse vs rateAPpostLight

queryOptoParameter =  'widthLightPulse';
queryRunFeature = 'rateAPpostLight';

selectedCellsIDs = [4];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao2') & cell2mat(runGroupsParametersTable.widthCurrentPulse) == 12000 & cell2mat(runGroupsParametersTable.amplitudeCurrentPulse) == 500 & cell2mat(runGroupsParametersTable.amplitudeLightPulse) == 5);
selectedRunGroupsAnalysisTable = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroups),:);

controlRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0') & cell2mat(runGroupsParametersTable.widthCurrentPulse) == 12000 & cell2mat(runGroupsParametersTable.amplitudeCurrentPulse) == 500);
controlRunGroupsAnalysisTable = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,controlRunGroups),:);

queryOptoParameterValues = zeros(numel(selectedRunGroups)+numel(controlRunGroups),1);
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

for iRunGroupControl = 1:numel(controlRunGroups)
    
    rowsControlRunGroup = find(controlRunGroupsAnalysisTable.runGroup(:) == controlRunGroups(iRunGroupControl));
    queryOptoParameterValues(groupsCounter) = controlRunGroupsAnalysisTable.optoParameters{rowsControlRunGroup(1)}.(queryOptoParameter);
    
    queryOptoParameterArray = [queryOptoParameterArray; repmat(queryOptoParameterValues(groupsCounter),numel(rowsControlRunGroup),1)];
    queryRunFeatureData{groupsCounter} = cell2mat(controlRunGroupsAnalysisTable.(queryRunFeature)(rowsControlRunGroup));
    cellNameData{groupsCounter} = cell2mat(controlRunGroupsAnalysisTable.cellName(rowsControlRunGroup));

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

%% widthLightPulse vs delayFirstAPpostLight

queryOptoParameter =  'widthLightPulse';
queryRunFeature = 'delayFirstAPpostLight';

selectedCellsIDs = [4, 8];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.widthCurrentPulse) == 12000 & cell2mat(runGroupsParametersTable.amplitudeCurrentPulse) == 500 & cell2mat(runGroupsParametersTable.amplitudeLightPulse) == 5);
selectedRunGroupsAnalysisTable = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroups),:);

controlRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0') & cell2mat(runGroupsParametersTable.widthCurrentPulse) == 12000 & cell2mat(runGroupsParametersTable.amplitudeCurrentPulse) == 500);
controlRunGroupsAnalysisTable = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,controlRunGroups),:);

queryOptoParameterValues = zeros(numel(selectedRunGroups)+numel(controlRunGroups),1);
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

for iRunGroupControl = 1:numel(controlRunGroups)
    
    rowsControlRunGroup = find(controlRunGroupsAnalysisTable.runGroup(:) == controlRunGroups(iRunGroupControl));
    queryOptoParameterValues(groupsCounter) = controlRunGroupsAnalysisTable.optoParameters{rowsControlRunGroup(1)}.(queryOptoParameter);
    
    queryOptoParameterArray = [queryOptoParameterArray; repmat(queryOptoParameterValues(groupsCounter),numel(rowsControlRunGroup),1)];
    queryRunFeatureData{groupsCounter} = cell2mat(controlRunGroupsAnalysisTable.(queryRunFeature)(rowsControlRunGroup));
    cellNameData{groupsCounter} = cell2mat(controlRunGroupsAnalysisTable.cellName(rowsControlRunGroup));

    groupsCounter = groupsCounter + 1;
    
end

allData = vertcat(queryRunFeatureData{:});
cellNameDataArray = vertcat(cellNameData{:});

meansQueryRunFeatureDataAllCells = nan(numel(unique(cellNameDataArray)),numel(unique(queryOptoParameterValues)));

figure;

for iQueryOptoParameter = 1:length(queryRunFeatureData)

    cellNames = unique(cellNameData{iQueryOptoParameter});
    
    for iCell = 1:numel(cellNames)
        
        cellRows = find(cellNameData{iQueryOptoParameter} == cellNames(iCell)); 
        meanQueryRunFeatureDataCell = nanmean(queryRunFeatureData{iQueryOptoParameter}(cellRows));
        positionOptoParameterValue = find(sort(queryOptoParameterValues) == queryOptoParameterValues(iQueryOptoParameter));
        scatter(positionOptoParameterValue, meanQueryRunFeatureDataCell,'k','filled');
        hold on;
        
        meansQueryRunFeatureDataAllCells(iCell,positionOptoParameterValue) = meanQueryRunFeatureDataCell;
        
    end
    
end

meanQueryRunFeatureDataAllCells = nanmean(meansQueryRunFeatureDataAllCells,1);
stdQueryRunFeatureDataAllCells = nanstd(meansQueryRunFeatureDataAllCells,1);
xValues = linspace(1,numel(queryOptoParameterValues),numel(queryOptoParameterValues));
bar(xValues, meanQueryRunFeatureDataAllCells,'FaceColor',[0 0.4470 0.7410],'BarWidth', 0.5);
errorbar(xValues,meanQueryRunFeatureDataAllCells,stdQueryRunFeatureDataAllCells,'.','Color',[0.3 0.3 0.3],'LineWidth', 1.2, 'CapSize', 14)

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
saveas(gcf, plotFullPath, 'pdf');
saveas(gcf, plotFullPath, 'fig');

%% widthLightPulse vs areaVtPulse

queryOptoParameter =  'widthLightPulse';
queryRunFeature = 'areaVtPulse';

selectedCellsIDs = [4, 8];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.widthCurrentPulse) == 12000 & cell2mat(runGroupsParametersTable.amplitudeCurrentPulse) == 500 & cell2mat(runGroupsParametersTable.amplitudeLightPulse) == 5);
selectedRunGroupsAnalysisTable = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroups),:);

controlRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0') & cell2mat(runGroupsParametersTable.widthCurrentPulse) == 12000 & cell2mat(runGroupsParametersTable.amplitudeCurrentPulse) == 500);
controlRunGroupsAnalysisTable = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,controlRunGroups),:);

queryOptoParameterValues = zeros(numel(selectedRunGroups)+numel(controlRunGroups),1);
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

for iRunGroupControl = 1:numel(controlRunGroups)
    
    rowsControlRunGroup = find(controlRunGroupsAnalysisTable.runGroup(:) == controlRunGroups(iRunGroupControl));
    queryOptoParameterValues(groupsCounter) = controlRunGroupsAnalysisTable.optoParameters{rowsControlRunGroup(1)}.(queryOptoParameter);
    
    queryOptoParameterArray = [queryOptoParameterArray; repmat(queryOptoParameterValues(groupsCounter),numel(rowsControlRunGroup),1)];
    queryRunFeatureData{groupsCounter} = cell2mat(controlRunGroupsAnalysisTable.(queryRunFeature)(rowsControlRunGroup));
    cellNameData{groupsCounter} = cell2mat(controlRunGroupsAnalysisTable.cellName(rowsControlRunGroup));

    groupsCounter = groupsCounter + 1;
    
end

allData = vertcat(queryRunFeatureData{:});
cellNameDataArray = vertcat(cellNameData{:});

meansQueryRunFeatureDataAllCells = nan(numel(unique(cellNameDataArray)),numel(unique(queryOptoParameterValues)));

figure;

for iQueryOptoParameter = 1:length(queryRunFeatureData)

    cellNames = unique(cellNameData{iQueryOptoParameter});
    
    for iCell = 1:numel(cellNames)
        
        cellRows = find(cellNameData{iQueryOptoParameter} == cellNames(iCell)); 
        meanQueryRunFeatureDataCell = nanmean(queryRunFeatureData{iQueryOptoParameter}(cellRows));
        positionOptoParameterValue = find(sort(queryOptoParameterValues) == queryOptoParameterValues(iQueryOptoParameter));
        scatter(positionOptoParameterValue, meanQueryRunFeatureDataCell,'k','filled');
        hold on;
        
        meansQueryRunFeatureDataAllCells(iCell,positionOptoParameterValue) = meanQueryRunFeatureDataCell;
        
    end
    
end

meanQueryRunFeatureDataAllCells = nanmean(meansQueryRunFeatureDataAllCells,1);
stdQueryRunFeatureDataAllCells = nanstd(meansQueryRunFeatureDataAllCells,1);
xValues = linspace(1,numel(queryOptoParameterValues),numel(queryOptoParameterValues));
bar(xValues, meanQueryRunFeatureDataAllCells,'FaceColor',[0 0.4470 0.7410],'BarWidth', 0.5);
errorbar(xValues,meanQueryRunFeatureDataAllCells,stdQueryRunFeatureDataAllCells,'.','Color',[0.3 0.3 0.3],'LineWidth', 1.2, 'CapSize', 14)

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
saveas(gcf, plotFullPath, 'pdf');
saveas(gcf, plotFullPath, 'fig');

%% amplitudeLightPulse vs delayFirstAPpostLight

queryOptoParameter =  'amplitudeLightPulse';
queryRunFeature = 'delayFirstAPpostLight';

selectedCellsIDs = [4, 8];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.widthCurrentPulse) == 12000 & cell2mat(runGroupsParametersTable.amplitudeCurrentPulse) == 500 & cell2mat(runGroupsParametersTable.widthLightPulse) == 2500);
selectedRunGroupsAnalysisTable = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroups),:);

controlRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0') & cell2mat(runGroupsParametersTable.widthCurrentPulse) == 12000 & cell2mat(runGroupsParametersTable.amplitudeCurrentPulse) == 500);
controlRunGroupsAnalysisTable = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,controlRunGroups),:);

queryOptoParameterValues = zeros(numel(selectedRunGroups)+numel(controlRunGroups),1);
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

for iRunGroupControl = 1:numel(controlRunGroups)
    
    rowsControlRunGroup = find(controlRunGroupsAnalysisTable.runGroup(:) == controlRunGroups(iRunGroupControl));
    queryOptoParameterValues(groupsCounter) = controlRunGroupsAnalysisTable.optoParameters{rowsControlRunGroup(1)}.(queryOptoParameter);
    
    queryOptoParameterArray = [queryOptoParameterArray; repmat(queryOptoParameterValues(groupsCounter),numel(rowsControlRunGroup),1)];
    queryRunFeatureData{groupsCounter} = cell2mat(controlRunGroupsAnalysisTable.(queryRunFeature)(rowsControlRunGroup));
    cellNameData{groupsCounter} = cell2mat(controlRunGroupsAnalysisTable.cellName(rowsControlRunGroup));

    groupsCounter = groupsCounter + 1;
    
end

allData = vertcat(queryRunFeatureData{:});
cellNameDataArray = vertcat(cellNameData{:});

meansQueryRunFeatureDataAllCells = nan(numel(unique(cellNameDataArray)),numel(unique(queryOptoParameterValues)));

figure;

for iQueryOptoParameter = 1:length(queryRunFeatureData)

    cellNames = unique(cellNameData{iQueryOptoParameter});
    
    for iCell = 1:numel(cellNames)
        
        cellRows = find(cellNameData{iQueryOptoParameter} == cellNames(iCell)); 
        meanQueryRunFeatureDataCell = nanmean(queryRunFeatureData{iQueryOptoParameter}(cellRows));
        positionOptoParameterValue = find(sort(queryOptoParameterValues) == queryOptoParameterValues(iQueryOptoParameter));
        scatter(positionOptoParameterValue, meanQueryRunFeatureDataCell,'k','filled');
        hold on;
        
        meansQueryRunFeatureDataAllCells(iCell,positionOptoParameterValue) = meanQueryRunFeatureDataCell;
        
    end
    
end

meanQueryRunFeatureDataAllCells = nanmean(meansQueryRunFeatureDataAllCells,1);
stdQueryRunFeatureDataAllCells = nanstd(meansQueryRunFeatureDataAllCells,1);
xValues = linspace(1,numel(queryOptoParameterValues),numel(queryOptoParameterValues));
bar(xValues, meanQueryRunFeatureDataAllCells,'FaceColor',[0 0.4470 0.7410],'BarWidth', 0.5);
errorbar(xValues,meanQueryRunFeatureDataAllCells,stdQueryRunFeatureDataAllCells,'.','Color',[0.3 0.3 0.3],'LineWidth', 1.2, 'CapSize', 14)

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
saveas(gcf, plotFullPath, 'pdf');
saveas(gcf, plotFullPath, 'fig');

%% amplitudeLightPulse vs nAPduringLight

analyzedMouseDataFolderName = 'mouseAnalysis2';
analyzedMouseDataFolderPath = [experimentDirectory filesep analyzedMouseDataFolderName];
mouseTableName = 'MouseAnalysisTable.mat';
mouseTableName = strrep(mouseTableName, ' ', '');
mouseParametersTableName = 'MouseParametersTable.mat';
mouseParametersTableName = strrep(mouseParametersTableName, ' ', '');
mouseAnalysisTablePath = fullfile(analyzedMouseDataFolderPath, mouseTableName);
mouseParametersTablePath = fullfile(analyzedMouseDataFolderPath, mouseParametersTableName);

mouseAnalysisTable = load(mouseAnalysisTablePath);
mouseAnalysisTable = mouseAnalysisTable.mouseAnalysisTable;

queryOptoParameter =  'amplitudeLightPulse';
queryRunFeature = 'nAPduringLight';

selectedCellsIDs = [4, 8];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.widthCurrentPulse) == 12000 & cell2mat(runGroupsParametersTable.amplitudeCurrentPulse) == 500 & cell2mat(runGroupsParametersTable.widthLightPulse) == 2500);
selectedRunGroupsAnalysisTable = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroups),:);

controlRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0') & cell2mat(runGroupsParametersTable.widthCurrentPulse) == 12000 & cell2mat(runGroupsParametersTable.amplitudeCurrentPulse) == 500);
controlRunGroupsAnalysisTable = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,controlRunGroups),:);

queryOptoParameterValues = zeros(numel(selectedRunGroups)+numel(controlRunGroups),1);
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

for iRunGroupControl = 1:numel(controlRunGroups)
    
    rowsControlRunGroup = find(controlRunGroupsAnalysisTable.runGroup(:) == controlRunGroups(iRunGroupControl));
    queryOptoParameterValues(groupsCounter) = controlRunGroupsAnalysisTable.optoParameters{rowsControlRunGroup(1)}.(queryOptoParameter);
    
    queryOptoParameterArray = [queryOptoParameterArray; repmat(queryOptoParameterValues(groupsCounter),numel(rowsControlRunGroup),1)];
    queryRunFeatureData{groupsCounter} = cell2mat(controlRunGroupsAnalysisTable.(queryRunFeature)(rowsControlRunGroup));
    cellNameData{groupsCounter} = cell2mat(controlRunGroupsAnalysisTable.cellName(rowsControlRunGroup));

    groupsCounter = groupsCounter + 1;
    
end

allData = vertcat(queryRunFeatureData{:});
cellNameDataArray = vertcat(cellNameData{:});

meansQueryRunFeatureDataAllCells = nan(numel(unique(cellNameDataArray)),numel(unique(queryOptoParameterValues)));

figure;

for iQueryOptoParameter = 1:length(queryRunFeatureData)

    cellNames = unique(cellNameData{iQueryOptoParameter});
    
    for iCell = 1:numel(cellNames)
        
        cellRows = find(cellNameData{iQueryOptoParameter} == cellNames(iCell)); 
        meanQueryRunFeatureDataCell = nanmean(queryRunFeatureData{iQueryOptoParameter}(cellRows));
        positionOptoParameterValue = find(sort(queryOptoParameterValues) == queryOptoParameterValues(iQueryOptoParameter));
        scatter(positionOptoParameterValue, meanQueryRunFeatureDataCell,'k','filled');
        hold on;
        
        meansQueryRunFeatureDataAllCells(iCell,positionOptoParameterValue) = meanQueryRunFeatureDataCell;
        
    end
    
end

meanQueryRunFeatureDataAllCells = nanmean(meansQueryRunFeatureDataAllCells,1);
stdQueryRunFeatureDataAllCells = nanstd(meansQueryRunFeatureDataAllCells,1);
xValues = linspace(1,numel(queryOptoParameterValues),numel(queryOptoParameterValues));
bar(xValues, meanQueryRunFeatureDataAllCells,'FaceColor',[0 0.4470 0.7410],'BarWidth', 0.5);
errorbar(xValues,meanQueryRunFeatureDataAllCells,stdQueryRunFeatureDataAllCells,'.','Color',[0.3 0.3 0.3],'LineWidth', 1.2, 'CapSize', 20)

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
saveas(gcf, plotFullPath, 'pdf');
saveas(gcf, plotFullPath, 'fig');

%% amplitudeCurrentPulse vs areaVtPulse -- withLightBeforeCurrent -- bar plot

typeCC = 'withLightBeforeCurrent';
queryOptoParameter =  'amplitudeCurrentPulse';
queryRunFeature = 'areaVtPulse';

selectedCellsIDs = [4, 8];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.widthCurrentPulse) == 6000 & cell2mat(runGroupsParametersTable.delayLightPulse) == 19000);
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

figure;

for iQueryOptoParameter = 1:length(queryRunFeatureData)

    cellNames = unique(cellNameData{iQueryOptoParameter});
    
    for iCell = 1:numel(cellNames)
        
        cellRows = find(cellNameData{iQueryOptoParameter} == cellNames(iCell)); 
        meanQueryRunFeatureDataCell = nanmean(queryRunFeatureData{iQueryOptoParameter}(cellRows));
        positionOptoParameterValue = find(sort(queryOptoParameterValues) == queryOptoParameterValues(iQueryOptoParameter));
        scatter(positionOptoParameterValue, meanQueryRunFeatureDataCell,'k','filled');
        hold on;
        
        meansQueryRunFeatureDataAllCells(iCell,positionOptoParameterValue) = meanQueryRunFeatureDataCell;
        
    end
    
end

meanQueryRunFeatureDataAllCells = nanmean(meansQueryRunFeatureDataAllCells,1);
stdQueryRunFeatureDataAllCells = nanstd(meansQueryRunFeatureDataAllCells,1);
xValues = linspace(1,numel(queryOptoParameterValues),numel(queryOptoParameterValues));
bar(xValues, meanQueryRunFeatureDataAllCells,'FaceColor',[0 0.4470 0.7410],'BarWidth', 0.5);
errorbar(xValues,meanQueryRunFeatureDataAllCells,stdQueryRunFeatureDataAllCells,'.','Color',[0.3 0.3 0.3],'LineWidth', 1.2, 'CapSize', 14)

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
saveas(gcf, plotFullPath, 'pdf');
saveas(gcf, plotFullPath, 'fig');

%% With light

selectedCellsIDs = [4];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao2') & cell2mat(runGroupsParametersTable.widthCurrentPulse) == 12000 & cell2mat(runGroupsParametersTable.widthLightPulse) == 1000);
selectedRunGroupsAnalysisTable = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroups),:);

traceWithLight = selectedRunGroupsAnalysisTable.trace{3,1};

figure;
plot(traceWithLight,'Color',[0 0.4470 0.7410],'LineWidth',1)

xLimMin = 2.0;
xLimMax = 3.8;
xticks = linspace(xLimMin,xLimMax, 10);
xlim([xLimMin*10000, xLimMax*10000]);
set(gca, 'XTickLabel', num2str((xticks'-xLimMin)))
xlabel('Time [s]');
ylabel('Membrane potential [mV]');
box on;
legend('with Light','Location', 'northwest', 'FontSize',12)
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');

plotName = [experimentName '_' typeCC '_trace'];
plotFullPath = fullfile(savePlotPath, plotName);
saveas(gcf, plotFullPath, 'pdf');
saveas(gcf, plotFullPath, 'fig');
