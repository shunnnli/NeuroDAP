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
experimentName = 'Pitx_GtACR_Ephys_1';

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

queryOptoParameter =  'amplitudeLightPulse';
queryRunFeature = 'delayFirstAPpostLight';

queryOptoParameterValuesOverall = [NaN, 0, 0.1, 0.3, 0.5, 1, 3, 5];
nNaNs = numel(find(isnan(queryOptoParameterValuesOverall)));

selectedCellsIDs = [9,11];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

% cell 3: 0.5-0.1 current 200
% cell 2:
% cell 5: 0.5-0.1 current 250

selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.widthCurrentPulse) == 12000 & cell2mat(runGroupsParametersTable.widthLightPulse) == 5000);
selectedRunGroupsAnalysisTable = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroups),:);

controlRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0') & cell2mat(runGroupsParametersTable.widthCurrentPulse) == 12000); %& cell2mat(runGroupsParametersTable.amplitudeCurrentPulse) == 250
controlRunGroupsAnalysisTable = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,controlRunGroups),:);

queryOptoParameterValues = nan(numel(selectedRunGroups)+numel(controlRunGroups),1);
queryRunFeatureData = {};
cellNameData = {};

groupsCounter = 1;

for iRunGroup = 1:numel(selectedRunGroups)
    
    rowsRunGroup = find(selectedRunGroupsAnalysisTable.runGroup(:) == selectedRunGroups(iRunGroup));
    
    if ~isempty(rowsRunGroup)
    
        queryOptoParameterValues(groupsCounter) = selectedRunGroupsAnalysisTable.optoParameters{rowsRunGroup(1)}.(queryOptoParameter);
        queryRunFeatureData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTable.(queryRunFeature)(rowsRunGroup));
        cellNameData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTable.cellName(rowsRunGroup))

    end
    
    groupsCounter = groupsCounter + 1;
    
end

for iRunGroupControl = 1:numel(controlRunGroups)
    
    rowsControlRunGroup = find(controlRunGroupsAnalysisTable.runGroup(:) == controlRunGroups(iRunGroupControl));
    
    if ~isempty(rowsControlRunGroup)
        
        queryOptoParameterValues(groupsCounter) = controlRunGroupsAnalysisTable.optoParameters{rowsControlRunGroup(1)}.(queryOptoParameter);
        queryRunFeatureData{groupsCounter} = cell2mat(controlRunGroupsAnalysisTable.(queryRunFeature)(rowsControlRunGroup));
        cellNameData{groupsCounter} = cell2mat(controlRunGroupsAnalysisTable.cellName(rowsControlRunGroup));

    end
    
    groupsCounter = groupsCounter + 1;
    
end

allData = vertcat(queryRunFeatureData{:});
cellNameDataArray = vertcat(cellNameData{:});

meansQueryRunFeatureDataAllCells = nan(numel(unique(cellNameDataArray)),numel(queryOptoParameterValuesOverall));

figure;

for iQueryOptoParameter = 1:length(queryRunFeatureData)

    cellNames = unique(cellNameData{iQueryOptoParameter})
    
    for iCell = 1:numel(cellNames)
        
        cellRows = find(cellNameData{iQueryOptoParameter} == cellNames(iCell)); 
        meanQueryRunFeatureDataCell = nanmean(queryRunFeatureData{iQueryOptoParameter}(cellRows));
        
        if ~isnan(queryOptoParameterValues(iQueryOptoParameter))
        
            positionOptoParameterValue = find(queryOptoParameterValuesOverall == queryOptoParameterValues(iQueryOptoParameter));
        
        else
            
            positionOptoParameterValue = find(isnan(queryOptoParameterValuesOverall));
            
        end

        scatter(positionOptoParameterValue-nNaNs, meanQueryRunFeatureDataCell,'k','filled');
        hold on;
        
        meansQueryRunFeatureDataAllCells(find(selectedCellsIDs == cellNames(iCell)),positionOptoParameterValue) = meanQueryRunFeatureDataCell;
        
    end
    
end

meansQueryRunFeatureDataAllCells1 = meansQueryRunFeatureDataAllCells;

% Indicate main directory and add it to path
matlabDirectory = ['/Users/' user '/Desktop/Thesis/Matlab/'];
mainDirectory = ['/Users/' user '/Desktop/Thesis/Matlab/Experiments/']; % including all Matlab toolboxes and necessary paths
addpath(genpath(matlabDirectory));

% Indicate experiment to analyze
experimentName = 'Pitx_GtACR_Ephys_2';

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

queryOptoParameter =  'amplitudeLightPulse';
queryRunFeature = 'delayFirstAPpostLight';

queryOptoParameterValuesOverall = [NaN, 0, 0.1, 0.3, 0.5, 1, 3, 5];
nNaNs = numel(find(isnan(queryOptoParameterValuesOverall)));

selectedCellsIDs = [4,5];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

% cell 3: 0.5-0.1 current 200
% cell 2:
% cell 5: 0.5-0.1 current 250

selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.widthCurrentPulse) == 12000 & cell2mat(runGroupsParametersTable.widthLightPulse) == 5000);
selectedRunGroupsAnalysisTable = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroups),:);

controlRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0') & cell2mat(runGroupsParametersTable.widthCurrentPulse) == 12000); %& cell2mat(runGroupsParametersTable.amplitudeCurrentPulse) == 250
controlRunGroupsAnalysisTable = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,controlRunGroups),:);

queryOptoParameterValues = nan(numel(selectedRunGroups)+numel(controlRunGroups),1);
queryRunFeatureData = {};
cellNameData = {};

groupsCounter = 1;

for iRunGroup = 1:numel(selectedRunGroups)
    
    rowsRunGroup = find(selectedRunGroupsAnalysisTable.runGroup(:) == selectedRunGroups(iRunGroup));
    
    if ~isempty(rowsRunGroup)
    
        queryOptoParameterValues(groupsCounter) = selectedRunGroupsAnalysisTable.optoParameters{rowsRunGroup(1)}.(queryOptoParameter);
        queryRunFeatureData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTable.(queryRunFeature)(rowsRunGroup));
        cellNameData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTable.cellName(rowsRunGroup))

    end
    
    groupsCounter = groupsCounter + 1;
    
end

for iRunGroupControl = 1:numel(controlRunGroups)
    
    rowsControlRunGroup = find(controlRunGroupsAnalysisTable.runGroup(:) == controlRunGroups(iRunGroupControl));
    
    if ~isempty(rowsControlRunGroup)
        
        queryOptoParameterValues(groupsCounter) = controlRunGroupsAnalysisTable.optoParameters{rowsControlRunGroup(1)}.(queryOptoParameter);
        queryRunFeatureData{groupsCounter} = cell2mat(controlRunGroupsAnalysisTable.(queryRunFeature)(rowsControlRunGroup));
        cellNameData{groupsCounter} = cell2mat(controlRunGroupsAnalysisTable.cellName(rowsControlRunGroup));

    end
    
    groupsCounter = groupsCounter + 1;
    
end

allData = vertcat(queryRunFeatureData{:});
cellNameDataArray = vertcat(cellNameData{:});

meansQueryRunFeatureDataAllCells = nan(numel(unique(cellNameDataArray)),numel(queryOptoParameterValuesOverall));


for iQueryOptoParameter = 1:length(queryRunFeatureData)

    cellNames = unique(cellNameData{iQueryOptoParameter})
    
    for iCell = 1:numel(cellNames)
        
        cellRows = find(cellNameData{iQueryOptoParameter} == cellNames(iCell)); 
        meanQueryRunFeatureDataCell = nanmean(queryRunFeatureData{iQueryOptoParameter}(cellRows));
        
        if ~isnan(queryOptoParameterValues(iQueryOptoParameter))
        
            positionOptoParameterValue = find(queryOptoParameterValuesOverall == queryOptoParameterValues(iQueryOptoParameter));
        
        else
            
            positionOptoParameterValue = find(isnan(queryOptoParameterValuesOverall));
            
        end

        scatter(positionOptoParameterValue-nNaNs, meanQueryRunFeatureDataCell,'k','filled');
        hold on;
        
        meansQueryRunFeatureDataAllCells(find(selectedCellsIDs == cellNames(iCell)),positionOptoParameterValue) = meanQueryRunFeatureDataCell;
        
    end
    
end

meansQueryRunFeatureDataAllCells = [meansQueryRunFeatureDataAllCells1; meansQueryRunFeatureDataAllCells]

meanQueryRunFeatureDataAllCells = nanmean(meansQueryRunFeatureDataAllCells,1);
stdQueryRunFeatureDataAllCells = nanstd(meansQueryRunFeatureDataAllCells,1);
xValues = linspace(1,numel(queryOptoParameterValuesOverall)-nNaNs,numel(queryOptoParameterValuesOverall)-nNaNs);
plot(xValues, meanQueryRunFeatureDataAllCells(2:end),'Color',[0 0.4470 0.7410],'LineWidth', 2);
errorbar(xValues,meanQueryRunFeatureDataAllCells(2:end),stdQueryRunFeatureDataAllCells(2:end),'.','Color',[0 0 0],'LineWidth', 1.2, 'CapSize', 14)

hScatter = findobj(gca, 'Type', 'scatter');
uistack(hScatter, 'top');

xticks(xValues)
xticklabels(queryOptoParameterValuesOverall(~isnan(queryOptoParameterValuesOverall)))
xlabel(queryOptoParameter);
ylabel(queryRunFeature);
box on;
title([queryOptoParameter,' vs ', queryRunFeature]);
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');

plotName = [experimentName '_' typeCC '_' queryOptoParameter,'_vs_', queryRunFeature, '_2mice'];
plotFullPath = fullfile(savePlotPath, plotName);
% saveas(gcf, plotFullPath, 'pdf');
% saveas(gcf, plotFullPath, 'fig');

%%

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
experimentName = 'Pitx_GtACR_Ephys_1';

% Set path
experimentDirectory = [mainDirectory experimentName];
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
runGroupsParametersTable = load(mouseParametersTablePath);
runGroupsParametersTable = runGroupsParametersTable.runGroupsParametersTable;

typeCC = 'withLight';
savePlotPath = [experimentDirectory filesep 'mouseAnalysis'];

queryOptoParameter =  'amplitudeLightPulse';
queryRunFeature = 'nAPduringLight';

queryOptoParameterValuesOverall = [NaN, 0, 0.1, 0.3, 0.5, 1, 3, 5];
nNaNs = numel(find(isnan(queryOptoParameterValuesOverall)));

selectedCellsIDs = [11];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

% cell 3: 0.5-0.1 current 200
% cell 2:
% cell 5: 0.5-0.1 current 250

selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.widthCurrentPulse) == 12000 & cell2mat(runGroupsParametersTable.widthLightPulse) == 5000);
selectedRunGroupsAnalysisTable = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroups),:);

controlRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0') & cell2mat(runGroupsParametersTable.widthCurrentPulse) == 12000); %& cell2mat(runGroupsParametersTable.amplitudeCurrentPulse) == 250
controlRunGroupsAnalysisTable = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,controlRunGroups),:);

queryOptoParameterValues = nan(numel(selectedRunGroups)+numel(controlRunGroups),1);
queryRunFeatureData = {};
cellNameData = {};

groupsCounter = 1;

for iRunGroup = 1:numel(selectedRunGroups)
    
    rowsRunGroup = find(selectedRunGroupsAnalysisTable.runGroup(:) == selectedRunGroups(iRunGroup));
    
    if ~isempty(rowsRunGroup)
    
        queryOptoParameterValues(groupsCounter) = selectedRunGroupsAnalysisTable.optoParameters{rowsRunGroup(1)}.(queryOptoParameter);
        queryRunFeatureData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTable.(queryRunFeature)(rowsRunGroup));
        cellNameData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTable.cellName(rowsRunGroup))

    end
    
    groupsCounter = groupsCounter + 1;
    
end

for iRunGroupControl = 1:numel(controlRunGroups)
    
    rowsControlRunGroup = find(controlRunGroupsAnalysisTable.runGroup(:) == controlRunGroups(iRunGroupControl));
    
    if ~isempty(rowsControlRunGroup)
        
        queryOptoParameterValues(groupsCounter) = controlRunGroupsAnalysisTable.optoParameters{rowsControlRunGroup(1)}.(queryOptoParameter);
        queryRunFeatureData{groupsCounter} = cell2mat(controlRunGroupsAnalysisTable.(queryRunFeature)(rowsControlRunGroup));
        cellNameData{groupsCounter} = cell2mat(controlRunGroupsAnalysisTable.cellName(rowsControlRunGroup));

    end
    
    groupsCounter = groupsCounter + 1;
    
end

allData = vertcat(queryRunFeatureData{:});
cellNameDataArray = vertcat(cellNameData{:});

meansQueryRunFeatureDataAllCells = nan(numel(unique(cellNameDataArray)),numel(queryOptoParameterValuesOverall));

figure;

for iQueryOptoParameter = 1:length(queryRunFeatureData)

    cellNames = unique(cellNameData{iQueryOptoParameter})
    
    for iCell = 1:numel(cellNames)
        
        cellRows = find(cellNameData{iQueryOptoParameter} == cellNames(iCell)); 
        meanQueryRunFeatureDataCell = nanmean(queryRunFeatureData{iQueryOptoParameter}(cellRows));
        
        if ~isnan(queryOptoParameterValues(iQueryOptoParameter))
        
            positionOptoParameterValue = find(queryOptoParameterValuesOverall == queryOptoParameterValues(iQueryOptoParameter));
        
        else
            
            positionOptoParameterValue = find(isnan(queryOptoParameterValuesOverall));
            
        end

        scatter(positionOptoParameterValue-nNaNs, meanQueryRunFeatureDataCell,'k','filled');
        hold on;
        
        meansQueryRunFeatureDataAllCells(find(selectedCellsIDs == cellNames(iCell)),positionOptoParameterValue) = meanQueryRunFeatureDataCell;
        
    end
    
end

meansQueryRunFeatureDataAllCells1 = meansQueryRunFeatureDataAllCells;

% Indicate main directory and add it to path
matlabDirectory = ['/Users/' user '/Desktop/Thesis/Matlab/'];
mainDirectory = ['/Users/' user '/Desktop/Thesis/Matlab/Experiments/']; % including all Matlab toolboxes and necessary paths
addpath(genpath(matlabDirectory));

% Indicate experiment to analyze
experimentName = 'Pitx_GtACR_Ephys_2';

% Set path
experimentDirectory = [mainDirectory experimentName];
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
runGroupsParametersTable = load(mouseParametersTablePath);
runGroupsParametersTable = runGroupsParametersTable.runGroupsParametersTable;

typeCC = 'withLight';
savePlotPath = [experimentDirectory filesep 'mouseAnalysis'];

queryOptoParameter =  'amplitudeLightPulse';
queryRunFeature = 'nAPduringLight';

queryOptoParameterValuesOverall = [NaN, 0, 0.1, 0.3, 0.5, 1, 3, 5];
nNaNs = numel(find(isnan(queryOptoParameterValuesOverall)));

selectedCellsIDs = [4,5];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

% cell 3: 0.5-0.1 current 200
% cell 2:
% cell 5: 0.5-0.1 current 250

selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.widthCurrentPulse) == 12000 & cell2mat(runGroupsParametersTable.widthLightPulse) == 5000);
selectedRunGroupsAnalysisTable = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroups),:);

controlRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0') & cell2mat(runGroupsParametersTable.widthCurrentPulse) == 12000); %& cell2mat(runGroupsParametersTable.amplitudeCurrentPulse) == 250
controlRunGroupsAnalysisTable = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,controlRunGroups),:);

queryOptoParameterValues = nan(numel(selectedRunGroups)+numel(controlRunGroups),1);
queryRunFeatureData = {};
cellNameData = {};

groupsCounter = 1;

for iRunGroup = 1:numel(selectedRunGroups)
    
    rowsRunGroup = find(selectedRunGroupsAnalysisTable.runGroup(:) == selectedRunGroups(iRunGroup));
    
    if ~isempty(rowsRunGroup)
    
        queryOptoParameterValues(groupsCounter) = selectedRunGroupsAnalysisTable.optoParameters{rowsRunGroup(1)}.(queryOptoParameter);
        queryRunFeatureData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTable.(queryRunFeature)(rowsRunGroup));
        cellNameData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTable.cellName(rowsRunGroup))

    end
    
    groupsCounter = groupsCounter + 1;
    
end

for iRunGroupControl = 1:numel(controlRunGroups)
    
    rowsControlRunGroup = find(controlRunGroupsAnalysisTable.runGroup(:) == controlRunGroups(iRunGroupControl));
    
    if ~isempty(rowsControlRunGroup)
        
        queryOptoParameterValues(groupsCounter) = controlRunGroupsAnalysisTable.optoParameters{rowsControlRunGroup(1)}.(queryOptoParameter);
        queryRunFeatureData{groupsCounter} = cell2mat(controlRunGroupsAnalysisTable.(queryRunFeature)(rowsControlRunGroup));
        cellNameData{groupsCounter} = cell2mat(controlRunGroupsAnalysisTable.cellName(rowsControlRunGroup));

    end
    
    groupsCounter = groupsCounter + 1;
    
end

allData = vertcat(queryRunFeatureData{:});
cellNameDataArray = vertcat(cellNameData{:});

meansQueryRunFeatureDataAllCells = nan(numel(unique(cellNameDataArray)),numel(queryOptoParameterValuesOverall));


for iQueryOptoParameter = 1:length(queryRunFeatureData)

    cellNames = unique(cellNameData{iQueryOptoParameter})
    
    for iCell = 1:numel(cellNames)
        
        cellRows = find(cellNameData{iQueryOptoParameter} == cellNames(iCell)); 
        meanQueryRunFeatureDataCell = nanmean(queryRunFeatureData{iQueryOptoParameter}(cellRows));
        
        if ~isnan(queryOptoParameterValues(iQueryOptoParameter))
        
            positionOptoParameterValue = find(queryOptoParameterValuesOverall == queryOptoParameterValues(iQueryOptoParameter));
        
        else
            
            positionOptoParameterValue = find(isnan(queryOptoParameterValuesOverall));
            
        end

        scatter(positionOptoParameterValue-nNaNs, meanQueryRunFeatureDataCell,'k','filled');
        hold on;
        
        meansQueryRunFeatureDataAllCells(find(selectedCellsIDs == cellNames(iCell)),positionOptoParameterValue) = meanQueryRunFeatureDataCell;
        
    end
    
end

meansQueryRunFeatureDataAllCells = [meansQueryRunFeatureDataAllCells1; meansQueryRunFeatureDataAllCells]

meanQueryRunFeatureDataAllCells = nanmean(meansQueryRunFeatureDataAllCells,1);
stdQueryRunFeatureDataAllCells = nanstd(meansQueryRunFeatureDataAllCells,1);
xValues = linspace(1,numel(queryOptoParameterValuesOverall)-nNaNs,numel(queryOptoParameterValuesOverall)-nNaNs);
plot(xValues, meanQueryRunFeatureDataAllCells(2:end),'Color',[0 0.4470 0.7410],'LineWidth', 2);
errorbar(xValues,meanQueryRunFeatureDataAllCells(2:end),stdQueryRunFeatureDataAllCells(2:end),'.','Color',[0 0 0],'LineWidth', 1.2, 'CapSize', 14)

hScatter = findobj(gca, 'Type', 'scatter');
uistack(hScatter, 'top');

xticks(xValues)
xticklabels(queryOptoParameterValuesOverall(~isnan(queryOptoParameterValuesOverall)))
xlabel(queryOptoParameter);
ylabel(queryRunFeature);
box on;
title([queryOptoParameter,' vs ', queryRunFeature]);
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');

plotName = [experimentName '_' typeCC '_' queryOptoParameter,'_vs_', queryRunFeature, '_2mice'];
plotFullPath = fullfile(savePlotPath, plotName);
saveas(gcf, plotFullPath, 'pdf');
saveas(gcf, plotFullPath, 'fig');

%%

%%

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
experimentName = 'Pitx_GtACR_Ephys_1';

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
runGroupsParametersTable = load(mouseParametersTablePath);
runGroupsParametersTable = runGroupsParametersTable.runGroupsParametersTable;

typeCC = 'withLight';
savePlotPath = [experimentDirectory filesep 'mouseAnalysis'];

queryOptoParameter =  'widthLightPulse';
queryRunFeature = 'areaVtPulse';

queryOptoParameterValuesOverall = [NaN, 0, 500, 1000, 2500, 5000];
nNaNs = numel(find(isnan(queryOptoParameterValuesOverall)));

selectedCellsIDs = [11];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

% cell 3: 0.5-0.1 current 200
% cell 2:
% cell 5: 0.5-0.1 current 250

selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.widthCurrentPulse) == 1200 & cell2mat(runGroupsParametersTable.amplitudeLightPulse) == 0.5);
selectedRunGroupsAnalysisTable = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroups),:);

controlRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0') & cell2mat(runGroupsParametersTable.widthCurrentPulse) == 1200); %& cell2mat(runGroupsParametersTable.amplitudeCurrentPulse) == 250
controlRunGroupsAnalysisTable = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,controlRunGroups),:);

queryOptoParameterValues = nan(numel(selectedRunGroups)+numel(controlRunGroups),1);
queryRunFeatureData = {};
cellNameData = {};

groupsCounter = 1;

for iRunGroup = 1:numel(selectedRunGroups)
    
    rowsRunGroup = find(selectedRunGroupsAnalysisTable.runGroup(:) == selectedRunGroups(iRunGroup));
    
    if ~isempty(rowsRunGroup)
    
        queryOptoParameterValues(groupsCounter) = selectedRunGroupsAnalysisTable.optoParameters{rowsRunGroup(1)}.(queryOptoParameter);
        queryRunFeatureData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTable.(queryRunFeature)(rowsRunGroup));
        cellNameData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTable.cellName(rowsRunGroup))

    end
    
    groupsCounter = groupsCounter + 1;
    
end

for iRunGroupControl = 1:numel(controlRunGroups)
    
    rowsControlRunGroup = find(controlRunGroupsAnalysisTable.runGroup(:) == controlRunGroups(iRunGroupControl));
    
    if ~isempty(rowsControlRunGroup)
        
        queryOptoParameterValues(groupsCounter) = controlRunGroupsAnalysisTable.optoParameters{rowsControlRunGroup(1)}.(queryOptoParameter);
        queryRunFeatureData{groupsCounter} = cell2mat(controlRunGroupsAnalysisTable.(queryRunFeature)(rowsControlRunGroup));
        cellNameData{groupsCounter} = cell2mat(controlRunGroupsAnalysisTable.cellName(rowsControlRunGroup));

    end
    
    groupsCounter = groupsCounter + 1;
    
end

allData = vertcat(queryRunFeatureData{:});
cellNameDataArray = vertcat(cellNameData{:});

meansQueryRunFeatureDataAllCells = nan(numel(unique(cellNameDataArray)),numel(queryOptoParameterValuesOverall));

figure;

for iQueryOptoParameter = 1:length(queryRunFeatureData)

    cellNames = unique(cellNameData{iQueryOptoParameter})
    
    for iCell = 1:numel(cellNames)
        
        cellRows = find(cellNameData{iQueryOptoParameter} == cellNames(iCell)); 
        meanQueryRunFeatureDataCell = nanmean(queryRunFeatureData{iQueryOptoParameter}(cellRows));
        
        if ~isnan(queryOptoParameterValues(iQueryOptoParameter))
        
            positionOptoParameterValue = find(queryOptoParameterValuesOverall == queryOptoParameterValues(iQueryOptoParameter));
        
        else
            
            positionOptoParameterValue = find(isnan(queryOptoParameterValuesOverall));
            
        end

        scatter(positionOptoParameterValue-nNaNs, meanQueryRunFeatureDataCell,'k','filled');
        hold on;
        
        meansQueryRunFeatureDataAllCells(find(selectedCellsIDs == cellNames(iCell)),positionOptoParameterValue) = meanQueryRunFeatureDataCell;
        
    end
    
end

meansQueryRunFeatureDataAllCells1 = meansQueryRunFeatureDataAllCells;

% Indicate main directory and add it to path
matlabDirectory = ['/Users/' user '/Desktop/Thesis/Matlab/'];
mainDirectory = ['/Users/' user '/Desktop/Thesis/Matlab/Experiments/']; % including all Matlab toolboxes and necessary paths
addpath(genpath(matlabDirectory));

% Indicate experiment to analyze
experimentName = 'Pitx_GtACR_Ephys_2';

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
runGroupsParametersTable = load(mouseParametersTablePath);
runGroupsParametersTable = runGroupsParametersTable.runGroupsParametersTable;

typeCC = 'withLight';
savePlotPath = [experimentDirectory filesep 'mouseAnalysis'];

queryOptoParameter =  'widthLightPulse';
queryRunFeature = 'areaVtPulse';

queryOptoParameterValuesOverall = [NaN, 0, 500, 1000, 2500, 5000];
nNaNs = numel(find(isnan(queryOptoParameterValuesOverall)));

selectedCellsIDs = [4,5];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

% cell 3: 0.5-0.1 current 200
% cell 2:
% cell 5: 0.5-0.1 current 250

selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.widthCurrentPulse) == 12000 & cell2mat(runGroupsParametersTable.amplitudeLightPulse) == 0.5);
selectedRunGroupsAnalysisTable = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroups),:);

controlRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0') & cell2mat(runGroupsParametersTable.widthCurrentPulse) == 12000); %& cell2mat(runGroupsParametersTable.amplitudeCurrentPulse) == 250
controlRunGroupsAnalysisTable = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,controlRunGroups),:);

queryOptoParameterValues = nan(numel(selectedRunGroups)+numel(controlRunGroups),1);
queryRunFeatureData = {};
cellNameData = {};

groupsCounter = 1;

for iRunGroup = 1:numel(selectedRunGroups)
    
    rowsRunGroup = find(selectedRunGroupsAnalysisTable.runGroup(:) == selectedRunGroups(iRunGroup));
    
    if ~isempty(rowsRunGroup)
    
        queryOptoParameterValues(groupsCounter) = selectedRunGroupsAnalysisTable.optoParameters{rowsRunGroup(1)}.(queryOptoParameter);
        queryRunFeatureData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTable.(queryRunFeature)(rowsRunGroup));
        cellNameData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTable.cellName(rowsRunGroup))

    end
    
    groupsCounter = groupsCounter + 1;
    
end

for iRunGroupControl = 1:numel(controlRunGroups)
    
    rowsControlRunGroup = find(controlRunGroupsAnalysisTable.runGroup(:) == controlRunGroups(iRunGroupControl));
    
    if ~isempty(rowsControlRunGroup)
        
        queryOptoParameterValues(groupsCounter) = controlRunGroupsAnalysisTable.optoParameters{rowsControlRunGroup(1)}.(queryOptoParameter);
        queryRunFeatureData{groupsCounter} = cell2mat(controlRunGroupsAnalysisTable.(queryRunFeature)(rowsControlRunGroup));
        cellNameData{groupsCounter} = cell2mat(controlRunGroupsAnalysisTable.cellName(rowsControlRunGroup));

    end
    
    groupsCounter = groupsCounter + 1;
    
end

allData = vertcat(queryRunFeatureData{:});
cellNameDataArray = vertcat(cellNameData{:});

meansQueryRunFeatureDataAllCells = nan(numel(unique(cellNameDataArray)),numel(queryOptoParameterValuesOverall));


for iQueryOptoParameter = 1:length(queryRunFeatureData)

    cellNames = unique(cellNameData{iQueryOptoParameter})
    
    for iCell = 1:numel(cellNames)
        
        cellRows = find(cellNameData{iQueryOptoParameter} == cellNames(iCell)); 
        meanQueryRunFeatureDataCell = nanmean(queryRunFeatureData{iQueryOptoParameter}(cellRows));
        
        if ~isnan(queryOptoParameterValues(iQueryOptoParameter))
        
            positionOptoParameterValue = find(queryOptoParameterValuesOverall == queryOptoParameterValues(iQueryOptoParameter));
        
        else
            
            positionOptoParameterValue = find(isnan(queryOptoParameterValuesOverall));
            
        end

        scatter(positionOptoParameterValue-nNaNs, meanQueryRunFeatureDataCell,'k','filled');
        hold on;
        
        meansQueryRunFeatureDataAllCells(find(selectedCellsIDs == cellNames(iCell)),positionOptoParameterValue) = meanQueryRunFeatureDataCell;
        
    end
    
end

meansQueryRunFeatureDataAllCells = [meansQueryRunFeatureDataAllCells1; meansQueryRunFeatureDataAllCells]

meanQueryRunFeatureDataAllCells = nanmean(meansQueryRunFeatureDataAllCells,1);
stdQueryRunFeatureDataAllCells = nanstd(meansQueryRunFeatureDataAllCells,1);
xValues = linspace(1,numel(queryOptoParameterValuesOverall)-nNaNs,numel(queryOptoParameterValuesOverall)-nNaNs);
bar(xValues, meanQueryRunFeatureDataAllCells(2:end),'FaceColor',[0 0.4470 0.7410],'BarWidth', 0.5);
errorbar(xValues,meanQueryRunFeatureDataAllCells(2:end),stdQueryRunFeatureDataAllCells(2:end),'.','Color',[0 0 0],'LineWidth', 1.2, 'CapSize', 14)

hScatter = findobj(gca, 'Type', 'scatter');
uistack(hScatter, 'top');

xticks(xValues)
xticklabels(queryOptoParameterValuesOverall(~isnan(queryOptoParameterValuesOverall)))
xlabel(queryOptoParameter);
ylabel(queryRunFeature);
box on;
title([queryOptoParameter,' vs ', queryRunFeature]);
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');

plotName = [experimentName '_' typeCC '_' queryOptoParameter,'_vs_', queryRunFeature, '_2mice'];
plotFullPath = fullfile(savePlotPath, plotName);
saveas(gcf, plotFullPath, 'pdf');
saveas(gcf, plotFullPath, 'fig');

%% New

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
experimentName = 'Pitx_GtACR_Ephys_1';

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
runGroupsParametersTable = load(mouseParametersTablePath);
runGroupsParametersTable = runGroupsParametersTable.runGroupsParametersTable;

typeCC = 'withLight';
savePlotPath = [experimentDirectory filesep 'mouseAnalysis'];

queryOptoParameter =  'amplitudeCurrentPulse';
queryRunFeature = 'nAPduringLight';

queryOptoParameterValuesOverall = [NaN, 0, 500, 1000, 2500, 5000];
nNaNs = numel(find(isnan(queryOptoParameterValuesOverall)));

selectedCellsIDs = [11];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

% cell 3: 0.5-0.1 current 200
% cell 2:
% cell 5: 0.5-0.1 current 250

selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.widthCurrentPulse) == 12000 & cell2mat(runGroupsParametersTable.amplitudeLightPulse) == 0.5);
selectedRunGroupsAnalysisTable = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroups),:);

controlRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0') & cell2mat(runGroupsParametersTable.widthCurrentPulse) == 12000); %& cell2mat(runGroupsParametersTable.amplitudeCurrentPulse) == 250
controlRunGroupsAnalysisTable = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,controlRunGroups),:);

queryOptoParameterValues = nan(numel(selectedRunGroups)+numel(controlRunGroups),1);
queryRunFeatureData = {};
cellNameData = {};

groupsCounter = 1;

for iRunGroup = 1:numel(selectedRunGroups)
    
    rowsRunGroup = find(selectedRunGroupsAnalysisTable.runGroup(:) == selectedRunGroups(iRunGroup));
    
    if ~isempty(rowsRunGroup)
    
        queryOptoParameterValues(groupsCounter) = selectedRunGroupsAnalysisTable.optoParameters{rowsRunGroup(1)}.(queryOptoParameter);
        queryRunFeatureData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTable.(queryRunFeature)(rowsRunGroup));
        cellNameData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTable.cellName(rowsRunGroup))

    end
    
    groupsCounter = groupsCounter + 1;
    
end

for iRunGroupControl = 1:numel(controlRunGroups)
    
    rowsControlRunGroup = find(controlRunGroupsAnalysisTable.runGroup(:) == controlRunGroups(iRunGroupControl));
    
    if ~isempty(rowsControlRunGroup)
        
        queryOptoParameterValues(groupsCounter) = controlRunGroupsAnalysisTable.optoParameters{rowsControlRunGroup(1)}.(queryOptoParameter);
        queryRunFeatureData{groupsCounter} = cell2mat(controlRunGroupsAnalysisTable.(queryRunFeature)(rowsControlRunGroup));
        cellNameData{groupsCounter} = cell2mat(controlRunGroupsAnalysisTable.cellName(rowsControlRunGroup));

    end
    
    groupsCounter = groupsCounter + 1;
    
end

allData = vertcat(queryRunFeatureData{:});
cellNameDataArray = vertcat(cellNameData{:});

meansQueryRunFeatureDataAllCells = nan(numel(unique(cellNameDataArray)),numel(queryOptoParameterValuesOverall));

figure;

for iQueryOptoParameter = 1:length(queryRunFeatureData)

    cellNames = unique(cellNameData{iQueryOptoParameter})
    
    for iCell = 1:numel(cellNames)
        
        cellRows = find(cellNameData{iQueryOptoParameter} == cellNames(iCell)); 
        meanQueryRunFeatureDataCell = nanmean(queryRunFeatureData{iQueryOptoParameter}(cellRows));
        
        if ~isnan(queryOptoParameterValues(iQueryOptoParameter))
        
            positionOptoParameterValue = find(queryOptoParameterValuesOverall == queryOptoParameterValues(iQueryOptoParameter));
        
        else
            
            positionOptoParameterValue = find(isnan(queryOptoParameterValuesOverall));
            
        end

        scatter(positionOptoParameterValue-nNaNs, meanQueryRunFeatureDataCell,'k','filled');
        hold on;
        
        meansQueryRunFeatureDataAllCells(find(selectedCellsIDs == cellNames(iCell)),positionOptoParameterValue) = meanQueryRunFeatureDataCell;
        
    end
    
end

meansQueryRunFeatureDataAllCells1 = meansQueryRunFeatureDataAllCells;

% Indicate main directory and add it to path
matlabDirectory = ['/Users/' user '/Desktop/Thesis/Matlab/'];
mainDirectory = ['/Users/' user '/Desktop/Thesis/Matlab/Experiments/']; % including all Matlab toolboxes and necessary paths
addpath(genpath(matlabDirectory));

% Indicate experiment to analyze
experimentName = 'Pitx_GtACR_Ephys_2';

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
runGroupsParametersTable = load(mouseParametersTablePath);
runGroupsParametersTable = runGroupsParametersTable.runGroupsParametersTable;

typeCC = 'withLight';
savePlotPath = [experimentDirectory filesep 'mouseAnalysis'];

queryOptoParameter =  'widthLightPulse';
queryRunFeature = 'nAPduringLight';

queryOptoParameterValuesOverall = [NaN, 0, 500, 1000, 2500, 5000];
nNaNs = numel(find(isnan(queryOptoParameterValuesOverall)));

selectedCellsIDs = [4,5];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

% cell 3: 0.5-0.1 current 200
% cell 2:
% cell 5: 0.5-0.1 current 250

selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.widthCurrentPulse) == 12000 & cell2mat(runGroupsParametersTable.amplitudeLightPulse) == 0.5);
selectedRunGroupsAnalysisTable = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroups),:);

controlRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0') & cell2mat(runGroupsParametersTable.widthCurrentPulse) == 12000); %& cell2mat(runGroupsParametersTable.amplitudeCurrentPulse) == 250
controlRunGroupsAnalysisTable = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,controlRunGroups),:);

queryOptoParameterValues = nan(numel(selectedRunGroups)+numel(controlRunGroups),1);
queryRunFeatureData = {};
cellNameData = {};

groupsCounter = 1;

for iRunGroup = 1:numel(selectedRunGroups)
    
    rowsRunGroup = find(selectedRunGroupsAnalysisTable.runGroup(:) == selectedRunGroups(iRunGroup));
    
    if ~isempty(rowsRunGroup)
    
        queryOptoParameterValues(groupsCounter) = selectedRunGroupsAnalysisTable.optoParameters{rowsRunGroup(1)}.(queryOptoParameter);
        queryRunFeatureData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTable.(queryRunFeature)(rowsRunGroup));
        cellNameData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTable.cellName(rowsRunGroup))

    end
    
    groupsCounter = groupsCounter + 1;
    
end

for iRunGroupControl = 1:numel(controlRunGroups)
    
    rowsControlRunGroup = find(controlRunGroupsAnalysisTable.runGroup(:) == controlRunGroups(iRunGroupControl));
    
    if ~isempty(rowsControlRunGroup)
        
        queryOptoParameterValues(groupsCounter) = controlRunGroupsAnalysisTable.optoParameters{rowsControlRunGroup(1)}.(queryOptoParameter);
        queryRunFeatureData{groupsCounter} = cell2mat(controlRunGroupsAnalysisTable.(queryRunFeature)(rowsControlRunGroup));
        cellNameData{groupsCounter} = cell2mat(controlRunGroupsAnalysisTable.cellName(rowsControlRunGroup));

    end
    
    groupsCounter = groupsCounter + 1;
    
end

allData = vertcat(queryRunFeatureData{:});
cellNameDataArray = vertcat(cellNameData{:});

meansQueryRunFeatureDataAllCells = nan(numel(unique(cellNameDataArray)),numel(queryOptoParameterValuesOverall));


for iQueryOptoParameter = 1:length(queryRunFeatureData)

    cellNames = unique(cellNameData{iQueryOptoParameter})
    
    for iCell = 1:numel(cellNames)
        
        cellRows = find(cellNameData{iQueryOptoParameter} == cellNames(iCell)); 
        meanQueryRunFeatureDataCell = nanmean(queryRunFeatureData{iQueryOptoParameter}(cellRows));
        
        if ~isnan(queryOptoParameterValues(iQueryOptoParameter))
        
            positionOptoParameterValue = find(queryOptoParameterValuesOverall == queryOptoParameterValues(iQueryOptoParameter));
        
        else
            
            positionOptoParameterValue = find(isnan(queryOptoParameterValuesOverall));
            
        end

        scatter(positionOptoParameterValue-nNaNs, meanQueryRunFeatureDataCell,'k','filled');
        hold on;
        
        meansQueryRunFeatureDataAllCells(find(selectedCellsIDs == cellNames(iCell)),positionOptoParameterValue) = meanQueryRunFeatureDataCell;
        
    end
    
end

meansQueryRunFeatureDataAllCells = [meansQueryRunFeatureDataAllCells1; meansQueryRunFeatureDataAllCells]

meanQueryRunFeatureDataAllCells = nanmean(meansQueryRunFeatureDataAllCells,1);
stdQueryRunFeatureDataAllCells = nanstd(meansQueryRunFeatureDataAllCells,1);
xValues = linspace(1,numel(queryOptoParameterValuesOverall)-nNaNs,numel(queryOptoParameterValuesOverall)-nNaNs);
bar(xValues, meanQueryRunFeatureDataAllCells(2:end),'FaceColor',[0 0.4470 0.7410],'BarWidth', 0.5);
errorbar(xValues,meanQueryRunFeatureDataAllCells(2:end),stdQueryRunFeatureDataAllCells(2:end),'.','Color',[0 0 0],'LineWidth', 1.2, 'CapSize', 14)

hScatter = findobj(gca, 'Type', 'scatter');
uistack(hScatter, 'top');

xticks(xValues)
xticklabels(queryOptoParameterValuesOverall(~isnan(queryOptoParameterValuesOverall)))
xlabel(queryOptoParameter);
ylabel(queryRunFeature);
box on;
title([queryOptoParameter,' vs ', queryRunFeature]);
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');

plotName = [experimentName '_' typeCC '_' queryOptoParameter,'_vs_', queryRunFeature, '_2mice'];
plotFullPath = fullfile(savePlotPath, plotName);
saveas(gcf, plotFullPath, 'pdf');
saveas(gcf, plotFullPath, 'fig');



