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
experimentName = 'vGAT_GtACR_Ephys_1_VC';

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

%% pulseWidthBlue vs heightControlPeak

queryOptoParameter =  'pulseWidthBlue';
queryRunFeature = 'heightControlPeak';

queryOptoParameterValuesOverall = [NaN, 250, 500, 2000, 5000, 10000];
nNaNs = numel(find(isnan(queryOptoParameterValuesOverall)));

selectedCellsIDs = [1,2,3];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.holdingVoltage) == 0 & cell2mat(runGroupsParametersTable.nPulsesBlue) == 1 & cell2mat(runGroupsParametersTable.amplitudeBlue) == 5);
selectedRunGroups = selectedRunGroups(selectedRunGroups~=15)
selectedRunGroups = selectedRunGroups(selectedRunGroups~=16)
selectedRunGroupsAnalysisTable = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroups),:);

% controlRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0') & cell2mat(runGroupsParametersTable.widthCurrentPulse) == 12000); %& cell2mat(runGroupsParametersTable.amplitudeCurrentPulse) == 250
% controlRunGroupsAnalysisTable = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,controlRunGroups),:);

queryOptoParameterValues = nan(numel(selectedRunGroups),1);
queryRunFeatureData = {};
cellNameData = {};

groupsCounter = 1;

for iRunGroup = 1:numel(selectedRunGroups)
    
    rowsRunGroup = find(selectedRunGroupsAnalysisTable.runGroup(:) == selectedRunGroups(iRunGroup));
    
    if ~isempty(rowsRunGroup)
    
        queryOptoParameterValues(groupsCounter) = selectedRunGroupsAnalysisTable.optoParameters{rowsRunGroup(1)}.(queryOptoParameter);
        queryRunFeatureData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTable.(queryRunFeature)(rowsRunGroup));
        cellNameData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTable.cellName(rowsRunGroup));

    end
    
    groupsCounter = groupsCounter + 1;
    
end

% for iRunGroupControl = 1:numel(controlRunGroups)
%     
%     rowsControlRunGroup = find(controlRunGroupsAnalysisTable.runGroup(:) == controlRunGroups(iRunGroupControl));
%     
%     if ~isempty(rowsControlRunGroup)
%         
%         queryOptoParameterValues(groupsCounter) = controlRunGroupsAnalysisTable.optoParameters{rowsControlRunGroup(1)}.(queryOptoParameter);
%         queryRunFeatureData{groupsCounter} = cell2mat(controlRunGroupsAnalysisTable.(queryRunFeature)(rowsControlRunGroup));
%         cellNameData{groupsCounter} = cell2mat(controlRunGroupsAnalysisTable.cellName(rowsControlRunGroup));
% 
%     end
%     
%     groupsCounter = groupsCounter + 1;
%     
% end

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

meanQueryRunFeatureDataAllCells = nanmean(meansQueryRunFeatureDataAllCells,1);
stdQueryRunFeatureDataAllCells = nanstd(meansQueryRunFeatureDataAllCells,1);
xValues = linspace(1,numel(queryOptoParameterValuesOverall)-nNaNs,numel(queryOptoParameterValuesOverall)-nNaNs);
bar(xValues, meanQueryRunFeatureDataAllCells(2:end),'FaceColor',[0 0.4470 0.7410],'BarWidth', 0.5);
errorbar(xValues,meanQueryRunFeatureDataAllCells(2:end),stdQueryRunFeatureDataAllCells(2:end),'.','Color',[0.3 0.3 0.3],'LineWidth', 1.2, 'CapSize', 14)

hScatter = findobj(gca, 'Type', 'scatter');
uistack(hScatter, 'top');

xticks(xValues)
xticklabels(queryOptoParameterValuesOverall(~isnan(queryOptoParameterValuesOverall)))
xlabel(queryOptoParameter);
ylabel(queryRunFeature);
box on;
title([queryOptoParameter,' vs ', queryRunFeature]);
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');

plotName = [experimentName '_' typeCC '_' queryOptoParameter,'_vs_', queryRunFeature];
plotFullPath = fullfile(savePlotPath, plotName);
saveas(gcf, plotFullPath, 'pdf');
saveas(gcf, plotFullPath, 'fig');

%% pulseWidthBlue vs heightPulsePeak

queryOptoParameter =  'pulseWidthBlue';
queryRunFeature = 'heightPulsePeak';

queryOptoParameterValuesOverall = [NaN, 250, 500, 2000, 5000, 10000];
nNaNs = numel(find(isnan(queryOptoParameterValuesOverall)));

selectedCellsIDs = [1,2,3];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.holdingVoltage) == 0 & cell2mat(runGroupsParametersTable.nPulsesBlue) == 1 & cell2mat(runGroupsParametersTable.amplitudeBlue) == 5);
selectedRunGroups = selectedRunGroups(selectedRunGroups~=15)
selectedRunGroups = selectedRunGroups(selectedRunGroups~=16)
selectedRunGroupsAnalysisTable = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroups),:);

% controlRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0') & cell2mat(runGroupsParametersTable.widthCurrentPulse) == 12000); %& cell2mat(runGroupsParametersTable.amplitudeCurrentPulse) == 250
% controlRunGroupsAnalysisTable = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,controlRunGroups),:);

queryOptoParameterValues = nan(numel(selectedRunGroups),1);
queryRunFeatureData = {};
cellNameData = {};

groupsCounter = 1;

for iRunGroup = 1:numel(selectedRunGroups)
    
    rowsRunGroup = find(selectedRunGroupsAnalysisTable.runGroup(:) == selectedRunGroups(iRunGroup));
    
    if ~isempty(rowsRunGroup)
    
        queryOptoParameterValues(groupsCounter) = selectedRunGroupsAnalysisTable.optoParameters{rowsRunGroup(1)}.(queryOptoParameter);
        queryRunFeatureData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTable.(queryRunFeature)(rowsRunGroup));
        cellNameData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTable.cellName(rowsRunGroup));

    end
    
    groupsCounter = groupsCounter + 1;
    
end

% for iRunGroupControl = 1:numel(controlRunGroups)
%     
%     rowsControlRunGroup = find(controlRunGroupsAnalysisTable.runGroup(:) == controlRunGroups(iRunGroupControl));
%     
%     if ~isempty(rowsControlRunGroup)
%         
%         queryOptoParameterValues(groupsCounter) = controlRunGroupsAnalysisTable.optoParameters{rowsControlRunGroup(1)}.(queryOptoParameter);
%         queryRunFeatureData{groupsCounter} = cell2mat(controlRunGroupsAnalysisTable.(queryRunFeature)(rowsControlRunGroup));
%         cellNameData{groupsCounter} = cell2mat(controlRunGroupsAnalysisTable.cellName(rowsControlRunGroup));
% 
%     end
%     
%     groupsCounter = groupsCounter + 1;
%     
% end

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

meanQueryRunFeatureDataAllCells = nanmean(meansQueryRunFeatureDataAllCells,1);
stdQueryRunFeatureDataAllCells = nanstd(meansQueryRunFeatureDataAllCells,1);
xValues = linspace(1,numel(queryOptoParameterValuesOverall)-nNaNs,numel(queryOptoParameterValuesOverall)-nNaNs);
bar(xValues, meanQueryRunFeatureDataAllCells(2:end),'FaceColor',[0 0.4470 0.7410],'BarWidth', 0.5);
errorbar(xValues,meanQueryRunFeatureDataAllCells(2:end),stdQueryRunFeatureDataAllCells(2:end),'.','Color',[0.3 0.3 0.3],'LineWidth', 1.2, 'CapSize', 14)

hScatter = findobj(gca, 'Type', 'scatter');
uistack(hScatter, 'top');

xticks(xValues)
xticklabels(queryOptoParameterValuesOverall(~isnan(queryOptoParameterValuesOverall)))
xlabel(queryOptoParameter);
ylabel(queryRunFeature);
box on;
title([queryOptoParameter,' vs ', queryRunFeature]);
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');

plotName = [experimentName '_' typeCC '_' queryOptoParameter,'_vs_', queryRunFeature];
plotFullPath = fullfile(savePlotPath, plotName);
saveas(gcf, plotFullPath, 'pdf');
saveas(gcf, plotFullPath, 'fig');

trace250 = selectedRunGroupsAnalysisTable.trace{44,1};
trace500 = selectedRunGroupsAnalysisTable.trace{50,1};  
trace2000 = selectedRunGroupsAnalysisTable.trace{56,1}; 
trace5000 = selectedRunGroupsAnalysisTable.trace{62,1}; 
trace10000 = selectedRunGroupsAnalysisTable.trace{68,1};  

%% pulseWidthBlue vs areaControl

queryOptoParameter =  'pulseWidthBlue';
queryRunFeature = 'areaControl';

queryOptoParameterValuesOverall = [NaN, 250, 500, 2000, 5000, 10000];
nNaNs = numel(find(isnan(queryOptoParameterValuesOverall)));

selectedCellsIDs = [1, 2, 3];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.holdingVoltage) == 0 & cell2mat(runGroupsParametersTable.nPulsesBlue) == 1 & cell2mat(runGroupsParametersTable.amplitudeBlue) == 5);
selectedRunGroups = selectedRunGroups(selectedRunGroups~=15)
selectedRunGroups = selectedRunGroups(selectedRunGroups~=16)
selectedRunGroupsAnalysisTable = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroups),:);

% controlRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0') & cell2mat(runGroupsParametersTable.widthCurrentPulse) == 12000); %& cell2mat(runGroupsParametersTable.amplitudeCurrentPulse) == 250
% controlRunGroupsAnalysisTable = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,controlRunGroups),:);

queryOptoParameterValues = nan(numel(selectedRunGroups),1);
queryRunFeatureData = {};
cellNameData = {};

groupsCounter = 1;

for iRunGroup = 1:numel(selectedRunGroups)
    
    rowsRunGroup = find(selectedRunGroupsAnalysisTable.runGroup(:) == selectedRunGroups(iRunGroup));
    
    if ~isempty(rowsRunGroup)
    
        queryOptoParameterValues(groupsCounter) = selectedRunGroupsAnalysisTable.optoParameters{rowsRunGroup(1)}.(queryOptoParameter);
        queryRunFeatureData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTable.(queryRunFeature)(rowsRunGroup));
        cellNameData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTable.cellName(rowsRunGroup));

    end
    
    groupsCounter = groupsCounter + 1;
    
end

% for iRunGroupControl = 1:numel(controlRunGroups)
%     
%     rowsControlRunGroup = find(controlRunGroupsAnalysisTable.runGroup(:) == controlRunGroups(iRunGroupControl));
%     
%     if ~isempty(rowsControlRunGroup)
%         
%         queryOptoParameterValues(groupsCounter) = controlRunGroupsAnalysisTable.optoParameters{rowsControlRunGroup(1)}.(queryOptoParameter);
%         queryRunFeatureData{groupsCounter} = cell2mat(controlRunGroupsAnalysisTable.(queryRunFeature)(rowsControlRunGroup));
%         cellNameData{groupsCounter} = cell2mat(controlRunGroupsAnalysisTable.cellName(rowsControlRunGroup));
% 
%     end
%     
%     groupsCounter = groupsCounter + 1;
%     
% end

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

meanQueryRunFeatureDataAllCells = nanmean(meansQueryRunFeatureDataAllCells,1);
stdQueryRunFeatureDataAllCells = nanstd(meansQueryRunFeatureDataAllCells,1);
xValues = linspace(1,numel(queryOptoParameterValuesOverall)-nNaNs,numel(queryOptoParameterValuesOverall)-nNaNs);
bar(xValues, meanQueryRunFeatureDataAllCells(2:end),'FaceColor',[0 0.4470 0.7410],'BarWidth', 0.5);
errorbar(xValues,meanQueryRunFeatureDataAllCells(2:end),stdQueryRunFeatureDataAllCells(2:end),'.','Color',[0.3 0.3 0.3],'LineWidth', 1.2, 'CapSize', 14)

hScatter = findobj(gca, 'Type', 'scatter');
uistack(hScatter, 'top');

xticks(xValues)
xticklabels(queryOptoParameterValuesOverall(~isnan(queryOptoParameterValuesOverall)))
xlabel(queryOptoParameter);
ylabel(queryRunFeature);
box on;
title([queryOptoParameter,' vs ', queryRunFeature]);
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');

plotName = [experimentName '_' typeCC '_' queryOptoParameter,'_vs_', queryRunFeature];
plotFullPath = fullfile(savePlotPath, plotName);
saveas(gcf, plotFullPath, 'pdf');
saveas(gcf, plotFullPath, 'fig');

%% pulseWidthBlue vs areaPulse

queryOptoParameter =  'pulseWidthBlue';
queryRunFeature = 'areaPulse';

queryOptoParameterValuesOverall = [NaN, 250, 500, 2000, 5000, 10000];
nNaNs = numel(find(isnan(queryOptoParameterValuesOverall)));

selectedCellsIDs = [1,2,3];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.holdingVoltage) == 0 & cell2mat(runGroupsParametersTable.nPulsesBlue) == 1 & cell2mat(runGroupsParametersTable.amplitudeBlue) == 5);
selectedRunGroups = selectedRunGroups(selectedRunGroups~=15)
selectedRunGroups = selectedRunGroups(selectedRunGroups~=16)
selectedRunGroupsAnalysisTable = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroups),:);

% controlRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0') & cell2mat(runGroupsParametersTable.widthCurrentPulse) == 12000); %& cell2mat(runGroupsParametersTable.amplitudeCurrentPulse) == 250
% controlRunGroupsAnalysisTable = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,controlRunGroups),:);

queryOptoParameterValues = nan(numel(selectedRunGroups),1);
queryRunFeatureData = {};
cellNameData = {};

groupsCounter = 1;

for iRunGroup = 1:numel(selectedRunGroups)
    
    rowsRunGroup = find(selectedRunGroupsAnalysisTable.runGroup(:) == selectedRunGroups(iRunGroup));
    
    if ~isempty(rowsRunGroup)
    
        queryOptoParameterValues(groupsCounter) = selectedRunGroupsAnalysisTable.optoParameters{rowsRunGroup(1)}.(queryOptoParameter);
        queryRunFeatureData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTable.(queryRunFeature)(rowsRunGroup));
        cellNameData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTable.cellName(rowsRunGroup));

    end
    
    groupsCounter = groupsCounter + 1;
    
end

% for iRunGroupControl = 1:numel(controlRunGroups)
%     
%     rowsControlRunGroup = find(controlRunGroupsAnalysisTable.runGroup(:) == controlRunGroups(iRunGroupControl));
%     
%     if ~isempty(rowsControlRunGroup)
%         
%         queryOptoParameterValues(groupsCounter) = controlRunGroupsAnalysisTable.optoParameters{rowsControlRunGroup(1)}.(queryOptoParameter);
%         queryRunFeatureData{groupsCounter} = cell2mat(controlRunGroupsAnalysisTable.(queryRunFeature)(rowsControlRunGroup));
%         cellNameData{groupsCounter} = cell2mat(controlRunGroupsAnalysisTable.cellName(rowsControlRunGroup));
% 
%     end
%     
%     groupsCounter = groupsCounter + 1;
%     
% end

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

meanQueryRunFeatureDataAllCells = nanmean(meansQueryRunFeatureDataAllCells,1);
stdQueryRunFeatureDataAllCells = nanstd(meansQueryRunFeatureDataAllCells,1);
xValues = linspace(1,numel(queryOptoParameterValuesOverall)-nNaNs,numel(queryOptoParameterValuesOverall)-nNaNs);
bar(xValues, meanQueryRunFeatureDataAllCells(2:end),'FaceColor',[0 0.4470 0.7410],'BarWidth', 0.5);
errorbar(xValues,meanQueryRunFeatureDataAllCells(2:end),stdQueryRunFeatureDataAllCells(2:end),'.','Color',[0.3 0.3 0.3],'LineWidth', 1.2, 'CapSize', 14)

hScatter = findobj(gca, 'Type', 'scatter');
uistack(hScatter, 'top');

xticks(xValues)
xticklabels(queryOptoParameterValuesOverall(~isnan(queryOptoParameterValuesOverall)))
xlabel(queryOptoParameter);
ylabel(queryRunFeature);
box on;
title([queryOptoParameter,' vs ', queryRunFeature]);
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');

plotName = [experimentName '_' typeCC '_' queryOptoParameter,'_vs_', queryRunFeature];
plotFullPath = fullfile(savePlotPath, plotName);
saveas(gcf, plotFullPath, 'pdf');
saveas(gcf, plotFullPath, 'fig');

%% amplitudeBlue vs heightControlPeak

queryOptoParameter =  'amplitudeBlue';
queryRunFeature = 'heightControlPeak';

queryOptoParameterValuesOverall = [NaN, 1, 3, 5];
nNaNs = numel(find(isnan(queryOptoParameterValuesOverall)));

selectedCellsIDs = [1,2,3];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.holdingVoltage) == 0 & cell2mat(runGroupsParametersTable.nPulsesBlue) == 1 & cell2mat(runGroupsParametersTable.pulseWidthBlue) == 5000);
selectedRunGroups = selectedRunGroups(selectedRunGroups~=15)
selectedRunGroups = selectedRunGroups(selectedRunGroups~=16)
selectedRunGroupsAnalysisTable = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroups),:);

% controlRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0') & cell2mat(runGroupsParametersTable.widthCurrentPulse) == 12000); %& cell2mat(runGroupsParametersTable.amplitudeCurrentPulse) == 250
% controlRunGroupsAnalysisTable = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,controlRunGroups),:);

queryOptoParameterValues = nan(numel(selectedRunGroups),1);
queryRunFeatureData = {};
cellNameData = {};

groupsCounter = 1;

for iRunGroup = 1:numel(selectedRunGroups)
    
    rowsRunGroup = find(selectedRunGroupsAnalysisTable.runGroup(:) == selectedRunGroups(iRunGroup));
    
    if ~isempty(rowsRunGroup)
    
        queryOptoParameterValues(groupsCounter) = selectedRunGroupsAnalysisTable.optoParameters{rowsRunGroup(1)}.(queryOptoParameter);
        queryRunFeatureData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTable.(queryRunFeature)(rowsRunGroup));
        cellNameData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTable.cellName(rowsRunGroup));

    end
    
    groupsCounter = groupsCounter + 1;
    
end

% for iRunGroupControl = 1:numel(controlRunGroups)
%     
%     rowsControlRunGroup = find(controlRunGroupsAnalysisTable.runGroup(:) == controlRunGroups(iRunGroupControl));
%     
%     if ~isempty(rowsControlRunGroup)
%         
%         queryOptoParameterValues(groupsCounter) = controlRunGroupsAnalysisTable.optoParameters{rowsControlRunGroup(1)}.(queryOptoParameter);
%         queryRunFeatureData{groupsCounter} = cell2mat(controlRunGroupsAnalysisTable.(queryRunFeature)(rowsControlRunGroup));
%         cellNameData{groupsCounter} = cell2mat(controlRunGroupsAnalysisTable.cellName(rowsControlRunGroup));
% 
%     end
%     
%     groupsCounter = groupsCounter + 1;
%     
% end

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

meanQueryRunFeatureDataAllCells = nanmean(meansQueryRunFeatureDataAllCells,1);
stdQueryRunFeatureDataAllCells = nanstd(meansQueryRunFeatureDataAllCells,1);
xValues = linspace(1,numel(queryOptoParameterValuesOverall)-nNaNs,numel(queryOptoParameterValuesOverall)-nNaNs);
bar(xValues, meanQueryRunFeatureDataAllCells(2:end),'FaceColor',[0 0.4470 0.7410],'BarWidth', 0.5);
errorbar(xValues,meanQueryRunFeatureDataAllCells(2:end),stdQueryRunFeatureDataAllCells(2:end),'.','Color',[0.3 0.3 0.3],'LineWidth', 1.2, 'CapSize', 14)

hScatter = findobj(gca, 'Type', 'scatter');
uistack(hScatter, 'top');

xticks(xValues)
xticklabels(queryOptoParameterValuesOverall(~isnan(queryOptoParameterValuesOverall)))
xlabel(queryOptoParameter);
ylabel(queryRunFeature);
box on;
title([queryOptoParameter,' vs ', queryRunFeature]);
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');

plotName = [experimentName '_' typeCC '_' queryOptoParameter,'_vs_', queryRunFeature];
plotFullPath = fullfile(savePlotPath, plotName);
saveas(gcf, plotFullPath, 'pdf');
saveas(gcf, plotFullPath, 'fig');

%% amplitudeBlue vs heightPulsePeak

queryOptoParameter =  'amplitudeBlue';
queryRunFeature = 'heightPulsePeak';

queryOptoParameterValuesOverall = [NaN, 1, 3, 5];
nNaNs = numel(find(isnan(queryOptoParameterValuesOverall)));

selectedCellsIDs = [1,2,3];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.holdingVoltage) == 0 & cell2mat(runGroupsParametersTable.nPulsesBlue) == 1 & cell2mat(runGroupsParametersTable.pulseWidthBlue) == 5000);
selectedRunGroups = selectedRunGroups(selectedRunGroups~=15)
selectedRunGroups = selectedRunGroups(selectedRunGroups~=16)
selectedRunGroupsAnalysisTable = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroups),:);

% controlRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0') & cell2mat(runGroupsParametersTable.widthCurrentPulse) == 12000); %& cell2mat(runGroupsParametersTable.amplitudeCurrentPulse) == 250
% controlRunGroupsAnalysisTable = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,controlRunGroups),:);

queryOptoParameterValues = nan(numel(selectedRunGroups),1);
queryRunFeatureData = {};
cellNameData = {};

groupsCounter = 1;

for iRunGroup = 1:numel(selectedRunGroups)
    
    rowsRunGroup = find(selectedRunGroupsAnalysisTable.runGroup(:) == selectedRunGroups(iRunGroup));
    
    if ~isempty(rowsRunGroup)
    
        queryOptoParameterValues(groupsCounter) = selectedRunGroupsAnalysisTable.optoParameters{rowsRunGroup(1)}.(queryOptoParameter);
        queryRunFeatureData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTable.(queryRunFeature)(rowsRunGroup));
        cellNameData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTable.cellName(rowsRunGroup));

    end
    
    groupsCounter = groupsCounter + 1;
    
end

% for iRunGroupControl = 1:numel(controlRunGroups)
%     
%     rowsControlRunGroup = find(controlRunGroupsAnalysisTable.runGroup(:) == controlRunGroups(iRunGroupControl));
%     
%     if ~isempty(rowsControlRunGroup)
%         
%         queryOptoParameterValues(groupsCounter) = controlRunGroupsAnalysisTable.optoParameters{rowsControlRunGroup(1)}.(queryOptoParameter);
%         queryRunFeatureData{groupsCounter} = cell2mat(controlRunGroupsAnalysisTable.(queryRunFeature)(rowsControlRunGroup));
%         cellNameData{groupsCounter} = cell2mat(controlRunGroupsAnalysisTable.cellName(rowsControlRunGroup));
% 
%     end
%     
%     groupsCounter = groupsCounter + 1;
%     
% end

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

meanQueryRunFeatureDataAllCells = nanmean(meansQueryRunFeatureDataAllCells,1);
stdQueryRunFeatureDataAllCells = nanstd(meansQueryRunFeatureDataAllCells,1);
xValues = linspace(1,numel(queryOptoParameterValuesOverall)-nNaNs,numel(queryOptoParameterValuesOverall)-nNaNs);
bar(xValues, meanQueryRunFeatureDataAllCells(2:end),'FaceColor',[0 0.4470 0.7410],'BarWidth', 0.5);
errorbar(xValues,meanQueryRunFeatureDataAllCells(2:end),stdQueryRunFeatureDataAllCells(2:end),'.','Color',[0.3 0.3 0.3],'LineWidth', 1.2, 'CapSize', 14)

hScatter = findobj(gca, 'Type', 'scatter');
uistack(hScatter, 'top');

xticks(xValues)
xticklabels(queryOptoParameterValuesOverall(~isnan(queryOptoParameterValuesOverall)))
xlabel(queryOptoParameter);
ylabel(queryRunFeature);
box on;
title([queryOptoParameter,' vs ', queryRunFeature]);
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');

plotName = [experimentName '_' typeCC '_' queryOptoParameter,'_vs_', queryRunFeature];
plotFullPath = fullfile(savePlotPath, plotName);
saveas(gcf, plotFullPath, 'pdf');
saveas(gcf, plotFullPath, 'fig');

trace5 = selectedRunGroupsAnalysisTable.trace{4,1};
trace3 = selectedRunGroupsAnalysisTable.trace{9,1};  
trace1 = selectedRunGroupsAnalysisTable.trace{14,1}; 

%% amplitudeBlue vs areaControl

queryOptoParameter =  'amplitudeBlue';
queryRunFeature = 'areaControl';

queryOptoParameterValuesOverall = [NaN, 1, 3, 5];
nNaNs = numel(find(isnan(queryOptoParameterValuesOverall)));

selectedCellsIDs = [1,2,3];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.holdingVoltage) == 0 & cell2mat(runGroupsParametersTable.nPulsesBlue) == 1 & cell2mat(runGroupsParametersTable.pulseWidthBlue) == 5000);
selectedRunGroups = selectedRunGroups(selectedRunGroups~=15)
selectedRunGroups = selectedRunGroups(selectedRunGroups~=16)
selectedRunGroupsAnalysisTable = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroups),:);

% controlRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0') & cell2mat(runGroupsParametersTable.widthCurrentPulse) == 12000); %& cell2mat(runGroupsParametersTable.amplitudeCurrentPulse) == 250
% controlRunGroupsAnalysisTable = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,controlRunGroups),:);

queryOptoParameterValues = nan(numel(selectedRunGroups),1);
queryRunFeatureData = {};
cellNameData = {};

groupsCounter = 1;

for iRunGroup = 1:numel(selectedRunGroups)
    
    rowsRunGroup = find(selectedRunGroupsAnalysisTable.runGroup(:) == selectedRunGroups(iRunGroup));
    
    if ~isempty(rowsRunGroup)
    
        queryOptoParameterValues(groupsCounter) = selectedRunGroupsAnalysisTable.optoParameters{rowsRunGroup(1)}.(queryOptoParameter);
        queryRunFeatureData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTable.(queryRunFeature)(rowsRunGroup));
        cellNameData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTable.cellName(rowsRunGroup));

    end
    
    groupsCounter = groupsCounter + 1;
    
end

% for iRunGroupControl = 1:numel(controlRunGroups)
%     
%     rowsControlRunGroup = find(controlRunGroupsAnalysisTable.runGroup(:) == controlRunGroups(iRunGroupControl));
%     
%     if ~isempty(rowsControlRunGroup)
%         
%         queryOptoParameterValues(groupsCounter) = controlRunGroupsAnalysisTable.optoParameters{rowsControlRunGroup(1)}.(queryOptoParameter);
%         queryRunFeatureData{groupsCounter} = cell2mat(controlRunGroupsAnalysisTable.(queryRunFeature)(rowsControlRunGroup));
%         cellNameData{groupsCounter} = cell2mat(controlRunGroupsAnalysisTable.cellName(rowsControlRunGroup));
% 
%     end
%     
%     groupsCounter = groupsCounter + 1;
%     
% end

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

meanQueryRunFeatureDataAllCells = nanmean(meansQueryRunFeatureDataAllCells,1);
stdQueryRunFeatureDataAllCells = nanstd(meansQueryRunFeatureDataAllCells,1);
xValues = linspace(1,numel(queryOptoParameterValuesOverall)-nNaNs,numel(queryOptoParameterValuesOverall)-nNaNs);
bar(xValues, meanQueryRunFeatureDataAllCells(2:end),'FaceColor',[0 0.4470 0.7410],'BarWidth', 0.5);
errorbar(xValues,meanQueryRunFeatureDataAllCells(2:end),stdQueryRunFeatureDataAllCells(2:end),'.','Color',[0.3 0.3 0.3],'LineWidth', 1.2, 'CapSize', 14)

hScatter = findobj(gca, 'Type', 'scatter');
uistack(hScatter, 'top');

xticks(xValues)
xticklabels(queryOptoParameterValuesOverall(~isnan(queryOptoParameterValuesOverall)))
xlabel(queryOptoParameter);
ylabel(queryRunFeature);
box on;
title([queryOptoParameter,' vs ', queryRunFeature]);
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');

plotName = [experimentName '_' typeCC '_' queryOptoParameter,'_vs_', queryRunFeature];
plotFullPath = fullfile(savePlotPath, plotName);
saveas(gcf, plotFullPath, 'pdf');
saveas(gcf, plotFullPath, 'fig');

%% amplitudeBlue vs areaPulse

queryOptoParameter =  'amplitudeBlue';
queryRunFeature = 'areaPulse';

queryOptoParameterValuesOverall = [NaN, 1, 3, 5];
nNaNs = numel(find(isnan(queryOptoParameterValuesOverall)));

selectedCellsIDs = [1,2,3];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.holdingVoltage) == 0 & cell2mat(runGroupsParametersTable.nPulsesBlue) == 1 & cell2mat(runGroupsParametersTable.pulseWidthBlue) == 5000);
selectedRunGroups = selectedRunGroups(selectedRunGroups~=15)
selectedRunGroups = selectedRunGroups(selectedRunGroups~=16)
selectedRunGroupsAnalysisTable = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroups),:);

% controlRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0') & cell2mat(runGroupsParametersTable.widthCurrentPulse) == 12000); %& cell2mat(runGroupsParametersTable.amplitudeCurrentPulse) == 250
% controlRunGroupsAnalysisTable = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,controlRunGroups),:);

queryOptoParameterValues = nan(numel(selectedRunGroups),1);
queryRunFeatureData = {};
cellNameData = {};

groupsCounter = 1;

for iRunGroup = 1:numel(selectedRunGroups)
    
    rowsRunGroup = find(selectedRunGroupsAnalysisTable.runGroup(:) == selectedRunGroups(iRunGroup));
    
    if ~isempty(rowsRunGroup)
    
        queryOptoParameterValues(groupsCounter) = selectedRunGroupsAnalysisTable.optoParameters{rowsRunGroup(1)}.(queryOptoParameter);
        queryRunFeatureData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTable.(queryRunFeature)(rowsRunGroup));
        cellNameData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTable.cellName(rowsRunGroup));

    end
    
    groupsCounter = groupsCounter + 1;
    
end

% for iRunGroupControl = 1:numel(controlRunGroups)
%     
%     rowsControlRunGroup = find(controlRunGroupsAnalysisTable.runGroup(:) == controlRunGroups(iRunGroupControl));
%     
%     if ~isempty(rowsControlRunGroup)
%         
%         queryOptoParameterValues(groupsCounter) = controlRunGroupsAnalysisTable.optoParameters{rowsControlRunGroup(1)}.(queryOptoParameter);
%         queryRunFeatureData{groupsCounter} = cell2mat(controlRunGroupsAnalysisTable.(queryRunFeature)(rowsControlRunGroup));
%         cellNameData{groupsCounter} = cell2mat(controlRunGroupsAnalysisTable.cellName(rowsControlRunGroup));
% 
%     end
%     
%     groupsCounter = groupsCounter + 1;
%     
% end

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

meanQueryRunFeatureDataAllCells = nanmean(meansQueryRunFeatureDataAllCells,1);
stdQueryRunFeatureDataAllCells = nanstd(meansQueryRunFeatureDataAllCells,1);
xValues = linspace(1,numel(queryOptoParameterValuesOverall)-nNaNs,numel(queryOptoParameterValuesOverall)-nNaNs);
bar(xValues, meanQueryRunFeatureDataAllCells(2:end),'FaceColor',[0 0.4470 0.7410],'BarWidth', 0.5);
errorbar(xValues,meanQueryRunFeatureDataAllCells(2:end),stdQueryRunFeatureDataAllCells(2:end),'.','Color',[0.3 0.3 0.3],'LineWidth', 1.2, 'CapSize', 14)

hScatter = findobj(gca, 'Type', 'scatter');
uistack(hScatter, 'top');

xticks(xValues)
xticklabels(queryOptoParameterValuesOverall(~isnan(queryOptoParameterValuesOverall)))
xlabel(queryOptoParameter);
ylabel(queryRunFeature);
box on;
title([queryOptoParameter,' vs ', queryRunFeature]);
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');

plotName = [experimentName '_' typeCC '_' queryOptoParameter,'_vs_', queryRunFeature];
plotFullPath = fullfile(savePlotPath, plotName);
saveas(gcf, plotFullPath, 'pdf');
saveas(gcf, plotFullPath, 'fig');

%% comparison traces widthLightPulse

typeCC = 'comparision_widthLightPulse_traces';

figure;
plot(trace250,'LineWidth',2)
hold on;
plot(trace500,'LineWidth',2)
hold on;
plot(trace2000,'LineWidth',2)
hold on;
plot(trace5000,'LineWidth',2)
hold on;
plot(trace10000,'LineWidth',2)

xLimMin = 0;
xLimMax = 2.5;
xticks = linspace(xLimMin,xLimMax, 6);
xlim([xLimMin*10000, xLimMax*10000]);
set(gca, 'XTickLabel', num2str((xticks'-xLimMin)))
xlabel('Time [s]');
ylabel('Current [pA]');
box on;
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');
legend('25 ms','50 ms','200 ms','500 ms','1000 ms','Location', 'northeast', 'FontSize',12)

plotName = [experimentName '_' typeCC];
plotFullPath = fullfile(savePlotPath, plotName);
saveas(gcf, plotFullPath, 'pdf');
saveas(gcf, plotFullPath, 'fig');

%% comparison traces amplitudeLightPulse

typeCC = 'comparision_amplitudeLightPulse_traces';

figure;
plot(trace1,'LineWidth',2)
hold on;
plot(trace3,'LineWidth',2)
hold on;
plot(trace5,'LineWidth',2)

xLimMin = 0;
xLimMax = 2;
xticks = linspace(xLimMin,xLimMax, 11);
xlim([xLimMin*10000, xLimMax*10000]);
set(gca, 'XTickLabel', num2str((xticks'-xLimMin)))
xlabel('Time [s]');
ylabel('Current [pA]');
box on;
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');
legend('1 mV','3 mV','5 mV','Location', 'northeast', 'FontSize',12)

plotName = [experimentName '_' typeCC];
plotFullPath = fullfile(savePlotPath, plotName);
saveas(gcf, plotFullPath, 'pdf');
saveas(gcf, plotFullPath, 'fig');