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
experimentName = 'vGAT_ZI_Ephys_11';

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

savePlotPath = [experimentDirectory filesep 'mouseAnalysis'];
savePlots = 0;

%% amplitudeRed vs heightPulsePeak

queryOptoParameter =  'amplitudeRed';
queryRunFeature = 'heightPulsePeak';

queryOptoParameterValuesOverall = [5];
nNaNs = numel(find(isnan(queryOptoParameterValuesOverall)));

selectedCellsIDs = [1,2,7,8];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao2') & cell2mat(runGroupsParametersTable.pulseWidthRed) == 100 & strcmp(runGroupsParametersTable.functionNameRed,'dmdOutputRect(1, 1, 1, 1, 2, 1)')& cell2mat(runGroupsParametersTable.nPulsesRed) == 1);
selectedRunGroupsAnalysisTable = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroups),:);

% desiredEpochsDrugsRows = find(cell2mat(selectedRunGroupsAnalysisTable.epoch) >= 36);
% selectedRunGroupsAnalysisTableDrugs = selectedRunGroupsAnalysisTable(desiredEpochsDrugsRows,:)

% desiredEpochsRows = find(cell2mat(selectedRunGroupsAnalysisTable.epoch) > 3 & cell2mat(selectedRunGroupsAnalysisTable.epoch) <= 13);
% selectedRunGroupsAnalysisTable = selectedRunGroupsAnalysisTable(desiredEpochsRows,:)

queryOptoParameterValues = nan(numel(selectedRunGroups),1);
queryRunFeatureData = {};
cellNameData = {};
allTraces = [];

groupsCounter = 1;

for iRunGroup = 1:numel(selectedRunGroups)
    
    rowsRunGroup = find(selectedRunGroupsAnalysisTable.runGroup(:) == selectedRunGroups(iRunGroup));
    
    if ~isempty(rowsRunGroup)
    
        queryOptoParameterValues(groupsCounter) = selectedRunGroupsAnalysisTable.optoParameters{rowsRunGroup(1)}.(queryOptoParameter);
        queryRunFeatureData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTable.(queryRunFeature)(rowsRunGroup));
        cellNameData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTable.cellName(rowsRunGroup));

        tracesGroup = selectedRunGroupsAnalysisTable.trace(rowsRunGroup)
        allTracesGroup = [];
        
        for iTrace = 1:numel(tracesGroup)
            
            trace = tracesGroup{iTrace};
            trace = preprocessSignalVCv2(trace);
            
            if numel(trace) < 10000
                
                trace = [zeros(1,5000) trace zeros(1,5000)];
                
            end
            
            allTracesGroup = [allTracesGroup; trace];
            
        end
        
        allTraces = [allTraces; allTracesGroup];
    
    end
    
    groupsCounter = groupsCounter + 1;
    
end

allData = vertcat(queryRunFeatureData{:});
meanAllData = nanmean(allData,1);
stdAllData = nanstd(allData,1);
cellNameDataArray = vertcat(cellNameData{:});

meansQueryRunFeatureDataAllCells = nan(numel(unique(cellNameDataArray)),numel(queryOptoParameterValuesOverall));

for iQueryOptoParameter = 1:length(queryRunFeatureData)

    cellNames = unique(cellNameData{iQueryOptoParameter})
    
    for iCell = 1:numel(cellNames)
        
        cellRows = find(cellNameData{iQueryOptoParameter} == cellNames(iCell)); 
        meanQueryRunFeatureDataCell = nanmean(queryRunFeatureData{iQueryOptoParameter}(cellRows))
        
        if ~isnan(queryOptoParameterValues(iQueryOptoParameter))
        
            positionOptoParameterValue = find(queryOptoParameterValuesOverall == queryOptoParameterValues(iQueryOptoParameter));
        
        else
            
            positionOptoParameterValue = find(isnan(queryOptoParameterValuesOverall));
            
        end
        
        meansQueryRunFeatureDataAllCells(find(selectedCellsIDs == cellNames(iCell)),positionOptoParameterValue) = meanQueryRunFeatureDataCell;
        
    end
    
end

meanQueryRunFeatureDataAllCells = nanmean(meansQueryRunFeatureDataAllCells,1);
stdQueryRunFeatureDataAllCells = nanstd(meansQueryRunFeatureDataAllCells,1);

selectedCellsIDsRest = [3,4,5,6];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDsRest));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

selectedRunGroupsRest = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao2') & cell2mat(runGroupsParametersTable.holdingVoltage) == 0 & cell2mat(runGroupsParametersTable.pulseWidthRed) == 100 & strcmp(runGroupsParametersTable.functionNameRed,'dmdOutputRect(1, 1, 1, 1, 2, 1)')& cell2mat(runGroupsParametersTable.nPulsesRed) == 1);
selectedRunGroupsAnalysisTableRest = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroupsRest),:);

queryRunFeatureData = {};
cellNameData = {};
allTracesRest = [];

groupsCounter = 1;

for iRunGroup = 1:numel(selectedRunGroupsRest)
    
    rowsRunGroup = find(selectedRunGroupsAnalysisTableRest.runGroup(:) == selectedRunGroupsRest(iRunGroup));
    
    if ~isempty(rowsRunGroup)
    
        queryOptoParameterValues(groupsCounter) = selectedRunGroupsAnalysisTableRest.optoParameters{rowsRunGroup(1)}.(queryOptoParameter);
        queryRunFeatureData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTableRest.(queryRunFeature)(rowsRunGroup));
        cellNameData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTableRest.cellName(rowsRunGroup));
        tracesGroup = selectedRunGroupsAnalysisTableRest.trace(rowsRunGroup)
        allTracesGroup = [];
        
        for iTrace = 1:numel(tracesGroup)
            
            trace = tracesGroup{iTrace};
            
            if numel(trace) < 10000
                
                trace = [zeros(1,5000) trace zeros(1,5000)];
                
            end
            
            allTracesGroup = [allTracesGroup; trace];
            
        end
        
        allTracesRest = [allTracesRest; allTracesGroup];
    
    end
    
    groupsCounter = groupsCounter + 1;
    
end

allDataRest = vertcat(queryRunFeatureData{:});
meanAllDataRest = nanmean(allDataRest,1);
stdAllDataRest = nanstd(allDataRest,1);
cellNameDataArray = vertcat(cellNameData{:});

meansQueryRunFeatureDataRestAllCells = nan(numel(unique(cellNameDataArray)),numel(queryOptoParameterValuesOverall));

for iQueryOptoParameter = 1:length(queryRunFeatureData)

    cellNames = unique(cellNameData{iQueryOptoParameter})
    
    for iCell = 1:numel(cellNames)
        
        cellRows = find(cellNameData{iQueryOptoParameter} == cellNames(iCell)); 
        meanQueryRunFeatureDataCell = nanmean(queryRunFeatureData{iQueryOptoParameter}(cellRows))
        
        if ~isnan(queryOptoParameterValues(iQueryOptoParameter))
        
            positionOptoParameterValue = find(queryOptoParameterValuesOverall == queryOptoParameterValues(iQueryOptoParameter));
        
        else
            
            positionOptoParameterValue = find(isnan(queryOptoParameterValuesOverall));
            
        end
        
        meansQueryRunFeatureDataRestAllCells(find(selectedCellsIDsRest == cellNames(iCell)),positionOptoParameterValue) = meanQueryRunFeatureDataCell;
        
    end
    
end

meanQueryRunFeatureDataRestAllCells = nanmean(meansQueryRunFeatureDataRestAllCells,1);
stdQueryRunFeatureDataRestAllCells = nanstd(meansQueryRunFeatureDataRestAllCells,1);

figure
scatter([0], meansQueryRunFeatureDataRestAllCells,'k','filled');
hold on
scatter([1], meansQueryRunFeatureDataAllCells,'k','filled');
hold on
bar([0,1], [meanQueryRunFeatureDataRestAllCells, meanQueryRunFeatureDataAllCells] ,'FaceColor',[0 0.4470 0.7410],'BarWidth', 0.5);
errorbar([0,1],[meanQueryRunFeatureDataRestAllCells, meanQueryRunFeatureDataAllCells],[stdQueryRunFeatureDataRestAllCells, stdQueryRunFeatureDataAllCells],'.','Color',[0.3 0.3 0.3],'LineWidth', 1.2, 'CapSize', 14)

hScatter = findobj(gca, 'Type', 'scatter');
uistack(hScatter, 'top');

xticks([0,1])
xticklabels({'Non-connected cells','Cells 1-2-7-8'})
xlabel(['Cells']);
ylabel(queryRunFeature);
ylim([-10,70])
box on;
title([queryOptoParameter,' vs ', queryRunFeature]);
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');

plotName = [experimentName '_', queryOptoParameter , '_vs_', queryRunFeature, '_', 'cellsComparison'];
plotFullPath = fullfile(savePlotPath, plotName);

if savePlots == 1
    saveas(gcf, plotFullPath, 'pdf')
    saveas(gcf, plotFullPath, 'fig')
end

figure;
plot(nanmean(allTraces,1),'Color',[.6, .1, .1],'LineWidth',2)
hold on

errorUp = nanmean(allTraces,1) + nanstd(allTraces,1);
errorDown = nanmean(allTraces,1) - nanstd(allTraces,1);
xLinspace = linspace(1,size(allTraces,2), size(allTraces,2));
plotpatch = patch([xLinspace fliplr(xLinspace)], [errorDown fliplr(errorUp)], [.6, .2, .2], 'FaceAlpha',0.5, 'EdgeColor','none');
plotpatch.Annotation.LegendInformation.IconDisplayStyle = 'off';

plot(nanmean(allTracesRest,1),'Color',[0.4 0.4 0.4],'LineWidth',2)

errorUp = nanmean(allTracesRest,1) + nanstd(allTracesRest,1);
errorDown = nanmean(allTracesRest,1) - nanstd(allTracesRest,1);
xLinspace = linspace(1,size(allTracesRest,2), size(allTracesRest,2));
plotpatch = patch([xLinspace fliplr(xLinspace)], [errorDown fliplr(errorUp)], [.8, .8, .8], 'FaceAlpha',0.5, 'EdgeColor','none');
plotpatch.Annotation.LegendInformation.IconDisplayStyle = 'off';

%xlim([4500,6500])
ylim([-10,40])
%xticks([4500,5500,6500])
xticklabels({'0','0.1','0.2'})
xlabel('Time [s]')
ylabel('Current [pA]')
box on;
title([queryOptoParameter,' vs ', queryRunFeature]);
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');
legend('Cells 1-2-7-8','Non-connected cells','Location', 'northeast', 'FontSize',12)

plotName = [experimentName '_', queryOptoParameter , '_vs_', queryRunFeature, '_', 'cellsComparison_', 'traces'];
plotFullPath = fullfile(savePlotPath, plotName);

if savePlots == 1
    saveas(gcf, plotFullPath, 'pdf')
    saveas(gcf, plotFullPath, 'fig')
end

%% amplitudeRed vs areaPulse

queryOptoParameter =  'amplitudeRed';
queryRunFeature = 'areaPulse';

queryOptoParameterValuesOverall = [5];
nNaNs = numel(find(isnan(queryOptoParameterValuesOverall)));

selectedCellsIDs = [1,2,7,8];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao2') & cell2mat(runGroupsParametersTable.holdingVoltage) == 0 & cell2mat(runGroupsParametersTable.pulseWidthRed) == 100 & strcmp(runGroupsParametersTable.functionNameRed,'dmdOutputRect(1, 1, 1, 1, 2, 1)')& cell2mat(runGroupsParametersTable.nPulsesRed) == 1);
selectedRunGroupsAnalysisTable = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroups),:);

% desiredEpochsDrugsRows = find(cell2mat(selectedRunGroupsAnalysisTable.epoch) >= 36);
% selectedRunGroupsAnalysisTableDrugs = selectedRunGroupsAnalysisTable(desiredEpochsDrugsRows,:)

% desiredEpochsRows = find(cell2mat(selectedRunGroupsAnalysisTable.epoch) > 3 & cell2mat(selectedRunGroupsAnalysisTable.epoch) <= 13);
% selectedRunGroupsAnalysisTable = selectedRunGroupsAnalysisTable(desiredEpochsRows,:)

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

allData = vertcat(queryRunFeatureData{:});
meanAllData = nanmean(allData,1);
stdAllData = nanstd(allData,1);
cellNameDataArray = vertcat(cellNameData{:});

meansQueryRunFeatureDataAllCells = nan(numel(unique(cellNameDataArray)),numel(queryOptoParameterValuesOverall));

for iQueryOptoParameter = 1:length(queryRunFeatureData)

    cellNames = unique(cellNameData{iQueryOptoParameter})
    
    for iCell = 1:numel(cellNames)
        
        cellRows = find(cellNameData{iQueryOptoParameter} == cellNames(iCell)); 
        meanQueryRunFeatureDataCell = nanmean(queryRunFeatureData{iQueryOptoParameter}(cellRows))
        
        if ~isnan(queryOptoParameterValues(iQueryOptoParameter))
        
            positionOptoParameterValue = find(queryOptoParameterValuesOverall == queryOptoParameterValues(iQueryOptoParameter));
        
        else
            
            positionOptoParameterValue = find(isnan(queryOptoParameterValuesOverall));
            
        end
        
        meansQueryRunFeatureDataAllCells(find(selectedCellsIDs == cellNames(iCell)),positionOptoParameterValue) = meanQueryRunFeatureDataCell;
        
    end
    
end

meanQueryRunFeatureDataAllCells = nanmean(meansQueryRunFeatureDataAllCells,1);
stdQueryRunFeatureDataAllCells = nanstd(meansQueryRunFeatureDataAllCells,1);

selectedCellsIDsRest = [3,4,5,6];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDsRest));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

selectedRunGroupsRest = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao2') & cell2mat(runGroupsParametersTable.holdingVoltage) == 0 & cell2mat(runGroupsParametersTable.pulseWidthRed) == 100 & strcmp(runGroupsParametersTable.functionNameRed,'dmdOutputRect(1, 1, 1, 1, 2, 1)')& cell2mat(runGroupsParametersTable.nPulsesRed) == 1);
selectedRunGroupsAnalysisTableRest = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroupsRest),:);

queryRunFeatureData = {};
cellNameData = {};

groupsCounter = 1;

for iRunGroup = 1:numel(selectedRunGroupsRest)
    
    rowsRunGroup = find(selectedRunGroupsAnalysisTableRest.runGroup(:) == selectedRunGroupsRest(iRunGroup));
    
    if ~isempty(rowsRunGroup)
    
        queryOptoParameterValues(groupsCounter) = selectedRunGroupsAnalysisTableRest.optoParameters{rowsRunGroup(1)}.(queryOptoParameter);
        queryRunFeatureData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTableRest.(queryRunFeature)(rowsRunGroup));
        cellNameData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTableRest.cellName(rowsRunGroup));

    end
    
    groupsCounter = groupsCounter + 1;
    
end

allDataRest = vertcat(queryRunFeatureData{:});
meanAllDataRest = nanmean(allDataRest,1);
stdAllDataRest = nanstd(allDataRest,1);
cellNameDataArray = vertcat(cellNameData{:});

meansQueryRunFeatureDataRestAllCells = nan(numel(unique(cellNameDataArray)),numel(queryOptoParameterValuesOverall));

for iQueryOptoParameter = 1:length(queryRunFeatureData)

    cellNames = unique(cellNameData{iQueryOptoParameter})
    
    for iCell = 1:numel(cellNames)
        
        cellRows = find(cellNameData{iQueryOptoParameter} == cellNames(iCell)); 
        meanQueryRunFeatureDataCell = nanmean(queryRunFeatureData{iQueryOptoParameter}(cellRows))
        
        if ~isnan(queryOptoParameterValues(iQueryOptoParameter))
        
            positionOptoParameterValue = find(queryOptoParameterValuesOverall == queryOptoParameterValues(iQueryOptoParameter));
        
        else
            
            positionOptoParameterValue = find(isnan(queryOptoParameterValuesOverall));
            
        end
        
        meansQueryRunFeatureDataRestAllCells(find(selectedCellsIDsRest == cellNames(iCell)),positionOptoParameterValue) = meanQueryRunFeatureDataCell;
        
    end
    
end

meanQueryRunFeatureDataRestAllCells = nanmean(meansQueryRunFeatureDataRestAllCells,1);
stdQueryRunFeatureDataRestAllCells = nanstd(meansQueryRunFeatureDataRestAllCells,1);

figure
scatter([0], meansQueryRunFeatureDataRestAllCells,'k','filled');
hold on
scatter([1], meansQueryRunFeatureDataAllCells,'k','filled');
hold on
bar([0,1], [meanQueryRunFeatureDataRestAllCells, meanQueryRunFeatureDataAllCells] ,'FaceColor',[0 0.4470 0.7410],'BarWidth', 0.5);
errorbar([0,1],[meanQueryRunFeatureDataRestAllCells, meanQueryRunFeatureDataAllCells],[stdQueryRunFeatureDataRestAllCells, stdQueryRunFeatureDataAllCells],'.','Color',[0.3 0.3 0.3],'LineWidth', 1.2, 'CapSize', 14)

hScatter = findobj(gca, 'Type', 'scatter');
uistack(hScatter, 'top');

xticks([0,1])
xticklabels({'Non-connected cells','Cells 1-2-7-8'})
xlabel(['Cells']);
ylabel(queryRunFeature);
ylim([-0.5,1])
box on;
title([queryOptoParameter,' vs ', queryRunFeature]);
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');

plotName = [experimentName '_', queryOptoParameter , '_vs_', queryRunFeature, '_', 'cellsComparison'];
plotFullPath = fullfile(savePlotPath, plotName);

if savePlots == 1
    saveas(gcf, plotFullPath, 'pdf')
    saveas(gcf, plotFullPath, 'fig')
end

%% amplitudeBlue vs heightPulsePeak

queryOptoParameter =  'amplitudeBlue';
queryRunFeature = 'heightPulsePeak';

queryOptoParameterValuesOverall = [5];
nNaNs = numel(find(isnan(queryOptoParameterValuesOverall)));

selectedCellsIDs = [3];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.holdingVoltage) == 0 & cell2mat(runGroupsParametersTable.pulseWidthBlue) == 100 & strcmp(runGroupsParametersTable.functionNameBlue,'dmdOutputRect(1, 1, 1, 1, 1, 1)') & cell2mat(runGroupsParametersTable.nPulsesBlue) == 1 & cell2mat(runGroupsParametersTable.amplitudeBlue) == 5);
selectedRunGroupsAnalysisTable = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroups),:);

% desiredEpochsDrugsRows = find(cell2mat(selectedRunGroupsAnalysisTable.epoch) >= 36);
% selectedRunGroupsAnalysisTableDrugs = selectedRunGroupsAnalysisTable(desiredEpochsDrugsRows,:)

% desiredEpochsRows = find(cell2mat(selectedRunGroupsAnalysisTable.epoch) > 3 & cell2mat(selectedRunGroupsAnalysisTable.epoch) <= 13);
% selectedRunGroupsAnalysisTable = selectedRunGroupsAnalysisTable(desiredEpochsRows,:)

queryOptoParameterValues = nan(numel(selectedRunGroups),1);
queryRunFeatureData = {};
cellNameData = {};
allTraces = [];

groupsCounter = 1;

for iRunGroup = 1:numel(selectedRunGroups)
    
    rowsRunGroup = find(selectedRunGroupsAnalysisTable.runGroup(:) == selectedRunGroups(iRunGroup));
    
    if ~isempty(rowsRunGroup)
    
        queryOptoParameterValues(groupsCounter) = selectedRunGroupsAnalysisTable.optoParameters{rowsRunGroup(1)}.(queryOptoParameter);
        queryRunFeatureData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTable.(queryRunFeature)(rowsRunGroup));
        cellNameData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTable.cellName(rowsRunGroup));

        tracesGroup = selectedRunGroupsAnalysisTable.trace(rowsRunGroup)
        allTracesGroup = [];
        
        for iTrace = 1:numel(tracesGroup)
            
            trace = tracesGroup{iTrace};
            
            if numel(trace) < 10000
                
                trace = [zeros(1,5000) trace zeros(1,5000)];
                
            end
            
            allTracesGroup = [allTracesGroup; trace];
            
        end
        
        allTraces = [allTraces; allTracesGroup];
    
    end
    
    groupsCounter = groupsCounter + 1;
    
end

allData = vertcat(queryRunFeatureData{:});
meanAllData = nanmean(allData,1);
stdAllData = nanstd(allData,1);
cellNameDataArray = vertcat(cellNameData{:});

meansQueryRunFeatureDataAllCells = nan(numel(unique(cellNameDataArray)),numel(queryOptoParameterValuesOverall));

for iQueryOptoParameter = 1:length(queryRunFeatureData)

    cellNames = unique(cellNameData{iQueryOptoParameter})
    
    for iCell = 1:numel(cellNames)
        
        cellRows = find(cellNameData{iQueryOptoParameter} == cellNames(iCell)); 
        meanQueryRunFeatureDataCell = nanmean(queryRunFeatureData{iQueryOptoParameter}(cellRows))
        
        if ~isnan(queryOptoParameterValues(iQueryOptoParameter))
        
            positionOptoParameterValue = find(queryOptoParameterValuesOverall == queryOptoParameterValues(iQueryOptoParameter));
        
        else
            
            positionOptoParameterValue = find(isnan(queryOptoParameterValuesOverall));
            
        end
        
        meansQueryRunFeatureDataAllCells(find(selectedCellsIDs == cellNames(iCell)),positionOptoParameterValue) = meanQueryRunFeatureDataCell;
        
    end
    
end

meanQueryRunFeatureDataAllCells = nanmean(meansQueryRunFeatureDataAllCells,1);
stdQueryRunFeatureDataAllCells = nanstd(meansQueryRunFeatureDataAllCells,1);

selectedCellsIDsRest = [1,2,4,5,6];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDsRest));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

selectedRunGroupsRest = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.holdingVoltage) == 0 & cell2mat(runGroupsParametersTable.pulseWidthBlue) == 100 & strcmp(runGroupsParametersTable.functionNameBlue,'dmdOutputRect(1, 1, 1, 1, 1, 1)')& cell2mat(runGroupsParametersTable.nPulsesBlue) == 1 & cell2mat(runGroupsParametersTable.amplitudeBlue) == 5);
selectedRunGroupsAnalysisTableRest = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroupsRest),:);

queryRunFeatureData = {};
cellNameData = {};
allTracesRest = [];

groupsCounter = 1;

for iRunGroup = 1:numel(selectedRunGroupsRest)
    
    rowsRunGroup = find(selectedRunGroupsAnalysisTableRest.runGroup(:) == selectedRunGroupsRest(iRunGroup));
    
    if ~isempty(rowsRunGroup)
    
        queryOptoParameterValues(groupsCounter) = selectedRunGroupsAnalysisTableRest.optoParameters{rowsRunGroup(1)}.(queryOptoParameter);
        queryRunFeatureData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTableRest.(queryRunFeature)(rowsRunGroup));
        cellNameData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTableRest.cellName(rowsRunGroup));
        tracesGroup = selectedRunGroupsAnalysisTableRest.trace(rowsRunGroup)
        allTracesGroup = [];
        
        for iTrace = 1:numel(tracesGroup)
            
            trace = tracesGroup{iTrace};
            
            if numel(trace) < 10000
                
                trace = [zeros(1,5000) trace zeros(1,5000)];
                
            end
            
            allTracesGroup = [allTracesGroup; trace];
            
        end
        
        allTracesRest = [allTracesRest; allTracesGroup];
    
    end
    
    groupsCounter = groupsCounter + 1;
    
end

allDataRest = vertcat(queryRunFeatureData{:});
meanAllDataRest = nanmean(allDataRest,1);
stdAllDataRest = nanstd(allDataRest,1);
cellNameDataArray = vertcat(cellNameData{:});

meansQueryRunFeatureDataRestAllCells = nan(numel(unique(cellNameDataArray)),numel(queryOptoParameterValuesOverall));

for iQueryOptoParameter = 1:length(queryRunFeatureData)

    cellNames = unique(cellNameData{iQueryOptoParameter})
    
    for iCell = 1:numel(cellNames)
        
        cellRows = find(cellNameData{iQueryOptoParameter} == cellNames(iCell)); 
        meanQueryRunFeatureDataCell = nanmean(queryRunFeatureData{iQueryOptoParameter}(cellRows))
        
        if ~isnan(queryOptoParameterValues(iQueryOptoParameter))
        
            positionOptoParameterValue = find(queryOptoParameterValuesOverall == queryOptoParameterValues(iQueryOptoParameter));
        
        else
            
            positionOptoParameterValue = find(isnan(queryOptoParameterValuesOverall));
            
        end
        
        meansQueryRunFeatureDataRestAllCells(find(selectedCellsIDsRest == cellNames(iCell)),positionOptoParameterValue) = meanQueryRunFeatureDataCell;
        
    end
    
end

meanQueryRunFeatureDataRestAllCells = nanmean(meansQueryRunFeatureDataRestAllCells,1);
stdQueryRunFeatureDataRestAllCells = nanstd(meansQueryRunFeatureDataRestAllCells,1);

figure
scatter([0], meansQueryRunFeatureDataRestAllCells,'k','filled');
hold on
scatter([1], meansQueryRunFeatureDataAllCells,'k','filled');
hold on
bar([0,1], [meanQueryRunFeatureDataRestAllCells, meanQueryRunFeatureDataAllCells] ,'FaceColor',[0 0.4470 0.7410],'BarWidth', 0.5);
errorbar([0,1],[meanQueryRunFeatureDataRestAllCells, meanQueryRunFeatureDataAllCells],[stdQueryRunFeatureDataRestAllCells, stdQueryRunFeatureDataAllCells],'.','Color',[0.3 0.3 0.3],'LineWidth', 1.2, 'CapSize', 14)

hScatter = findobj(gca, 'Type', 'scatter');
uistack(hScatter, 'top');

xticks([0,1])
xticklabels({'Non-connected cells','Cells 1-2-7-8'})
xlabel(['Cells']);
ylabel(queryRunFeature);
ylim([-10,70])
box on;
title([queryOptoParameter,' vs ', queryRunFeature]);
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');

plotName = [experimentName '_', queryOptoParameter , '_vs_', queryRunFeature, '_', 'cellsComparison'];
plotFullPath = fullfile(savePlotPath, plotName);

if savePlots == 1
    saveas(gcf, plotFullPath, 'pdf')
    saveas(gcf, plotFullPath, 'fig')
end

figure;
plot(preprocessSignalVCv2(nanmean(allTraces,1)),'Color',[0 0.4470 0.7410],'LineWidth',2)
hold on

errorUp = nanmean(allTraces,1) + nanstd(allTraces,1);
errorDown = nanmean(allTraces,1) - nanstd(allTraces,1);
xLinspace = linspace(1,size(allTraces,2), size(allTraces,2));
% plotpatch = patch([xLinspace fliplr(xLinspace)], [errorDown fliplr(errorUp)], [.63, .79, .95], 'FaceAlpha',0.5, 'EdgeColor','none');
% plotpatch.Annotation.LegendInformation.IconDisplayStyle = 'off';

plot(preprocessSignalVCv2(nanmean(allTracesRest,1)),'Color',[0.4 0.4 0.4],'LineWidth',2)

errorUp = nanmean(allTracesRest,1) + nanstd(allTracesRest,1);
errorDown = nanmean(allTracesRest,1) - nanstd(allTracesRest,1);
xLinspace = linspace(1,size(allTracesRest,2), size(allTracesRest,2));
% plotpatch = patch([xLinspace fliplr(xLinspace)], [errorDown fliplr(errorUp)], [.8, .8, .8], 'FaceAlpha',0.5, 'EdgeColor','none');
% plotpatch.Annotation.LegendInformation.IconDisplayStyle = 'off';


xlim([4500,6500])
ylim([-10,40])
xticks([4500,5500,6500])
xticklabels({'0','0.1','0.2'})
xlabel('Time [s]')
ylabel('Current [pA]')
box on;
title([queryOptoParameter,' vs ', queryRunFeature]);
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');
legend('Connected cell','Non-connected cells','Location', 'northeast', 'FontSize',12)

plotName = [experimentName '_', queryOptoParameter , '_vs_', queryRunFeature, '_', 'cellsComparison_', 'traces'];
plotFullPath = fullfile(savePlotPath, plotName);

if savePlots == 1
    saveas(gcf, plotFullPath, 'pdf')
    saveas(gcf, plotFullPath, 'fig')
end

%% amplitudeBlue vs areaPulse

queryOptoParameter =  'amplitudeBlue';
queryRunFeature = 'areaPulse';

queryOptoParameterValuesOverall = [5];
nNaNs = numel(find(isnan(queryOptoParameterValuesOverall)));

selectedCellsIDs = [3];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.holdingVoltage) == 0 & cell2mat(runGroupsParametersTable.pulseWidthBlue) == 100 & strcmp(runGroupsParametersTable.functionNameBlue,'dmdOutputRect(1, 1, 1, 1, 1, 1)')& cell2mat(runGroupsParametersTable.nPulsesBlue) == 1 & cell2mat(runGroupsParametersTable.amplitudeBlue) == 5);
selectedRunGroupsAnalysisTable = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroups),:);

% desiredEpochsDrugsRows = find(cell2mat(selectedRunGroupsAnalysisTable.epoch) >= 36);
% selectedRunGroupsAnalysisTableDrugs = selectedRunGroupsAnalysisTable(desiredEpochsDrugsRows,:)

% desiredEpochsRows = find(cell2mat(selectedRunGroupsAnalysisTable.epoch) > 3 & cell2mat(selectedRunGroupsAnalysisTable.epoch) <= 13);
% selectedRunGroupsAnalysisTable = selectedRunGroupsAnalysisTable(desiredEpochsRows,:)

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

allData = vertcat(queryRunFeatureData{:});
meanAllData = nanmean(allData,1);
stdAllData = nanstd(allData,1);
cellNameDataArray = vertcat(cellNameData{:});

meansQueryRunFeatureDataAllCells = nan(numel(unique(cellNameDataArray)),numel(queryOptoParameterValuesOverall));

for iQueryOptoParameter = 1:length(queryRunFeatureData)

    cellNames = unique(cellNameData{iQueryOptoParameter})
    
    for iCell = 1:numel(cellNames)
        
        cellRows = find(cellNameData{iQueryOptoParameter} == cellNames(iCell)); 
        meanQueryRunFeatureDataCell = nanmean(queryRunFeatureData{iQueryOptoParameter}(cellRows))
        
        if ~isnan(queryOptoParameterValues(iQueryOptoParameter))
        
            positionOptoParameterValue = find(queryOptoParameterValuesOverall == queryOptoParameterValues(iQueryOptoParameter));
        
        else
            
            positionOptoParameterValue = find(isnan(queryOptoParameterValuesOverall));
            
        end
        
        meansQueryRunFeatureDataAllCells(find(selectedCellsIDs == cellNames(iCell)),positionOptoParameterValue) = meanQueryRunFeatureDataCell;
        
    end
    
end

meanQueryRunFeatureDataAllCells = nanmean(meansQueryRunFeatureDataAllCells,1);
stdQueryRunFeatureDataAllCells = nanstd(meansQueryRunFeatureDataAllCells,1);

selectedCellsIDsRest = [1,2,4,5,6];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDsRest));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

selectedRunGroupsRest = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.holdingVoltage) == 0 & cell2mat(runGroupsParametersTable.pulseWidthBlue) == 100 & strcmp(runGroupsParametersTable.functionNameBlue,'dmdOutputRect(1, 1, 1, 1, 1, 1)')& cell2mat(runGroupsParametersTable.nPulsesBlue) == 1 & cell2mat(runGroupsParametersTable.amplitudeBlue) == 5);
selectedRunGroupsAnalysisTableRest = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroupsRest),:);

queryRunFeatureData = {};
cellNameData = {};

groupsCounter = 1;

for iRunGroup = 1:numel(selectedRunGroupsRest)
    
    rowsRunGroup = find(selectedRunGroupsAnalysisTableRest.runGroup(:) == selectedRunGroupsRest(iRunGroup));
    
    if ~isempty(rowsRunGroup)
    
        queryOptoParameterValues(groupsCounter) = selectedRunGroupsAnalysisTableRest.optoParameters{rowsRunGroup(1)}.(queryOptoParameter);
        queryRunFeatureData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTableRest.(queryRunFeature)(rowsRunGroup));
        cellNameData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTableRest.cellName(rowsRunGroup));

    end
    
    groupsCounter = groupsCounter + 1;
    
end

allDataRest = vertcat(queryRunFeatureData{:});
meanAllDataRest = nanmean(allDataRest,1);
stdAllDataRest = nanstd(allDataRest,1);
cellNameDataArray = vertcat(cellNameData{:});

meansQueryRunFeatureDataRestAllCells = nan(numel(unique(cellNameDataArray)),numel(queryOptoParameterValuesOverall));

for iQueryOptoParameter = 1:length(queryRunFeatureData)

    cellNames = unique(cellNameData{iQueryOptoParameter})
    
    for iCell = 1:numel(cellNames)
        
        cellRows = find(cellNameData{iQueryOptoParameter} == cellNames(iCell)); 
        meanQueryRunFeatureDataCell = nanmean(queryRunFeatureData{iQueryOptoParameter}(cellRows))
        
        if ~isnan(queryOptoParameterValues(iQueryOptoParameter))
        
            positionOptoParameterValue = find(queryOptoParameterValuesOverall == queryOptoParameterValues(iQueryOptoParameter));
        
        else
            
            positionOptoParameterValue = find(isnan(queryOptoParameterValuesOverall));
            
        end
        
        meansQueryRunFeatureDataRestAllCells(find(selectedCellsIDsRest == cellNames(iCell)),positionOptoParameterValue) = meanQueryRunFeatureDataCell;
        
    end
    
end

meanQueryRunFeatureDataRestAllCells = nanmean(meansQueryRunFeatureDataRestAllCells,1);
stdQueryRunFeatureDataRestAllCells = nanstd(meansQueryRunFeatureDataRestAllCells,1);

figure
scatter([0], meansQueryRunFeatureDataRestAllCells,'k','filled');
hold on
scatter([1], meansQueryRunFeatureDataAllCells,'k','filled');
hold on
bar([0,1], [meanQueryRunFeatureDataRestAllCells, meanQueryRunFeatureDataAllCells] ,'FaceColor',[0 0.4470 0.7410],'BarWidth', 0.5);
errorbar([0,1],[meanQueryRunFeatureDataRestAllCells, meanQueryRunFeatureDataAllCells],[stdQueryRunFeatureDataRestAllCells, stdQueryRunFeatureDataAllCells],'.','Color',[0.3 0.3 0.3],'LineWidth', 1.2, 'CapSize', 14)

hScatter = findobj(gca, 'Type', 'scatter');
uistack(hScatter, 'top');

xticks([0,1])
xticklabels({'Non-connected cells','Connected cells'})
xlabel(['Cells']);
ylabel(queryRunFeature);
% ylim([-0.5,1])
box on;
title([queryOptoParameter,' vs ', queryRunFeature]);
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');

plotName = [experimentName '_', queryOptoParameter , '_vs_', queryRunFeature, '_', 'cellsComparison'];
plotFullPath = fullfile(savePlotPath, plotName);

if savePlots == 1
    saveas(gcf, plotFullPath, 'pdf')
    saveas(gcf, plotFullPath, 'fig')
end