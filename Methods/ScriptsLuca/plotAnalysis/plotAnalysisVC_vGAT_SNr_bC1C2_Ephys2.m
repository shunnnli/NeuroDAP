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
experimentName = 'vGAT_SNr_bC1C2_Ephys_2';

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

%% amplitudeBlue vs heightPulsePeak

queryOptoParameter =  'amplitudeBlue';
queryRunFeature = 'heightPulsePeak';

queryOptoParameterValuesOverall = [5];
nNaNs = numel(find(isnan(queryOptoParameterValuesOverall)));

selectedCellsIDs = [2,5];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.holdingVoltage) == 0 & cell2mat(runGroupsParametersTable.pulseWidthBlue) == 100 & cell2mat(runGroupsParametersTable.nPulsesBlue) == 1);
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

selectedCellsIDsRest = [4];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDsRest));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

selectedRunGroupsRest = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.holdingVoltage) == 0 & cell2mat(runGroupsParametersTable.pulseWidthBlue) == 100 & cell2mat(runGroupsParametersTable.nPulsesBlue) == 1);
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

selectedCellsIDsWeak = [1,3];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDsWeak));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

selectedRunGroupsWeak = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.holdingVoltage) == 0 & cell2mat(runGroupsParametersTable.pulseWidthBlue) == 100 & cell2mat(runGroupsParametersTable.nPulsesBlue) == 1);
selectedRunGroupsAnalysisTableWeak = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroupsWeak),:);

queryRunFeatureData = {};
cellNameData = {};
allTracesWeak = [];

groupsCounter = 1;

for iRunGroup = 1:numel(selectedRunGroupsWeak)
    
    rowsRunGroup = find(selectedRunGroupsAnalysisTableWeak.runGroup(:) == selectedRunGroupsWeak(iRunGroup));
    
    if ~isempty(rowsRunGroup)
    
        queryOptoParameterValues(groupsCounter) = selectedRunGroupsAnalysisTableWeak.optoParameters{rowsRunGroup(1)}.(queryOptoParameter);
        queryRunFeatureData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTableWeak.(queryRunFeature)(rowsRunGroup));
        cellNameData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTableWeak.cellName(rowsRunGroup));
        tracesGroup = selectedRunGroupsAnalysisTableWeak.trace(rowsRunGroup)
        allTracesGroup = [];
        
        for iTrace = 1:numel(tracesGroup)
            
            trace = tracesGroup{iTrace};
            
            if numel(trace) < 10000
                
                trace = [zeros(1,5000) trace zeros(1,5000)];
                
            end
            
            allTracesGroup = [allTracesGroup; trace];
            
        end
        
        allTracesWeak = [allTracesWeak; allTracesGroup];
    
    end
    
    groupsCounter = groupsCounter + 1;
    
end

allDataWeak = vertcat(queryRunFeatureData{:});
meanAllDataWeak = nanmean(allDataWeak,1);
stdAllDataWeak = nanstd(allDataWeak,1);
cellNameDataArray = vertcat(cellNameData{:});

meansQueryRunFeatureDataWeakAllCells = nan(numel(unique(cellNameDataArray)),numel(queryOptoParameterValuesOverall));

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
        
        meansQueryRunFeatureDataWeakAllCells(find(selectedCellsIDsWeak == cellNames(iCell)),positionOptoParameterValue) = meanQueryRunFeatureDataCell;
        
    end
    
end

meanQueryRunFeatureDataWeakAllCells = nanmean(meansQueryRunFeatureDataWeakAllCells,1);
stdQueryRunFeatureDataWeakAllCells = nanstd(meansQueryRunFeatureDataWeakAllCells,1);


figure
scatter([0], meansQueryRunFeatureDataRestAllCells,'k','filled');
hold on
scatter([1], meansQueryRunFeatureDataWeakAllCells,'k','filled');
hold on
scatter([2], meansQueryRunFeatureDataAllCells,'k','filled');
hold on
bar([0,1,2], [meanQueryRunFeatureDataRestAllCells, meanQueryRunFeatureDataWeakAllCells, meanQueryRunFeatureDataAllCells] ,'FaceColor',[0 0.4470 0.7410],'BarWidth', 0.5);
errorbar([0,1,2],[meanQueryRunFeatureDataRestAllCells, meanQueryRunFeatureDataWeakAllCells, meanQueryRunFeatureDataAllCells],[stdQueryRunFeatureDataRestAllCells, stdQueryRunFeatureDataWeakAllCells, stdQueryRunFeatureDataAllCells],'.','Color',[0.3 0.3 0.3],'LineWidth', 1.2, 'CapSize', 14)

hScatter = findobj(gca, 'Type', 'scatter');
uistack(hScatter, 'top');

xticks([0,1,2])
xticklabels({'Non-connected','Weakly connected','Strongly connected'})
xlabel(['Cells']);
ylabel(queryRunFeature);
ylim([0,600])
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
plot(nanmean(allTraces,1),'Color',[0 0.4470 0.7410],'LineWidth',2)
hold on

errorUp = nanmean(allTraces,1) + nanstd(allTraces,1);
errorDown = nanmean(allTraces,1) - nanstd(allTraces,1);
xLinspace = linspace(1,size(allTraces,2), size(allTraces,2));
%plotpatch = patch([xLinspace fliplr(xLinspace)], [errorDown fliplr(errorUp)], [.6, .2, .2], 'FaceAlpha',0.5, 'EdgeColor','none');
%plotpatch.Annotation.LegendInformation.IconDisplayStyle = 'off';

plot(nanmean(allTracesWeak,1),'Color',[0 0.7 0.9],'LineWidth',2)
hold on

errorUp = nanmean(allTracesWeak,1) + nanstd(allTracesWeak,1);
errorDown = nanmean(allTracesWeak,1) - nanstd(allTracesWeak,1);
xLinspace = linspace(1,size(allTracesWeak,2), size(allTracesWeak,2));
%plotpatch = patch([xLinspace fliplr(xLinspace)], [errorDown fliplr(errorUp)], [.6, .2, .2], 'FaceAlpha',0.5, 'EdgeColor','none');
%plotpatch.Annotation.LegendInformation.IconDisplayStyle = 'off';

plot(nanmean(allTracesRest,1),'Color',[0.4 0.4 0.4],'LineWidth',2)

errorUp = nanmean(allTracesRest,1) + nanstd(allTracesRest,1);
errorDown = nanmean(allTracesRest,1) - nanstd(allTracesRest,1);
xLinspace = linspace(1,size(allTracesRest,2), size(allTracesRest,2));
%plotpatch = patch([xLinspace fliplr(xLinspace)], [errorDown fliplr(errorUp)], [.8, .8, .8], 'FaceAlpha',0.5, 'EdgeColor','none');
%plotpatch.Annotation.LegendInformation.IconDisplayStyle = 'off';

xlim([4000,7000])
ylim([-50,350])
%xticks([4500,5500,6500])
xticklabels({'0','0.1','0.2'})
xlabel('Time [s]')
ylabel('Current [pA]')
box on;
title([queryOptoParameter,' vs ', queryRunFeature]);
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');
legend('Strongly connected','Weakly connected', 'Non-connected cells','Location', 'northeast', 'FontSize',12)

plotName = [experimentName '_', queryOptoParameter , '_vs_', queryRunFeature, '_', 'cellsComparison_', 'traces'];
plotFullPath = fullfile(savePlotPath, plotName);

if savePlots == 1
    saveas(gcf, plotFullPath, 'pdf')
    saveas(gcf, plotFullPath, 'fig')
end

meansStrongBlueHeightPulsePeak  = meansQueryRunFeatureDataAllCells;
tracesStrongBlue = allTraces;

%% amplitudeBlue vs areaPulse

queryOptoParameter =  'amplitudeBlue';
queryRunFeature = 'areaPulse';

queryOptoParameterValuesOverall = [5];
nNaNs = numel(find(isnan(queryOptoParameterValuesOverall)));

selectedCellsIDs = [2,5];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.holdingVoltage) == 0 & cell2mat(runGroupsParametersTable.pulseWidthBlue) == 100 & cell2mat(runGroupsParametersTable.nPulsesBlue) == 1);
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

selectedCellsIDsRest = [4];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDsRest));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

selectedRunGroupsRest = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.holdingVoltage) == 0 & cell2mat(runGroupsParametersTable.pulseWidthBlue) == 100 & cell2mat(runGroupsParametersTable.nPulsesBlue) == 1);
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

selectedCellsIDsWeak = [1,3];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDsWeak));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

selectedRunGroupsWeak = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.holdingVoltage) == 0 & cell2mat(runGroupsParametersTable.pulseWidthBlue) == 100 & cell2mat(runGroupsParametersTable.nPulsesBlue) == 1);
selectedRunGroupsAnalysisTableWeak = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroupsWeak),:);

queryRunFeatureData = {};
cellNameData = {};
allTracesWeak = [];

groupsCounter = 1;

for iRunGroup = 1:numel(selectedRunGroupsWeak)
    
    rowsRunGroup = find(selectedRunGroupsAnalysisTableWeak.runGroup(:) == selectedRunGroupsWeak(iRunGroup));
    
    if ~isempty(rowsRunGroup)
    
        queryOptoParameterValues(groupsCounter) = selectedRunGroupsAnalysisTableWeak.optoParameters{rowsRunGroup(1)}.(queryOptoParameter);
        queryRunFeatureData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTableWeak.(queryRunFeature)(rowsRunGroup));
        cellNameData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTableWeak.cellName(rowsRunGroup));
        tracesGroup = selectedRunGroupsAnalysisTableWeak.trace(rowsRunGroup)
        allTracesGroup = [];
        
        for iTrace = 1:numel(tracesGroup)
            
            trace = tracesGroup{iTrace};
            
            if numel(trace) < 10000
                
                trace = [zeros(1,5000) trace zeros(1,5000)];
                
            end
            
            allTracesGroup = [allTracesGroup; trace];
            
        end
        
        allTracesWeak = [allTracesWeak; allTracesGroup];
    
    end
    
    groupsCounter = groupsCounter + 1;
    
end

allDataWeak = vertcat(queryRunFeatureData{:});
meanAllDataWeak = nanmean(allDataWeak,1);
stdAllDataWeak = nanstd(allDataWeak,1);
cellNameDataArray = vertcat(cellNameData{:});

meansQueryRunFeatureDataWeakAllCells = nan(numel(unique(cellNameDataArray)),numel(queryOptoParameterValuesOverall));

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
        
        meansQueryRunFeatureDataWeakAllCells(find(selectedCellsIDsWeak == cellNames(iCell)),positionOptoParameterValue) = meanQueryRunFeatureDataCell;
        
    end
    
end

meanQueryRunFeatureDataWeakAllCells = nanmean(meansQueryRunFeatureDataWeakAllCells,1);
stdQueryRunFeatureDataWeakAllCells = nanstd(meansQueryRunFeatureDataWeakAllCells,1);


figure
scatter([0], meansQueryRunFeatureDataRestAllCells,'k','filled');
hold on
scatter([1], meansQueryRunFeatureDataWeakAllCells,'k','filled');
hold on
scatter([2], meansQueryRunFeatureDataAllCells,'k','filled');
hold on
bar([0,1,2], [meanQueryRunFeatureDataRestAllCells, meanQueryRunFeatureDataWeakAllCells, meanQueryRunFeatureDataAllCells] ,'FaceColor',[0 0.4470 0.7410],'BarWidth', 0.5);
errorbar([0,1,2],[meanQueryRunFeatureDataRestAllCells, meanQueryRunFeatureDataWeakAllCells, meanQueryRunFeatureDataAllCells],[stdQueryRunFeatureDataRestAllCells, stdQueryRunFeatureDataWeakAllCells, stdQueryRunFeatureDataAllCells],'.','Color',[0.3 0.3 0.3],'LineWidth', 1.2, 'CapSize', 14)

hScatter = findobj(gca, 'Type', 'scatter');
uistack(hScatter, 'top');

xticks([0,1,2])
xticklabels({'Non-connected','Weakly connected','Strongly connected'})
xlabel(['Cells']);
ylabel(queryRunFeature);
ylim([-2,14])
box on;
title([queryOptoParameter,' vs ', queryRunFeature]);
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');

plotName = [experimentName '_', queryOptoParameter , '_vs_', queryRunFeature, '_', 'cellsComparison'];
plotFullPath = fullfile(savePlotPath, plotName);

if savePlots == 1
    saveas(gcf, plotFullPath, 'pdf')
    saveas(gcf, plotFullPath, 'fig')
end

meansStrongBlueAreaPulse  = meansQueryRunFeatureDataAllCells;

%% amplitudeRed vs heightPulsePeak

queryOptoParameter =  'amplitudeRed';
queryRunFeature = 'heightPulsePeak';

queryOptoParameterValuesOverall = [5];
nNaNs = numel(find(isnan(queryOptoParameterValuesOverall)));

selectedCellsIDs = [2,5];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao2') & cell2mat(runGroupsParametersTable.holdingVoltage) == 0 & cell2mat(runGroupsParametersTable.pulseWidthRed) == 100 & cell2mat(runGroupsParametersTable.nPulsesRed) == 1);
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

selectedCellsIDsRest = [4];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDsRest));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

selectedRunGroupsRest = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao2') & cell2mat(runGroupsParametersTable.holdingVoltage) == 0 & cell2mat(runGroupsParametersTable.pulseWidthRed) == 100 & cell2mat(runGroupsParametersTable.nPulsesRed) == 1);
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

selectedCellsIDsWeak = [1,3];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDsWeak));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

selectedRunGroupsWeak = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao2') & cell2mat(runGroupsParametersTable.holdingVoltage) == 0 & cell2mat(runGroupsParametersTable.pulseWidthRed) == 100 & cell2mat(runGroupsParametersTable.nPulsesRed) == 1);
selectedRunGroupsAnalysisTableWeak = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroupsWeak),:);

queryRunFeatureData = {};
cellNameData = {};
allTracesWeak = [];

groupsCounter = 1;

for iRunGroup = 1:numel(selectedRunGroupsWeak)
    
    rowsRunGroup = find(selectedRunGroupsAnalysisTableWeak.runGroup(:) == selectedRunGroupsWeak(iRunGroup));
    
    if ~isempty(rowsRunGroup)
    
        queryOptoParameterValues(groupsCounter) = selectedRunGroupsAnalysisTableWeak.optoParameters{rowsRunGroup(1)}.(queryOptoParameter);
        queryRunFeatureData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTableWeak.(queryRunFeature)(rowsRunGroup));
        cellNameData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTableWeak.cellName(rowsRunGroup));
        tracesGroup = selectedRunGroupsAnalysisTableWeak.trace(rowsRunGroup)
        allTracesGroup = [];
        
        for iTrace = 1:numel(tracesGroup)
            
            trace = tracesGroup{iTrace};
            
            if numel(trace) < 10000
                
                trace = [zeros(1,5000) trace zeros(1,5000)];
                
            end
            
            allTracesGroup = [allTracesGroup; trace];
            
        end
        
        allTracesWeak = [allTracesWeak; allTracesGroup];
    
    end
    
    groupsCounter = groupsCounter + 1;
    
end

allDataWeak = vertcat(queryRunFeatureData{:});
meanAllDataWeak = nanmean(allDataWeak,1);
stdAllDataWeak = nanstd(allDataWeak,1);
cellNameDataArray = vertcat(cellNameData{:});

meansQueryRunFeatureDataWeakAllCells = nan(numel(unique(cellNameDataArray)),numel(queryOptoParameterValuesOverall));

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
        
        meansQueryRunFeatureDataWeakAllCells(find(selectedCellsIDsWeak == cellNames(iCell)),positionOptoParameterValue) = meanQueryRunFeatureDataCell;
        
    end
    
end

meanQueryRunFeatureDataWeakAllCells = nanmean(meansQueryRunFeatureDataWeakAllCells,1);
stdQueryRunFeatureDataWeakAllCells = nanstd(meansQueryRunFeatureDataWeakAllCells,1);


figure
scatter([0], meansQueryRunFeatureDataRestAllCells,'k','filled');
hold on
scatter([1], meansQueryRunFeatureDataWeakAllCells,'k','filled');
hold on
scatter([2], meansQueryRunFeatureDataAllCells,'k','filled');
hold on
bar([0,1,2], [meanQueryRunFeatureDataRestAllCells, meanQueryRunFeatureDataWeakAllCells, meanQueryRunFeatureDataAllCells] ,'FaceColor',[0 0.4470 0.7410],'BarWidth', 0.5);
errorbar([0,1,2],[meanQueryRunFeatureDataRestAllCells, meanQueryRunFeatureDataWeakAllCells, meanQueryRunFeatureDataAllCells],[stdQueryRunFeatureDataRestAllCells, stdQueryRunFeatureDataWeakAllCells, stdQueryRunFeatureDataAllCells],'.','Color',[0.3 0.3 0.3],'LineWidth', 1.2, 'CapSize', 14)

hScatter = findobj(gca, 'Type', 'scatter');
uistack(hScatter, 'top');

xticks([0,1,2])
xticklabels({'Non-connected','Weakly connected','Strongly connected'})
xlabel(['Cells']);
ylabel(queryRunFeature);
ylim([0,600])
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
plot(nanmean(allTraces,1),'Color',[0 0.4470 0.7410],'LineWidth',2)
hold on

errorUp = nanmean(allTraces,1) + nanstd(allTraces,1);
errorDown = nanmean(allTraces,1) - nanstd(allTraces,1);
xLinspace = linspace(1,size(allTraces,2), size(allTraces,2));
%plotpatch = patch([xLinspace fliplr(xLinspace)], [errorDown fliplr(errorUp)], [.6, .2, .2], 'FaceAlpha',0.5, 'EdgeColor','none');
%plotpatch.Annotation.LegendInformation.IconDisplayStyle = 'off';

plot(nanmean(allTracesWeak,1),'Color',[0 0.7 0.9],'LineWidth',2)
hold on

errorUp = nanmean(allTracesWeak,1) + nanstd(allTracesWeak,1);
errorDown = nanmean(allTracesWeak,1) - nanstd(allTracesWeak,1);
xLinspace = linspace(1,size(allTracesWeak,2), size(allTracesWeak,2));
%plotpatch = patch([xLinspace fliplr(xLinspace)], [errorDown fliplr(errorUp)], [.6, .2, .2], 'FaceAlpha',0.5, 'EdgeColor','none');
%plotpatch.Annotation.LegendInformation.IconDisplayStyle = 'off';

plot(nanmean(allTracesRest,1),'Color',[0.4 0.4 0.4],'LineWidth',2)

errorUp = nanmean(allTracesRest,1) + nanstd(allTracesRest,1);
errorDown = nanmean(allTracesRest,1) - nanstd(allTracesRest,1);
xLinspace = linspace(1,size(allTracesRest,2), size(allTracesRest,2));
%plotpatch = patch([xLinspace fliplr(xLinspace)], [errorDown fliplr(errorUp)], [.8, .8, .8], 'FaceAlpha',0.5, 'EdgeColor','none');
%plotpatch.Annotation.LegendInformation.IconDisplayStyle = 'off';

xlim([4500,6500])
ylim([-50,350])
xticks([4500,5500,6500])
xticklabels({'0','0.1','0.2'})
xlabel('Time [s]')
ylabel('Current [pA]')
box on;
title([queryOptoParameter,' vs ', queryRunFeature]);
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');
legend('Strongly connected','Weakly connected', 'Non-connected cells','Location', 'northeast', 'FontSize',12)

plotName = [experimentName '_', queryOptoParameter , '_vs_', queryRunFeature, '_', 'cellsComparison_', 'traces'];
plotFullPath = fullfile(savePlotPath, plotName);

if savePlots == 1
    saveas(gcf, plotFullPath, 'pdf')
    saveas(gcf, plotFullPath, 'fig')
end

meansStrongRedHeightPulsePeak  = meansQueryRunFeatureDataAllCells;
tracesStrongRed = allTraces;

%% amplitudeRed vs areaPulse

queryOptoParameter =  'amplitudeRed';
queryRunFeature = 'areaPulse';

queryOptoParameterValuesOverall = [5];
nNaNs = numel(find(isnan(queryOptoParameterValuesOverall)));

selectedCellsIDs = [2,5];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao2') & cell2mat(runGroupsParametersTable.holdingVoltage) == 0 & cell2mat(runGroupsParametersTable.pulseWidthRed) == 100 & cell2mat(runGroupsParametersTable.nPulsesRed) == 1);
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

selectedCellsIDsRest = [4];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDsRest));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

selectedRunGroupsRest = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao2') & cell2mat(runGroupsParametersTable.holdingVoltage) == 0 & cell2mat(runGroupsParametersTable.pulseWidthRed) == 100 & cell2mat(runGroupsParametersTable.nPulsesRed) == 1);
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

selectedCellsIDsWeak = [1,3];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDsWeak));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

selectedRunGroupsWeak = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao2') & cell2mat(runGroupsParametersTable.holdingVoltage) == 0 & cell2mat(runGroupsParametersTable.pulseWidthRed) == 100 & cell2mat(runGroupsParametersTable.nPulsesRed) == 1);
selectedRunGroupsAnalysisTableWeak = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroupsWeak),:);

queryRunFeatureData = {};
cellNameData = {};
allTracesWeak = [];

groupsCounter = 1;

for iRunGroup = 1:numel(selectedRunGroupsWeak)
    
    rowsRunGroup = find(selectedRunGroupsAnalysisTableWeak.runGroup(:) == selectedRunGroupsWeak(iRunGroup));
    
    if ~isempty(rowsRunGroup)
    
        queryOptoParameterValues(groupsCounter) = selectedRunGroupsAnalysisTableWeak.optoParameters{rowsRunGroup(1)}.(queryOptoParameter);
        queryRunFeatureData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTableWeak.(queryRunFeature)(rowsRunGroup));
        cellNameData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTableWeak.cellName(rowsRunGroup));
        tracesGroup = selectedRunGroupsAnalysisTableWeak.trace(rowsRunGroup)
        allTracesGroup = [];
        
        for iTrace = 1:numel(tracesGroup)
            
            trace = tracesGroup{iTrace};
            
            if numel(trace) < 10000
                
                trace = [zeros(1,5000) trace zeros(1,5000)];
                
            end
            
            allTracesGroup = [allTracesGroup; trace];
            
        end
        
        allTracesWeak = [allTracesWeak; allTracesGroup];
    
    end
    
    groupsCounter = groupsCounter + 1;
    
end

allDataWeak = vertcat(queryRunFeatureData{:});
meanAllDataWeak = nanmean(allDataWeak,1);
stdAllDataWeak = nanstd(allDataWeak,1);
cellNameDataArray = vertcat(cellNameData{:});

meansQueryRunFeatureDataWeakAllCells = nan(numel(unique(cellNameDataArray)),numel(queryOptoParameterValuesOverall));

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
        
        meansQueryRunFeatureDataWeakAllCells(find(selectedCellsIDsWeak == cellNames(iCell)),positionOptoParameterValue) = meanQueryRunFeatureDataCell;
        
    end
    
end

meanQueryRunFeatureDataWeakAllCells = nanmean(meansQueryRunFeatureDataWeakAllCells,1);
stdQueryRunFeatureDataWeakAllCells = nanstd(meansQueryRunFeatureDataWeakAllCells,1);

figure
scatter([0], meansQueryRunFeatureDataRestAllCells,'k','filled');
hold on
scatter([1], meansQueryRunFeatureDataWeakAllCells,'k','filled');
hold on
scatter([2], meansQueryRunFeatureDataAllCells,'k','filled');
hold on
bar([0,1,2], [meanQueryRunFeatureDataRestAllCells, meanQueryRunFeatureDataWeakAllCells, meanQueryRunFeatureDataAllCells] ,'FaceColor',[0 0.4470 0.7410],'BarWidth', 0.5);
errorbar([0,1,2],[meanQueryRunFeatureDataRestAllCells, meanQueryRunFeatureDataWeakAllCells, meanQueryRunFeatureDataAllCells],[stdQueryRunFeatureDataRestAllCells, stdQueryRunFeatureDataWeakAllCells, stdQueryRunFeatureDataAllCells],'.','Color',[0.3 0.3 0.3],'LineWidth', 1.2, 'CapSize', 14)

hScatter = findobj(gca, 'Type', 'scatter');
uistack(hScatter, 'top');

xticks([0,1,2])
xticklabels({'Non-connected','Weakly connected','Strongly connected'})
xlabel(['Cells']);
ylabel(queryRunFeature);
ylim([-2,14])
box on;
title([queryOptoParameter,' vs ', queryRunFeature]);
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');

plotName = [experimentName '_', queryOptoParameter , '_vs_', queryRunFeature, '_', 'cellsComparison'];
plotFullPath = fullfile(savePlotPath, plotName);

if savePlots == 1
    saveas(gcf, plotFullPath, 'pdf')
    saveas(gcf, plotFullPath, 'fig')
end

meansStrongRedAreaPulse  = meansQueryRunFeatureDataAllCells;

%%

queryOptoParameter =  'amplitude';
queryRunFeature = 'heightPulsePeak';

meanStrongBlueHeightPulsePeak = nanmean(meansStrongBlueHeightPulsePeak,1);
meanStrongRedHeightPulsePeak = nanmean(meansStrongRedHeightPulsePeak,1);
stdStrongBlueHeightPulsePeak = nanstd(meansStrongBlueHeightPulsePeak,1);
stdStrongRedHeightPulsePeak = nanstd(meansStrongRedHeightPulsePeak,1);

figure
scatter([0], meansStrongBlueHeightPulsePeak,'k','filled');
hold on
scatter([1], meansStrongRedHeightPulsePeak,'k','filled');
hold on
bar([0,1], [meanStrongBlueHeightPulsePeak, meanStrongRedHeightPulsePeak] ,'FaceColor',[0 0.4470 0.7410],'BarWidth', 0.5);
errorbar([0,1],[meanStrongBlueHeightPulsePeak, meanStrongRedHeightPulsePeak],[stdStrongBlueHeightPulsePeak, stdStrongRedHeightPulsePeak],'.','Color',[0.3 0.3 0.3],'LineWidth', 1.2, 'CapSize', 14)

hScatter = findobj(gca, 'Type', 'scatter');
uistack(hScatter, 'top');

xticks([0,1])
xticklabels({'Blue LED','Red LED'})
xlabel(['Activation channel']);
ylabel(queryRunFeature);
ylim([0,600])
box on;
title([queryOptoParameter,' vs ', queryRunFeature]);
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');

plotName = [experimentName '_', queryOptoParameter , '_vs_', queryRunFeature, '_', 'colorComparison'];
plotFullPath = fullfile(savePlotPath, plotName);

if savePlots == 1
    saveas(gcf, plotFullPath, 'pdf')
    saveas(gcf, plotFullPath, 'fig')
end

figure;
plot(nanmean(tracesStrongBlue,1),'Color',[0 0.4470 0.7410],'LineWidth',2)
hold on

errorUp = nanmean(tracesStrongBlue,1) + nanstd(tracesStrongBlue,1);
errorDown = nanmean(tracesStrongBlue,1) - nanstd(tracesStrongBlue,1);
xLinspace = linspace(1,size(tracesStrongBlue,2), size(tracesStrongBlue,2));
%plotpatch = patch([xLinspace fliplr(xLinspace)], [errorDown fliplr(errorUp)], [.6, .2, .2], 'FaceAlpha',0.5, 'EdgeColor','none');
%plotpatch.Annotation.LegendInformation.IconDisplayStyle = 'off';

plot(nanmean(tracesStrongRed,1),'Color',[.6, .1, .1],'LineWidth',2)
hold on

errorUp = nanmean(tracesStrongRed,1) + nanstd(tracesStrongRed,1);
errorDown = nanmean(tracesStrongRed,1) - nanstd(tracesStrongRed,1);
xLinspace = linspace(1,size(tracesStrongRed,2), size(tracesStrongRed,2));
%plotpatch = patch([xLinspace fliplr(xLinspace)], [errorDown fliplr(errorUp)], [.6, .2, .2], 'FaceAlpha',0.5, 'EdgeColor','none');
%plotpatch.Annotation.LegendInformation.IconDisplayStyle = 'off';


xlim([4500,6500])
ylim([-50,350])
xticks([4500,5500,6500])
xticklabels({'0','0.1','0.2'})
xlabel('Time [s]')
ylabel('Current [pA]')
box on;
title([queryOptoParameter,' vs ', queryRunFeature]);
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');
legend('Blue LED','Red LED','Location', 'northeast', 'FontSize',12)

plotName = [experimentName '_', queryOptoParameter , '_vs_', queryRunFeature, '_', 'colorComparison_', 'traces'];
plotFullPath = fullfile(savePlotPath, plotName);

if savePlots == 1
    saveas(gcf, plotFullPath, 'pdf')
    saveas(gcf, plotFullPath, 'fig')
end

%%

queryOptoParameter =  'amplitude';
queryRunFeature = 'areaPulse';

meanStrongBlueAreaPulse = nanmean(meansStrongBlueAreaPulse,1);
meanStrongRedAreaPulse = nanmean(meansStrongRedAreaPulse,1);
stdStrongBlueAreaPulse = nanstd(meansStrongBlueAreaPulse,1);
stdStrongRedAreaPulse = nanstd(meansStrongRedAreaPulse,1);

figure
scatter([0], meansStrongBlueAreaPulse,'k','filled');
hold on
scatter([1], meansStrongRedAreaPulse,'k','filled');
hold on
bar([0,1], [meanStrongBlueAreaPulse, meanStrongRedAreaPulse] ,'FaceColor',[0 0.4470 0.7410],'BarWidth', 0.5);
errorbar([0,1],[meanStrongBlueAreaPulse, meanStrongRedAreaPulse],[stdStrongBlueAreaPulse, stdStrongRedAreaPulse],'.','Color',[0.3 0.3 0.3],'LineWidth', 1.2, 'CapSize', 14)

hScatter = findobj(gca, 'Type', 'scatter');
uistack(hScatter, 'top');

xticks([0,1])
xticklabels({'Blue LED','Red LED'})
xlabel(['Activation channel']);
ylabel(queryRunFeature);
ylim([-2,14])
box on;
title([queryOptoParameter,' vs ', queryRunFeature]);
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');

plotName = [experimentName '_', queryOptoParameter , '_vs_', queryRunFeature, '_', 'colorComparison'];
plotFullPath = fullfile(savePlotPath, plotName);

if savePlots == 1
    saveas(gcf, plotFullPath, 'pdf')
    saveas(gcf, plotFullPath, 'fig')
end
