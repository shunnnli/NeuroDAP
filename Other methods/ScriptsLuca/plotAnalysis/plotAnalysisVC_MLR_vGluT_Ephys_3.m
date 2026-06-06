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
experimentName = 'vGluT_MLR_Ephys_3';

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

%% amplitudeBlue vs areaPulse - drugs

queryOptoParameter =  'amplitudeBlue';
queryRunFeature = 'areaPulse';

queryOptoParameterValuesOverall = [NaN, 5];
nNaNs = numel(find(isnan(queryOptoParameterValuesOverall)));

selectedCellsIDs = [5];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.holdingVoltage) == -70 & cell2mat(runGroupsParametersTable.pulseWidthBlue) == 100 & strcmp(runGroupsParametersTable.functionNameBlue,'dmdOutputRect(1, 1, 1, 1, 1, 1)'));
selectedRunGroupsAnalysisTable = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroups),:);

desiredEpochsDrugsRows = find(cell2mat(selectedRunGroupsAnalysisTable.epoch) >= 36);
selectedRunGroupsAnalysisTableDrugs = selectedRunGroupsAnalysisTable(desiredEpochsDrugsRows,:)

desiredEpochsRows = find(cell2mat(selectedRunGroupsAnalysisTable.epoch) > 3 & cell2mat(selectedRunGroupsAnalysisTable.epoch) <= 13);
selectedRunGroupsAnalysisTable = selectedRunGroupsAnalysisTable(desiredEpochsRows,:)

selectedCellsIDs = [2,3,4];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

controlRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.holdingVoltage) == -70 & cell2mat(runGroupsParametersTable.pulseWidthBlue) == 100 & strcmp(runGroupsParametersTable.functionNameBlue,'dmdOutputRect(1, 1, 1, 1, 1, 1)') & cell2mat(runGroupsParametersTable.nPulsesBlue) == 1);
controlRunGroupsAnalysisTable = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,controlRunGroups),:);

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

queryRunFeatureData = {};

for iRunGroup = 1:numel(selectedRunGroups)
    
    rowsRunGroup = find(selectedRunGroupsAnalysisTableDrugs.runGroup(:) == selectedRunGroups(iRunGroup));
    
    if ~isempty(rowsRunGroup)
    
        queryOptoParameterValues(groupsCounter) = selectedRunGroupsAnalysisTableDrugs.optoParameters{rowsRunGroup(1)}.(queryOptoParameter);
        queryRunFeatureData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTableDrugs.(queryRunFeature)(rowsRunGroup));
        cellNameData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTableDrugs.cellName(rowsRunGroup));

    end
    
    groupsCounter = groupsCounter + 1;
    
end

allDataDrugs = vertcat(queryRunFeatureData{:});
meanAllDataDrugs = nanmean(allDataDrugs,1);
stdAllDataDrugs = nanstd(allDataDrugs,1);

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

figure
scatter([0], allData,'k','filled');
hold on
scatter([1], allDataDrugs,'k','filled');
hold on
bar([0,1], [meanAllData, meanAllDataDrugs] ,'FaceColor',[0 0.4470 0.7410],'BarWidth', 0.5);
errorbar([0,1],[meanAllData, meanAllDataDrugs],[stdAllData, stdAllDataDrugs],'.','Color',[0.3 0.3 0.3],'LineWidth', 1.2, 'CapSize', 14)

hScatter = findobj(gca, 'Type', 'scatter');
uistack(hScatter, 'top');

xticks([0,1])
xticklabels({'No drug','NBQX + CPP'})
xlabel(['Condition']);
ylabel(queryRunFeature);
box on;
title([queryOptoParameter,' vs ', queryRunFeature]);
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');

plotName = [experimentName '_', queryRunFeature, '_', 'drugComparison'];
plotFullPath = fullfile(savePlotPath, plotName);

if savePlots == 1
    saveas(gcf, plotFullPath, 'pdf')
    saveas(gcf, plotFullPath, 'fig')
end

%% amplitudeBlue vs heightPulsePeak - drugs

queryOptoParameter =  'amplitudeBlue';
queryRunFeature = 'heightPulsePeak';

queryOptoParameterValuesOverall = [NaN, 5];
nNaNs = numel(find(isnan(queryOptoParameterValuesOverall)));

selectedCellsIDs = [5];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.holdingVoltage) == -70 & cell2mat(runGroupsParametersTable.pulseWidthBlue) == 100 & strcmp(runGroupsParametersTable.functionNameBlue,'dmdOutputRect(1, 1, 1, 1, 1, 1)'));
selectedRunGroupsAnalysisTable = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroups),:);

desiredEpochsDrugsRows = find(cell2mat(selectedRunGroupsAnalysisTable.epoch) >= 36);
selectedRunGroupsAnalysisTableDrugs = selectedRunGroupsAnalysisTable(desiredEpochsDrugsRows,:);

desiredEpochsRows = find(cell2mat(selectedRunGroupsAnalysisTable.epoch) > 3 & cell2mat(selectedRunGroupsAnalysisTable.epoch) <= 13);
selectedRunGroupsAnalysisTable = selectedRunGroupsAnalysisTable(desiredEpochsRows,:);

selectedCellsIDs = [2,3,4];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

controlRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.holdingVoltage) == -70 & cell2mat(runGroupsParametersTable.pulseWidthBlue) == 100 & strcmp(runGroupsParametersTable.functionNameBlue,'dmdOutputRect(1, 1, 1, 1, 1, 1)') & cell2mat(runGroupsParametersTable.nPulsesBlue) == 1);
controlRunGroupsAnalysisTable = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,controlRunGroups),:);

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
                
                trace = [nan(1,4000) trace nan(1,4000)];
                
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

queryRunFeatureData = {};
allTracesDrugs = [];

for iRunGroup = 1:numel(selectedRunGroups)
    
    rowsRunGroup = find(selectedRunGroupsAnalysisTableDrugs.runGroup(:) == selectedRunGroups(iRunGroup));
    
    if ~isempty(rowsRunGroup)
    
        queryOptoParameterValues(groupsCounter) = selectedRunGroupsAnalysisTableDrugs.optoParameters{rowsRunGroup(1)}.(queryOptoParameter);
        queryRunFeatureData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTableDrugs.(queryRunFeature)(rowsRunGroup));
        cellNameData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTableDrugs.cellName(rowsRunGroup));
        tracesGroup = selectedRunGroupsAnalysisTableDrugs.trace(rowsRunGroup)
        allTracesGroup = [];
        
        for iTrace = 1:numel(tracesGroup)
            
            trace = tracesGroup{iTrace};
            
            if numel(trace) < 10000
                
                trace = [zeros(1,5000) trace zeros(1,5000)];
                
            end
            
            allTracesGroup = [allTracesGroup; trace];
            
        end
        
        allTracesDrugs = [allTracesDrugs; allTracesGroup];
    
    end
    
    groupsCounter = groupsCounter + 1;
    
end

allDataDrugs = vertcat(queryRunFeatureData{:});
meanAllDataDrugs = nanmean(allDataDrugs,1);
stdAllDataDrugs = nanstd(allDataDrugs,1);

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

figure
scatter([0], allData,'k','filled');
hold on
scatter([1], allDataDrugs,'k','filled');
hold on
bar([0,1], [meanAllData, meanAllDataDrugs] ,'FaceColor',[0 0.4470 0.7410],'BarWidth', 0.5);
errorbar([0,1],[meanAllData, meanAllDataDrugs],[stdAllData, stdAllDataDrugs],'.','Color',[0.3 0.3 0.3],'LineWidth', 1.2, 'CapSize', 14)

hScatter = findobj(gca, 'Type', 'scatter');
uistack(hScatter, 'top');

xticks([0,1])
xticklabels({'No drug','NBQX + CPP'})
xlabel(['Condition']);
ylabel(queryRunFeature);
ylim([-2000,200])
box on;
title([queryOptoParameter,' vs ', queryRunFeature]);
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');

plotName = [experimentName '_', queryRunFeature, '_', 'drugComparison'];
plotFullPath = fullfile(savePlotPath, plotName);

if savePlots == 1
    saveas(gcf, plotFullPath, 'pdf')
    saveas(gcf, plotFullPath, 'fig')
end

figure;
plot(nanmean(allTraces,1),'LineWidth',2)
hold on

errorUp = nanmean(allTraces,1) + nanstd(allTraces,1);
errorDown = nanmean(allTraces,1) - nanstd(allTraces,1);
xLinspace = linspace(1,size(allTraces,2), size(allTraces,2));
plotpatch = patch([xLinspace fliplr(xLinspace)], [errorDown fliplr(errorUp)], [.63, .79, .95], 'FaceAlpha',0.5, 'EdgeColor','none');
plotpatch.Annotation.LegendInformation.IconDisplayStyle = 'off';

plot(nanmean(allTracesDrugs,1),'Color',[0.4 0.4 0.4],'LineWidth',2)

errorUp = nanmean(allTracesDrugs,1) + nanstd(allTracesDrugs,1);
errorDown = nanmean(allTracesDrugs,1) - nanstd(allTracesDrugs,1);
xLinspace = linspace(1,size(allTracesDrugs,2), size(allTracesDrugs,2));
plotpatch = patch([xLinspace fliplr(xLinspace)], [errorDown fliplr(errorUp)], [.8, .8, .8], 'FaceAlpha',0.5, 'EdgeColor','none');
plotpatch.Annotation.LegendInformation.IconDisplayStyle = 'off';


xlim([4000,7000])
ylim([-2000,200])
xticks([4000,5000,6000,7000])
xticklabels({'0.1','0.2','0.3','0.4'})
xlabel('Time [s]')
ylabel('Current [pA]')
box on;
title([queryOptoParameter,' vs ', queryRunFeature]);
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');
legend('No drugs','NBQX + CPP','Location', 'southeast', 'FontSize',12)

plotName = [experimentName '_', queryRunFeature, '_', 'drugComparison_' 'traces'];
plotFullPath = fullfile(savePlotPath, plotName);

if savePlots == 1
    saveas(gcf, plotFullPath, 'pdf')
    saveas(gcf, plotFullPath, 'fig')
end

%% amplitudeBlue vs areaPulse - cells

queryOptoParameter =  'amplitudeBlue';
queryRunFeature = 'areaPulse';

queryOptoParameterValuesOverall = [NaN, 5];
nNaNs = numel(find(isnan(queryOptoParameterValuesOverall)));

selectedCellsIDs = [5];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.holdingVoltage) == -70 & cell2mat(runGroupsParametersTable.pulseWidthBlue) == 100 & strcmp(runGroupsParametersTable.functionNameBlue,'dmdOutputRect(1, 1, 1, 1, 1, 1)'));
selectedRunGroupsAnalysisTable = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroups),:);

desiredEpochsRows = find(cell2mat(selectedRunGroupsAnalysisTable.epoch) > 3 & cell2mat(selectedRunGroupsAnalysisTable.epoch) <= 13);
selectedRunGroupsAnalysisTable = selectedRunGroupsAnalysisTable(desiredEpochsRows,:);

selectedCellsIDs = [2,3,4];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTableControl = mouseAnalysisTable(selectedCellsRows,:);

controlRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.holdingVoltage) == -70 & cell2mat(runGroupsParametersTable.pulseWidthBlue) == 100 & strcmp(runGroupsParametersTable.functionNameBlue,'dmdOutputRect(1, 1, 1, 1, 1, 1)') & cell2mat(runGroupsParametersTable.nPulsesBlue) == 1);
controlRunGroupsAnalysisTable = mouseSelectedCellsAnalysisTableControl(ismember(mouseSelectedCellsAnalysisTable.runGroup,controlRunGroups),:);

selectedCellsIDs = [1];
selectedCellsRestRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTableRest = mouseAnalysisTable(selectedCellsRestRows,:);

selectedRestRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.holdingVoltage) == -70 & cell2mat(runGroupsParametersTable.pulseWidthBlue) == 100 & strcmp(runGroupsParametersTable.functionNameBlue,'dmdOutputRect(1, 1, 1, 1, 1, 1)') & cell2mat(runGroupsParametersTable.nPulsesBlue) == 1);
selectedRunGroupsAnalysisTableRest = mouseSelectedCellsAnalysisTableRest(ismember(mouseSelectedCellsAnalysisTableRest.runGroup,selectedRestRunGroups),:);

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

queryRunFeatureData = {};

for iRunGroup = 1:numel(selectedRestRunGroups)
    
    rowsRunGroup = find(selectedRunGroupsAnalysisTableRest.runGroup(:) == selectedRestRunGroups(iRunGroup));
    
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

queryRunFeatureData = {};

for iRunGroup = 1:numel(controlRunGroups)
    
    rowsRunGroup = find(controlRunGroupsAnalysisTable.runGroup(:) == controlRunGroups(iRunGroup));
    
    if ~isempty(rowsRunGroup)
    
        queryOptoParameterValues(groupsCounter) = controlRunGroupsAnalysisTable.optoParameters{rowsRunGroup(1)}.(queryOptoParameter);
        queryRunFeatureData{groupsCounter} = cell2mat(controlRunGroupsAnalysisTable.(queryRunFeature)(rowsRunGroup));
        cellNameData{groupsCounter} = cell2mat(controlRunGroupsAnalysisTable.cellName(rowsRunGroup));
    
    end
    
    groupsCounter = groupsCounter + 1;
    
end

allDataControl = vertcat(queryRunFeatureData{:});
meanAllDataControl = nanmean(allDataControl,1);
stdAllDataControl = nanstd(allDataControl,1);

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

figure
scatter([0], allDataControl,'k','filled');
hold on
scatter([1], allDataRest,'k','filled');
hold on
scatter([2], allData,'k','filled');
hold on
bar([0,1,2], [meanAllDataControl, meanAllDataRest, meanAllData] ,'FaceColor',[0 0.4470 0.7410],'BarWidth', 0.5);
errorbar([0,1,2],[meanAllDataControl, meanAllDataRest, meanAllData],[stdAllDataControl, stdAllDataRest, stdAllData],'.','Color',[0.3 0.3 0.3],'LineWidth', 1.2, 'CapSize', 14)

hScatter = findobj(gca, 'Type', 'scatter');
uistack(hScatter, 'top');

xticks([0,1,2])
xticklabels({'Non-connected cells','Cell 1','Cell 5'})
xlabel(['Cells']);
ylabel(queryRunFeature);
box on;
title([queryOptoParameter,' vs ', queryRunFeature]);
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');

plotName = [experimentName '_', queryRunFeature '_' , 'cellComparison'];
plotFullPath = fullfile(savePlotPath, plotName);

if savePlots == 1
    saveas(gcf, plotFullPath, 'pdf')
    saveas(gcf, plotFullPath, 'fig')
end

%% amplitudeBlue vs heightPulsePeak - cells

queryOptoParameter =  'amplitudeBlue';
queryRunFeature = 'heightPulsePeak';

queryOptoParameterValuesOverall = [NaN, 5];
nNaNs = numel(find(isnan(queryOptoParameterValuesOverall)));

selectedCellsIDs = [5];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.holdingVoltage) == -70 & cell2mat(runGroupsParametersTable.pulseWidthBlue) == 100 & strcmp(runGroupsParametersTable.functionNameBlue,'dmdOutputRect(1, 1, 1, 1, 1, 1)'));
selectedRunGroupsAnalysisTable = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroups),:);

desiredEpochsRows = find(cell2mat(selectedRunGroupsAnalysisTable.epoch) > 3 & cell2mat(selectedRunGroupsAnalysisTable.epoch) <= 13);
selectedRunGroupsAnalysisTable = selectedRunGroupsAnalysisTable(desiredEpochsRows,:);

selectedCellsIDs = [2,3,4];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTableControl = mouseAnalysisTable(selectedCellsRows,:);

controlRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.holdingVoltage) == -70 & cell2mat(runGroupsParametersTable.pulseWidthBlue) == 100 & strcmp(runGroupsParametersTable.functionNameBlue,'dmdOutputRect(1, 1, 1, 1, 1, 1)') & cell2mat(runGroupsParametersTable.nPulsesBlue) == 1);
controlRunGroupsAnalysisTable = mouseSelectedCellsAnalysisTableControl(ismember(mouseSelectedCellsAnalysisTable.runGroup,controlRunGroups),:);

selectedCellsIDs = [1];
selectedCellsRestRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTableRest = mouseAnalysisTable(selectedCellsRestRows,:);

selectedRestRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.holdingVoltage) == -70 & cell2mat(runGroupsParametersTable.pulseWidthBlue) == 100 & strcmp(runGroupsParametersTable.functionNameBlue,'dmdOutputRect(1, 1, 1, 1, 1, 1)') & cell2mat(runGroupsParametersTable.nPulsesBlue) == 1);
selectedRunGroupsAnalysisTableRest = mouseSelectedCellsAnalysisTableRest(ismember(mouseSelectedCellsAnalysisTableRest.runGroup,selectedRestRunGroups),:);

queryOptoParameterValues = nan(numel(selectedRunGroups),1);
queryRunFeatureData = {};
cellNameData = {};
allTracesCell5 = [];

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
                
                trace = [zeros(1,4000) trace zeros(1,4000)];
                
            end
            
            allTracesGroup = [allTracesGroup; trace];
            
        end
        
        allTracesCell5 = [allTracesCell5; allTracesGroup];
    
    end
    
    groupsCounter = groupsCounter + 1;

end

allData = vertcat(queryRunFeatureData{:});
meanAllData = nanmean(allData,1);
stdAllData = nanstd(allData,1);

queryRunFeatureData = {};
allTracesRest = [];

for iRunGroup = 1:numel(selectedRestRunGroups)
    
    rowsRunGroup = find(selectedRunGroupsAnalysisTableRest.runGroup(:) == selectedRestRunGroups(iRunGroup));
    
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

queryRunFeatureData = {};
allTracesControl = [];

for iRunGroup = 1:numel(controlRunGroups)
    
    rowsRunGroup = find(controlRunGroupsAnalysisTable.runGroup(:) == controlRunGroups(iRunGroup));
    
    if ~isempty(rowsRunGroup)
    
        queryOptoParameterValues(groupsCounter) = controlRunGroupsAnalysisTable.optoParameters{rowsRunGroup(1)}.(queryOptoParameter);
        queryRunFeatureData{groupsCounter} = cell2mat(controlRunGroupsAnalysisTable.(queryRunFeature)(rowsRunGroup));
        cellNameData{groupsCounter} = cell2mat(controlRunGroupsAnalysisTable.cellName(rowsRunGroup));
        tracesGroup = controlRunGroupsAnalysisTable.trace(rowsRunGroup)
        allTracesGroup = [];
        
        for iTrace = 1:numel(tracesGroup)
            
            trace = tracesGroup{iTrace};
            
            if numel(trace) < 10000
                
                trace = [zeros(1,5000) trace zeros(1,5000)];
                
            end
            
            allTracesGroup = [allTracesGroup; trace];
            
        end
        
        allTracesControl = [allTracesControl; allTracesGroup];
    
    end
    
    groupsCounter = groupsCounter + 1;
    
end

allDataControl = vertcat(queryRunFeatureData{:});
meanAllDataControl = nanmean(allDataControl,1);
stdAllDataControl = nanstd(allDataControl,1);

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

figure
scatter([0], allDataControl,'k','filled');
hold on
scatter([1], allDataRest,'k','filled');
hold on
scatter([2], allData,'k','filled');
hold on
bar([0,1,2], [meanAllDataControl, meanAllDataRest, meanAllData] ,'FaceColor',[0 0.4470 0.7410],'BarWidth', 0.5);
errorbar([0,1,2],[meanAllDataControl, meanAllDataRest, meanAllData],[stdAllDataControl, stdAllDataRest, stdAllData],'.','Color',[0.3 0.3 0.3],'LineWidth', 1.2, 'CapSize', 14)

hScatter = findobj(gca, 'Type', 'scatter');
uistack(hScatter, 'top');

xticks([0,1,2])
xticklabels({'Non-connected cells','Cell 1','Cell 5'})
xlabel(['Cells']);
ylabel(queryRunFeature);
box on;
title([queryOptoParameter,' vs ', queryRunFeature]);
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');

plotName = [experimentName '_', queryRunFeature '_' , 'cellComparison'];
plotFullPath = fullfile(savePlotPath, plotName);

if savePlots == 1
    saveas(gcf, plotFullPath, 'pdf')
    saveas(gcf, plotFullPath, 'fig')
end

figure;

plot(nanmean(allTracesControl,1),'Color',[0.4 0.4 0.4],'LineWidth',2)
hold on
errorUp = nanmean(allTracesControl,1) + nanstd(allTracesControl,1);
errorDown = nanmean(allTracesControl,1) - nanstd(allTracesControl,1);
xLinspace = linspace(1,size(allTracesControl,2), size(allTracesControl,2));
% plotpatch = patch([xLinspace fliplr(xLinspace)], [errorDown fliplr(errorUp)], [.8, .8, .8], 'FaceAlpha',0.5, 'EdgeColor','none');
% plotpatch.Annotation.LegendInformation.IconDisplayStyle = 'off';

plot(nanmean(allTracesRest,1),'LineWidth',2)
hold on

errorUp = nanmean(allTracesRest,1) + nanstd(allTracesRest,1);
errorDown = nanmean(allTracesRest,1) - nanstd(allTracesRest,1);
xLinspace = linspace(1,size(allTracesRest,2), size(allTracesRest,2));
% plotpatch = patch([xLinspace fliplr(xLinspace)], [errorDown fliplr(errorUp)], [.6, .1, .1], 'FaceAlpha',0.5, 'EdgeColor','none');
% plotpatch.Annotation.LegendInformation.IconDisplayStyle = 'off';

plot(nanmean(allTracesCell5,1),'Color',[0 0.4470 0.7410],'LineWidth',2)
hold on

errorUp = nanmean(allTracesCell5,1) + nanstd(allTracesCell5,1);
errorDown = nanmean(allTracesCell5,1) - nanstd(allTracesCell5,1);
xLinspace = linspace(1,size(allTracesCell5,2), size(allTracesCell5,2));
% plotpatch = patch([xLinspace fliplr(xLinspace)], [errorDown fliplr(errorUp)], [.63, .79, .95], 'FaceAlpha',0.5, 'EdgeColor','none');
% plotpatch.Annotation.LegendInformation.IconDisplayStyle = 'off';

xlim([4000,7000])
ylim([-2000,200])
xticks([4000,5000,6000,7000])
xticklabels({'0.1','0.2','0.3','0.4'})
xlabel('Time [s]')
ylabel('Current [pA]')
box on;
title([queryOptoParameter,' vs ', queryRunFeature]);
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');
legend('Non-connected cells','Cell 1','Cell 5','Location', 'southeast', 'FontSize',12)

plotName = [experimentName '_' queryRunFeature, '_', 'cellComparison_','traces'];
plotFullPath = fullfile(savePlotPath, plotName);

if savePlots == 1
    saveas(gcf, plotFullPath, 'pdf')
    saveas(gcf, plotFullPath, 'fig')
end

%% amplitudeBlue vs heightPulsePeak - colors

queryOptoParameter =  'amplitudeBlue';
queryRunFeature = 'heightPulsePeak';

queryOptoParameterValuesOverall = [NaN, 5];
nNaNs = numel(find(isnan(queryOptoParameterValuesOverall)));

selectedCellsIDs = [5];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.holdingVoltage) == -70 & cell2mat(runGroupsParametersTable.pulseWidthBlue) == 100 & strcmp(runGroupsParametersTable.functionNameBlue,'dmdOutputRect(1, 1, 1, 1, 1, 1)'));
selectedRunGroupsAnalysisTable = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroups),:);

desiredEpochsDrugsRows = find(cell2mat(selectedRunGroupsAnalysisTable.epoch) >= 36);
selectedRunGroupsAnalysisTableDrugs = selectedRunGroupsAnalysisTable(desiredEpochsDrugsRows,:);

desiredEpochsRows = find(cell2mat(selectedRunGroupsAnalysisTable.epoch) > 3 & cell2mat(selectedRunGroupsAnalysisTable.epoch) <= 13);
selectedRunGroupsAnalysisTable = selectedRunGroupsAnalysisTable(desiredEpochsRows,:);

selectedCellsIDs = [5];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

controlRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao2') & cell2mat(runGroupsParametersTable.holdingVoltage) == -70 & cell2mat(runGroupsParametersTable.pulseWidthRed) == 100 & strcmp(runGroupsParametersTable.functionNameRed,'dmdOutputRect(1, 1, 1, 1, 2, 1)') & cell2mat(runGroupsParametersTable.nPulsesRed) == 1);
controlRunGroupsAnalysisTable = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,controlRunGroups),:);

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
                
                trace = [nan(1,4000) trace nan(1,4000)];
                
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

queryRunFeatureData = {};
allTracesRed = [];

for iRunGroup = 1:numel(selectedRunGroups)
    
    rowsRunGroup = find(selectedRunGroupsAnalysisTableDrugs.runGroup(:) == selectedRunGroups(iRunGroup));
    
    if ~isempty(rowsRunGroup)
    
        queryOptoParameterValues(groupsCounter) = selectedRunGroupsAnalysisTableDrugs.optoParameters{rowsRunGroup(1)}.(queryOptoParameter);
        queryRunFeatureData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTableDrugs.(queryRunFeature)(rowsRunGroup));
        cellNameData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTableDrugs.cellName(rowsRunGroup));
        tracesGroup = selectedRunGroupsAnalysisTableDrugs.trace(rowsRunGroup)
        allTracesGroup = [];
        
        for iTrace = 1:numel(tracesGroup)
            
            trace = tracesGroup{iTrace};
            
            if numel(trace) < 10000
                
                trace = [zeros(1,5000) trace zeros(1,5000)];
                
            end
            
            allTracesGroup = [allTracesGroup; trace];
            
        end
        
        allTracesRed = [allTracesRed; allTracesGroup];
    
    end
    
    groupsCounter = groupsCounter + 1;
    
end

allDataRed = vertcat(queryRunFeatureData{:});
meanAllDataRed = nanmean(allDataRed,1);
stdAllDataRed = nanstd(allDataRed,1);

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

figure
scatter([0], allData,'k','filled');
hold on
scatter([1], allDataRed,'k','filled');
hold on
bar([0,1], [meanAllData, meanAllDataRed] ,'FaceColor',[0 0.4470 0.7410],'BarWidth', 0.5);
errorbar([0,1],[meanAllData, meanAllDataRed],[stdAllData, stdAllDataRed],'.','Color',[0.3 0.3 0.3],'LineWidth', 1.2, 'CapSize', 14)

hScatter = findobj(gca, 'Type', 'scatter');
uistack(hScatter, 'top');

xticks([0,1])
xticklabels({'Blue LED','Red LED'})
xlabel(['Condition']);
ylabel(queryRunFeature);
ylim([-2000,200])
box on;
title([queryOptoParameter,' vs ', queryRunFeature]);
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');

plotName = [experimentName '_', queryRunFeature, '_', 'colorComparison'];
plotFullPath = fullfile(savePlotPath, plotName);

if savePlots == 1
    saveas(gcf, plotFullPath, 'pdf')
    saveas(gcf, plotFullPath, 'fig')
end

figure;
plot(nanmean(allTraces,1),'LineWidth',2)
hold on

errorUp = nanmean(allTraces,1) + nanstd(allTraces,1);
errorDown = nanmean(allTraces,1) - nanstd(allTraces,1);
xLinspace = linspace(1,size(allTraces,2), size(allTraces,2));
plotpatch = patch([xLinspace fliplr(xLinspace)], [errorDown fliplr(errorUp)], [.63, .79, .95], 'FaceAlpha',0.5, 'EdgeColor','none');
plotpatch.Annotation.LegendInformation.IconDisplayStyle = 'off';

plot(nanmean(allTracesRed,1),'LineWidth',2)

errorUp = nanmean(allTracesRed,1) + nanstd(allTracesRed,1);
errorDown = nanmean(allTracesRed,1) - nanstd(allTracesRed,1);
xLinspace = linspace(1,size(allTracesRed,2), size(allTracesRed,2));
plotpatch = patch([xLinspace fliplr(xLinspace)], [errorDown fliplr(errorUp)], [.6, .1, .1], 'FaceAlpha',0.5, 'EdgeColor','none');
plotpatch.Annotation.LegendInformation.IconDisplayStyle = 'off';


xlim([4000,7000])
ylim([-2000,200])
xticks([4000,5000,6000,7000])
xticklabels({'0.1','0.2','0.3','0.4'})
xlabel('Time [s]')
ylabel('Current [pA]')
box on;
title([queryOptoParameter,' vs ', queryRunFeature]);
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');
legend('Blue LED','Red LED','Location', 'southeast', 'FontSize',12)

plotName = [experimentName '_', queryRunFeature, '_', 'colorComparison_' 'traces'];
plotFullPath = fullfile(savePlotPath, plotName);

if savePlots == 1
    saveas(gcf, plotFullPath, 'pdf')
    saveas(gcf, plotFullPath, 'fig')
end


%% amplitudeBlue vs areaPulse - colors

queryOptoParameter =  'amplitudeBlue';
queryRunFeature = 'areaPulse';

queryOptoParameterValuesOverall = [NaN, 5];
nNaNs = numel(find(isnan(queryOptoParameterValuesOverall)));

selectedCellsIDs = [5];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.holdingVoltage) == -70 & cell2mat(runGroupsParametersTable.pulseWidthBlue) == 100 & strcmp(runGroupsParametersTable.functionNameBlue,'dmdOutputRect(1, 1, 1, 1, 1, 1)'));
selectedRunGroupsAnalysisTable = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroups),:);

desiredEpochsDrugsRows = find(cell2mat(selectedRunGroupsAnalysisTable.epoch) >= 36);
selectedRunGroupsAnalysisTableDrugs = selectedRunGroupsAnalysisTable(desiredEpochsDrugsRows,:);

desiredEpochsRows = find(cell2mat(selectedRunGroupsAnalysisTable.epoch) > 3 & cell2mat(selectedRunGroupsAnalysisTable.epoch) <= 13);
selectedRunGroupsAnalysisTable = selectedRunGroupsAnalysisTable(desiredEpochsRows,:);

selectedCellsIDs = [5];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

controlRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao2') & cell2mat(runGroupsParametersTable.holdingVoltage) == -70 & cell2mat(runGroupsParametersTable.pulseWidthRed) == 100 & strcmp(runGroupsParametersTable.functionNameRed,'dmdOutputRect(1, 1, 1, 1, 2, 1)') & cell2mat(runGroupsParametersTable.nPulsesRed) == 1);
controlRunGroupsAnalysisTable = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,controlRunGroups),:);

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
                
                trace = [nan(1,4000) trace nan(1,4000)];
                
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

queryRunFeatureData = {};
allTracesRed = [];

for iRunGroup = 1:numel(selectedRunGroups)
    
    rowsRunGroup = find(selectedRunGroupsAnalysisTableDrugs.runGroup(:) == selectedRunGroups(iRunGroup));
    
    if ~isempty(rowsRunGroup)
    
        queryOptoParameterValues(groupsCounter) = selectedRunGroupsAnalysisTableDrugs.optoParameters{rowsRunGroup(1)}.(queryOptoParameter);
        queryRunFeatureData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTableDrugs.(queryRunFeature)(rowsRunGroup));
        cellNameData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTableDrugs.cellName(rowsRunGroup));
        tracesGroup = selectedRunGroupsAnalysisTableDrugs.trace(rowsRunGroup)
        allTracesGroup = [];
        
        for iTrace = 1:numel(tracesGroup)
            
            trace = tracesGroup{iTrace};
            
            if numel(trace) < 10000
                
                trace = [zeros(1,5000) trace zeros(1,5000)];
                
            end
            
            allTracesGroup = [allTracesGroup; trace];
            
        end
        
        allTracesRed = [allTracesRed; allTracesGroup];
    
    end
    
    groupsCounter = groupsCounter + 1;
    
end

allDataRed = vertcat(queryRunFeatureData{:});
meanAllDataRed = nanmean(allDataRed,1);
stdAllDataRed = nanstd(allDataRed,1);

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

figure
scatter([0], allData,'k','filled');
hold on
scatter([1], allDataRed,'k','filled');
hold on
bar([0,1], [meanAllData, meanAllDataRed] ,'FaceColor',[0 0.4470 0.7410],'BarWidth', 0.5);
errorbar([0,1],[meanAllData, meanAllDataRed],[stdAllData, stdAllDataRed],'.','Color',[0.3 0.3 0.3],'LineWidth', 1.2, 'CapSize', 14)

hScatter = findobj(gca, 'Type', 'scatter');
uistack(hScatter, 'top');

xticks([0,1])
xticklabels({'Blue LED','Red LED'})
xlabel(['Condition']);
ylabel(queryRunFeature);
box on;
title([queryOptoParameter,' vs ', queryRunFeature]);
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');

plotName = [experimentName '_', queryRunFeature, '_', 'colorComparison'];
plotFullPath = fullfile(savePlotPath, plotName);

if savePlots == 1
    saveas(gcf, plotFullPath, 'pdf')
    saveas(gcf, plotFullPath, 'fig')
end
