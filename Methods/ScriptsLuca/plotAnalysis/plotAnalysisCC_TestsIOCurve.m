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
experimentName = 'Pitx_GtACR_Ephys_2';%'vGAT_GtACR_Ephys_1_240319';

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

%% amplitudeCurrentPulse vs areaVtPulse -- withLightBeforeCurrent -- line plot

queryOptoParameter =  'amplitudeCurrentPulse';
queryRunFeature = 'areaVtPulse';

selectedCellsIDs = [1,2,3,5];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.widthCurrentPulse) == 6000 & cell2mat(runGroupsParametersTable.delayLightPulse) == 19000 & cell2mat(runGroupsParametersTable.amplitudeLightPulse) == 5);
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
        %plotscatter = scatter(queryOptoParameterValues(iQueryOptoParameter), meanQueryRunFeatureDataCell, 40, 'MarkerFaceColor', [0 0.4 0.7], 'MarkerEdgeColor', [0 0.4 0.7]);
        %plotscatter.Annotation.LegendInformation.IconDisplayStyle = 'off';
        %hold on;
        
        meansQueryRunFeatureDataAllCells(iCell,positionOptoParameterValue) = meanQueryRunFeatureDataCell;
        
    end
    
end

meanQueryRunFeatureDataAllCells = nanmean(meansQueryRunFeatureDataAllCells,1);
stdQueryRunFeatureDataAllCells = nanstd(meansQueryRunFeatureDataAllCells,1);
xValues = linspace(1,numel(queryOptoParameterValues),numel(queryOptoParameterValues));
plot(queryOptoParameterValues, meanQueryRunFeatureDataAllCells,'Color',[0 0.4470 0.7410],'LineWidth',2);

errorUp = meanQueryRunFeatureDataAllCells + stdQueryRunFeatureDataAllCells;
errorDown = meanQueryRunFeatureDataAllCells - stdQueryRunFeatureDataAllCells;
plotpatch = patch([queryOptoParameterValues' fliplr(queryOptoParameterValues')], [errorDown fliplr(errorUp)], [.63, .79, .95], 'FaceAlpha',0.5, 'EdgeColor','none');
plotpatch.Annotation.LegendInformation.IconDisplayStyle = 'off';


% ploterrorbar = errorbar(queryOptoParameterValues,meanQueryRunFeatureDataAllCells,stdQueryRunFeatureDataAllCells,'.','Color',[0 0 0],'LineWidth', 1.2, 'CapSize', 20);
% ploterrorbar.Annotation.LegendInformation.IconDisplayStyle = 'off';
hold on;

%hScatter = findobj(gca, 'Type', 'scatter');
%uistack(hScatter, 'top');

% xticks(xValues)
% xticklabels([unique(queryOptoParameterValues)])
xlabel(queryOptoParameter);
ylabel(queryRunFeature);
box on;
title([queryOptoParameter,' vs ', queryRunFeature]);
%set(gca, 'XScale', 'log');
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');

%traceWithLight = selectedRunGroupsAnalysisTable.trace{81,1};
traceWithLight = selectedRunGroupsAnalysisTable.trace{39,1};%62

%%

queryOptoParameter =  'amplitudeCurrentPulse';
queryRunFeature = 'areaVtPulse';

selectedCellsIDs = [1,2,3,5];
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
        %plotscatter = scatter(queryOptoParameterValues(iQueryOptoParameter), meanQueryRunFeatureDataCell, 40, 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', [0.5 0.5 0.5]);
        %plotscatter.Annotation.LegendInformation.IconDisplayStyle = 'off';
        %hold on;
        
        meansQueryRunFeatureDataAllCells(iCell,positionOptoParameterValue) = meanQueryRunFeatureDataCell;
        
    end
    
end

meanQueryRunFeatureDataAllCells = nanmean(meansQueryRunFeatureDataAllCells,1);
stdQueryRunFeatureDataAllCells = nanstd(meansQueryRunFeatureDataAllCells,1);
xValues = linspace(1,numel(queryOptoParameterValues),numel(queryOptoParameterValues));

errorUp = meanQueryRunFeatureDataAllCells + stdQueryRunFeatureDataAllCells;
errorDown = meanQueryRunFeatureDataAllCells - stdQueryRunFeatureDataAllCells;
plotpatch = patch([queryOptoParameterValues' fliplr(queryOptoParameterValues')], [errorDown fliplr(errorUp)], [.8, .8, .8], 'FaceAlpha',0.5, 'EdgeColor','none');
plotpatch.Annotation.LegendInformation.IconDisplayStyle = 'off';


plot(queryOptoParameterValues, meanQueryRunFeatureDataAllCells,'Color',[0.7 0.7 0.7],'LineWidth',2);


% ploterrorbar = errorbar(queryOptoParameterValues,meanQueryRunFeatureDataAllCells,stdQueryRunFeatureDataAllCells,'.','Color',[0 0 0],'LineWidth', 1.2, 'CapSize', 20);
% ploterrorbar.Annotation.LegendInformation.IconDisplayStyle = 'off';

% hScatter = findobj(gca, 'Type', 'scatter');
% uistack(hScatter, 'top');

% xticks(xValues)
% xticklabels([unique(queryOptoParameterValues)])
%xlim([0.5, 9.5])
xlim([-200, 1200])
xlabel(queryOptoParameter);
ylabel(queryRunFeature);
box on;
title([queryOptoParameter,' vs ', queryRunFeature]);
legend('with Light','no Light','Location', 'northwest', 'FontSize',12)
%set(gca, 'XScale', 'log');
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');

plotName = [experimentName '_' typeCC '_' queryOptoParameter,'_vs_', queryRunFeature];
plotFullPath = fullfile(savePlotPath, plotName);
% saveas(gcf, plotFullPath, 'pdf');
% saveas(gcf, plotFullPath, 'fig');

%traceNoLight = selectedRunGroupsAnalysisTable.trace{81,1};
traceNoLight = selectedRunGroupsAnalysisTable.trace{163,1};%117

%%

figure;
plot(traceWithLight,'Color',[0 0.4470 0.7410],'LineWidth',2)
hold on;
plot(traceNoLight,'Color',[0.7 0.7 0.7],'LineWidth',2)

xticks = get(gca, 'XTick');
set(gca, 'XTickLabel', num2str(xticks'/10000));
xlim([18000, 32000]);
xlabel('Time [s]');
ylabel('Membrane potential [mV]');
box on;
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');
legend('with Light','no Light','Location', 'northwest', 'FontSize',12)

plotName = [experimentName '_' typeCC '_traces'];
plotFullPath = fullfile(savePlotPath, plotName);
saveas(gcf, plotFullPath, 'pdf');
saveas(gcf, plotFullPath, 'fig');