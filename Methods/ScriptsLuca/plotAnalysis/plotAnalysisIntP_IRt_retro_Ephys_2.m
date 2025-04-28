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
experimentName = 'IRt_retro_Ephys_2';

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
savePlots = 1;

%% amplitudeCurrentPulse

queryOptoParameter =  'amplitudeCurrentPulse';

selectedCellsIDs = [1,2,3];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0') & cell2mat(runGroupsParametersTable.widthCurrentPulse) == 6000);
selectedRunGroupsAnalysisTable = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroups),:);

queryOptoParameterValues = zeros(numel(selectedRunGroups),1);
queryOptoParameterArray = [];
nAPallData = {};
timesAPallData = {};
timeFirstAPallData = {};
timeLastAPallData = {};
peakVoltageAPallData = {};
thresholdAPallData = {};
cellNameData= {};

groupsCounter = 1;

for iRunGroup = 1:numel(selectedRunGroups)
    
    rowsRunGroup = find(selectedRunGroupsAnalysisTable.runGroup(:) == selectedRunGroups(iRunGroup));
    queryOptoParameterValues(groupsCounter) = selectedRunGroupsAnalysisTable.optoParameters{rowsRunGroup(1)}.(queryOptoParameter);
    
    queryOptoParameterArray = [queryOptoParameterArray; repmat(queryOptoParameterValues(groupsCounter),numel(rowsRunGroup),1)];
    pulseAProwsGroup = selectedRunGroupsAnalysisTable.pulseAP(rowsRunGroup);
    
    nAProwsGroup = {};
    timesAProwsGroup = {};
    timeFirstAProwsGroup = {};
    timeLastAProwsGroup = {};
    peakVoltageAProwsGroup = {};
    thresholdAProwsGroup = {};
      
    for iRun = 1:numel(rowsRunGroup)
    
      if ~isempty(pulseAProwsGroup{iRun})
          
          peakVoltageAP = pulseAProwsGroup{iRun}.AP_peak_V;
          indicesCorrectPeaks = find(peakVoltageAP>-15);
    
          nAP = numel(indicesCorrectPeaks);
          timesAP = pulseAProwsGroup{iRun}.AP_peak_time(indicesCorrectPeaks);
          timeFirstAP = timesAP(1);
          timeLastAP = timesAP(end);
          peakVoltageAP = pulseAProwsGroup{iRun}.AP_peak_V(indicesCorrectPeaks);
          thresholdAP = pulseAProwsGroup{iRun}.AP_thresh_V(indicesCorrectPeaks);
    
%           nAP = pulseAProwsGroup{iRun}.nAP;
%           timesAP = pulseAProwsGroup{iRun}.AP_peak_time;
%           timeFirstAP = pulseAProwsGroup{iRun}.AP_peak_time(1);
%           timeLastAP = pulseAProwsGroup{iRun}.AP_peak_time(end);
%           peakVoltageAP = pulseAProwsGroup{iRun}.AP_peak_V;
%           thresholdAP = pulseAProwsGroup{iRun}.AP_thresh_V;
    
      else
    
          nAP = 0;
          timesAP = nan;
          timeFirstAP = nan;
          timeLastAP = nan;
          peakVoltageAP = nan;
          thresholdAP = nan;
    
      end
    
      nAProwsGroup{iRun} = nAP;
      timesAProwsGroup{iRun} = timesAP;
      timeFirstAProwsGroup{iRun} = timeFirstAP;
      timeLastAProwsGroup{iRun} = timeLastAP;
      peakVoltageAProwsGroup{iRun} = peakVoltageAP;
      thresholdAProwsGroup{iRun} = thresholdAP;
    
    end
    
    nAPallData{groupsCounter} = nAProwsGroup;
    timesAPallData{groupsCounter} = timesAProwsGroup;
    timeFirstAPallData{groupsCounter} = timeFirstAProwsGroup;
    timeLastAPallData{groupsCounter} = timeLastAProwsGroup;
    peakVoltageAPallData{groupsCounter} = peakVoltageAProwsGroup;
    thresholdAPallData{groupsCounter} = thresholdAProwsGroup;
    cellNameData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTable.cellName(rowsRunGroup))';
    
    groupsCounter = groupsCounter + 1;
    
end

%% timeFirstAP

queryOptoParameter =  'amplitudeCurrentPulse';
queryRunFeature = 'timeFirstAP';
allData = timeFirstAPallData;

meansAllDataCells = nan(numel(selectedCellsIDs),numel(selectedRunGroups));
stdsAllDataCells = nan(numel(selectedCellsIDs),numel(selectedRunGroups));
figure;
    
for iRunGroup = 1:numel(selectedRunGroups)
    
    cellNames = unique(cellNameData{1,iRunGroup});
    
    for iCell = 1:numel(cellNames)
        
        groupCells = vertcat(cellNameData{1,iRunGroup});
        groupCellRows = find(groupCells == cellNames(iCell));
        
        allDataCells = vertcat(allData{1,iRunGroup}{:});
        allDataCell = allDataCells(groupCellRows);
        
        positionOptoParameterValue = find(sort(queryOptoParameterValues) == queryOptoParameterValues(iRunGroup));
        
        meansAllDataCells(find(selectedCellsIDs == cellNames(iCell)),positionOptoParameterValue) = nanmean(allDataCell);
        stdsAllDataCells(find(selectedCellsIDs == cellNames(iCell)),positionOptoParameterValue) = nanstd(allDataCell);
       
        hscatter = scatter(positionOptoParameterValue, nanmean(allDataCell),'k','filled');
        hold on;
    
    end

end

meanAllDataCells = nanmean(meansAllDataCells,1);
stdAllDataCells = nanstd(meansAllDataCells,1);
xValues = linspace(1,numel(queryOptoParameterValues),numel(queryOptoParameterValues));
bar(xValues,meanAllDataCells,'FaceColor',[0 0.4470 0.7410],'BarWidth', 0.5)

if numel(selectedCellsIDs) > 1
    errorbar(xValues,meanAllDataCells,stdAllDataCells,'.','Color',[0 0 0],'LineWidth', 1.2, 'CapSize', 14)
end

hScatter = findobj(gca, 'Type', 'scatter');
uistack(hScatter, 'top');

xticks(xValues)
xticklabels([unique(queryOptoParameterValues)])
xlabel(queryOptoParameter);
ylabel(queryRunFeature);
box on;
title([queryOptoParameter,' vs ', queryRunFeature]);
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');
ylim([0,600])

plotName = ['IntP' '_' queryRunFeature,'_vs_', queryOptoParameter '_' experimentName];
plotFullPath = fullfile(savePlotPath, plotName);

if savePlots == 1
    saveas(gcf, plotFullPath, 'pdf')
    saveas(gcf, plotFullPath, 'fig')
end

%% timeLastAP

queryOptoParameter =  'amplitudeCurrentPulse';
queryRunFeature = 'timeLastAP';
allData = timeLastAPallData;

meansAllDataCells = nan(numel(selectedCellsIDs),numel(selectedRunGroups));
stdsAllDataCells = nan(numel(selectedCellsIDs),numel(selectedRunGroups));
figure;
    
for iRunGroup = 1:numel(selectedRunGroups)
    
    cellNames = unique(cellNameData{1,iRunGroup});
    
    for iCell = 1:numel(cellNames)
        
        groupCells = vertcat(cellNameData{1,iRunGroup});
        groupCellRows = find(groupCells == cellNames(iCell));
        
        allDataCells = vertcat(allData{1,iRunGroup}{:});
        allDataCell = allDataCells(groupCellRows);
        
        positionOptoParameterValue = find(sort(queryOptoParameterValues) == queryOptoParameterValues(iRunGroup));
        
        meansAllDataCells(find(selectedCellsIDs == cellNames(iCell)),positionOptoParameterValue) = nanmean(allDataCell);
        stdsAllDataCells(find(selectedCellsIDs == cellNames(iCell)),positionOptoParameterValue) = nanstd(allDataCell);
       
        hscatter = scatter(positionOptoParameterValue, nanmean(allDataCell),'k','filled');
        hold on;
    
    end

end

meanAllDataCells = nanmean(meansAllDataCells,1);
stdAllDataCells = nanstd(meansAllDataCells,1);
xValues = linspace(1,numel(queryOptoParameterValues),numel(queryOptoParameterValues));
bar(xValues,meanAllDataCells,'FaceColor',[0 0.4470 0.7410],'BarWidth', 0.5)

if numel(selectedCellsIDs) > 1
    errorbar(xValues,meanAllDataCells,stdAllDataCells,'.','Color',[0 0 0],'LineWidth', 1.2, 'CapSize', 14)
end

hScatter = findobj(gca, 'Type', 'scatter');
uistack(hScatter, 'top');

xticks(xValues)
xticklabels([unique(queryOptoParameterValues)])
xlabel(queryOptoParameter);
ylabel(queryRunFeature);
box on;
title([queryOptoParameter,' vs ', queryRunFeature]);
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');
ylim([0,700])

plotName = ['IntP' '_' queryRunFeature,'_vs_', queryOptoParameter '_' experimentName];
plotFullPath = fullfile(savePlotPath, plotName);

if savePlots == 1
    saveas(gcf, plotFullPath, 'pdf')
    saveas(gcf, plotFullPath, 'fig')
end

%% nAP

queryOptoParameter =  'amplitudeCurrentPulse';
queryRunFeature = 'nAP';
allData = nAPallData;

meansAllDataCells = nan(numel(selectedCellsIDs),numel(selectedRunGroups));
stdsAllDataCells = nan(numel(selectedCellsIDs),numel(selectedRunGroups));
figure;
    
for iRunGroup = 1:numel(selectedRunGroups)
    
    cellNames = unique(cellNameData{1,iRunGroup});
    
    for iCell = 1:numel(cellNames)
        
        groupCells = vertcat(cellNameData{1,iRunGroup});
        groupCellRows = find(groupCells == cellNames(iCell));
        
        allDataCells = vertcat(allData{1,iRunGroup}{:});
        allDataCell = allDataCells(groupCellRows);
        
        positionOptoParameterValue = find(sort(queryOptoParameterValues) == queryOptoParameterValues(iRunGroup));
        
        meansAllDataCells(find(selectedCellsIDs == cellNames(iCell)),positionOptoParameterValue) = nanmean(allDataCell);
        stdsAllDataCells(find(selectedCellsIDs == cellNames(iCell)),positionOptoParameterValue) = nanstd(allDataCell);
       
        hscatter = scatter(positionOptoParameterValue, nanmean(allDataCell),'k','filled');
        hold on;
    
    end

end

meanAllDataCells = nanmean(meansAllDataCells,1);
stdAllDataCells = nanstd(meansAllDataCells,1);
xValues = linspace(1,numel(queryOptoParameterValues),numel(queryOptoParameterValues));
bar(xValues,meanAllDataCells,'FaceColor',[0 0.4470 0.7410],'BarWidth', 0.5)
if numel(selectedCellsIDs) > 1
    errorbar(xValues,meanAllDataCells,stdAllDataCells,'.','Color',[0 0 0],'LineWidth', 1.2, 'CapSize', 14)
end

hScatter = findobj(gca, 'Type', 'scatter');
uistack(hScatter, 'top');

xticks(xValues)
xticklabels([unique(queryOptoParameterValues)])
xlabel(queryOptoParameter);
ylabel(queryRunFeature);
box on;
title([queryOptoParameter,' vs ', queryRunFeature]);
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');
ylim([0,80])

plotName = ['IntP' '_' queryRunFeature,'_vs_', queryOptoParameter '_' experimentName];
plotFullPath = fullfile(savePlotPath, plotName);

if savePlots == 1
    saveas(gcf, plotFullPath, 'pdf')
    saveas(gcf, plotFullPath, 'fig')
end

%% nAP line

queryOptoParameter =  'amplitudeCurrentPulse';
queryRunFeature = 'nAP';
allData = nAPallData;

meansAllDataCells = nan(numel(selectedCellsIDs),numel(selectedRunGroups));
stdsAllDataCells = nan(numel(selectedCellsIDs),numel(selectedRunGroups));
figure;
    
for iRunGroup = 1:numel(selectedRunGroups)
    
    cellNames = unique(cellNameData{1,iRunGroup});
    
    for iCell = 1:numel(cellNames)
        
        groupCells = vertcat(cellNameData{1,iRunGroup});
        groupCellRows = find(groupCells == cellNames(iCell));
        
        allDataCells = vertcat(allData{1,iRunGroup}{:});
        allDataCell = allDataCells(groupCellRows);
        
        positionOptoParameterValue = find(sort(queryOptoParameterValues) == queryOptoParameterValues(iRunGroup));
        
        meansAllDataCells(find(selectedCellsIDs == cellNames(iCell)),positionOptoParameterValue) = nanmean(allDataCell);
        stdsAllDataCells(find(selectedCellsIDs == cellNames(iCell)),positionOptoParameterValue) = nanstd(allDataCell);
       
        hscatter = scatter(positionOptoParameterValue, nanmean(allDataCell),'k','filled');
        hold on;
    
    end

end

meanAllDataCells = nanmean(meansAllDataCells,1);
stdAllDataCells = nanstd(meansAllDataCells,1);
xValues = linspace(1,numel(queryOptoParameterValues),numel(queryOptoParameterValues));
plot(xValues,meanAllDataCells,'Color',[0 0.4470 0.7410],'LineWidth', 2)

if numel(selectedCellsIDs) > 1
    errorbar(xValues,meanAllDataCells,stdAllDataCells,'.','Color',[0 0 0],'LineWidth', 1.2, 'CapSize', 14)
end

hScatter = findobj(gca, 'Type', 'scatter');
uistack(hScatter, 'top');

xticks(xValues)
xticklabels([unique(queryOptoParameterValues)])
xlabel(queryOptoParameter);
ylabel(queryRunFeature);
box on;
title([queryOptoParameter,' vs ', queryRunFeature]);
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');
ylim([0,80])

plotName = ['IntP' '_' queryRunFeature,'_vs_', queryOptoParameter '_line_' experimentName];
plotFullPath = fullfile(savePlotPath, plotName);

if savePlots == 1
    saveas(gcf, plotFullPath, 'pdf')
    saveas(gcf, plotFullPath, 'fig')
end

%%

queryOptoParameter =  'amplitudeCurrentPulse';
queryRunFeature = 'timesAP';
allData = timesAPallData;

edges = [0    100   200   300   400   500   600];

meanAllData = [];
stdAllData = [];

figure;
cMap = flip(hsv(numel(selectedRunGroups)));

for iRunGroup = 1:numel(selectedRunGroups)
    
    cellNames = unique(cellNameData{1,iRunGroup});
    meansAllDataCells = [];
    
    for iCell = 1:numel(cellNames)
        
        groupCells = vertcat(cellNameData{1,iRunGroup});
        groupCellRows = find(groupCells == cellNames(iCell));
        
        allDataCell = horzcat(allData{1,iRunGroup}{1,groupCellRows})';

        [dataCount,~] = histcounts(allDataCell,6);
        binCenters = edges(2:end) - (edges(2)-edges(1))/2;
        
        meansAllDataCells = [meansAllDataCells; dataCount/numel(groupCellRows)];
        
    end
    
    meanAllDataCells = nanmean(meansAllDataCells,1)/0.1;
    stdAllDataCells = nanmean(meansAllDataCells,1);
    plot(binCenters, meanAllDataCells, 'Color', cMap(iRunGroup,:), 'LineWidth', 2)
    
    if numel(selectedCellsIDs) > 1
        %ploterrorbar = errorbar(binCenters,meanAllDataCells,stdAllDataCells,'.','Color',[0 0 0],'LineWidth', 1.2, 'CapSize', 14)
        %ploterrorbar.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
    
    hold on
    cMap(iRunGroup,:)
    
end

xlim([0,600])
ylabel('Firing frequency [Hz]')
box on;
title([queryOptoParameter,' vs ', queryRunFeature]);
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');
legend('-100 pA','-50 pA','0 pA','50 pA','100 pA','250 pA','500 pA','750 pA','1000 pA','Location', 'northeast', 'FontSize',12)
ylim([0,120])

plotName = ['IntP' '_' queryRunFeature,'_vs_', queryOptoParameter '_' experimentName];
plotFullPath = fullfile(savePlotPath, plotName);

if savePlots == 1
    saveas(gcf, plotFullPath, 'pdf')
    saveas(gcf, plotFullPath, 'fig')
end

%% amplitudeCurrentPulse vs areaVtPulse

queryOptoParameter =  'amplitudeCurrentPulse';
queryRunFeature = 'areaVtPulse';

selectedCellsIDs = [1,2,3];
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
plot(xValues, meanQueryRunFeatureDataAllCells,'Color',[0 0.4470 0.7410],'LineWidth', 2);
errorbar(xValues,meanQueryRunFeatureDataAllCells,stdQueryRunFeatureDataAllCells,'.','Color',[0 0 0],'LineWidth', 1.2, 'CapSize', 14)

hScatter = findobj(gca, 'Type', 'scatter');
uistack(hScatter, 'top');

xticks(xValues)
xticklabels([unique(queryOptoParameterValues)])
xlabel(queryOptoParameter);
ylabel(queryRunFeature);
box on;
title([queryOptoParameter,' vs ', queryRunFeature]);
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');

plotName = ['IntP' '_' queryRunFeature,'_vs_', queryOptoParameter '_' experimentName];
plotFullPath = fullfile(savePlotPath, plotName);
if savePlots == 1
    saveas(gcf, plotFullPath, 'pdf')
    saveas(gcf, plotFullPath, 'fig')
end
%% trace

queryOptoParameter =  'amplitudeCurrentPulse';
queryRunFeature = 'areaVtPulse';

selectedCellsIDs = [1,2,3];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0') & cell2mat(runGroupsParametersTable.widthCurrentPulse) == 6000);
selectedRunGroupsAnalysisTable = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroups),:);

trace1 = selectedRunGroupsAnalysisTable.trace{56,1};
trace1 = trace1(1:28048);
trace = selectedRunGroupsAnalysisTable.trace{48,1};
trace = trace(28000:end);

trace = [trace1, trace];
trace = trace(1:50000);

figure;
plot(trace,'Color',[0.7 0.7 0.7],'LineWidth',2)

xLimMin = 1.8;
xLimMax = 3.2;
xticks = linspace(xLimMin,xLimMax, 8);
xlim([xLimMin*10000, xLimMax*10000]);
ylim([-80, 60])
set(gca, 'XTickLabel', num2str((xticks'-xLimMin)))
xlabel('Time [s]');
ylabel('Membrane potential [mV]');
box on;
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');

plotName = ['IntP' '_', 'trace', '_' experimentName];
plotFullPath = fullfile(savePlotPath, plotName);

if savePlots == 1
    saveas(gcf, plotFullPath, 'pdf')
    saveas(gcf, plotFullPath, 'fig')
end
