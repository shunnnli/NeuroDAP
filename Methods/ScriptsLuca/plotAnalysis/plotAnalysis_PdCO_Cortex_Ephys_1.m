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
experimentName = 'PdCO_Cortex_Ephys_1_240409';

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

%% spotOnCell - heightPulsePeak

typePattern = 'spotOnCell';
queryRunFeature = 'heightPulsePeak';

selectedCellsIDs = [3;4];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

selectedRunGroupsRed = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao2') & cell2mat(runGroupsParametersTable.pulseWidthRed) == 20 & cell2mat(runGroupsParametersTable.amplitudeRed) == 5 & strcmp(runGroupsParametersTable.functionNameRed,'[onCellPattern, outsideCellPattern] = activationPatternForPdCO(3,1,0,200,20)'));
selectedRunGroupsAnalysisTableRed = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroupsRed),:);

selectedRunGroupsBlueRed = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1, ao2') & cell2mat(runGroupsParametersTable.pulseWidthRed) == 20 & cell2mat(runGroupsParametersTable.delayPulseBlue) == 28500 & cell2mat(runGroupsParametersTable.amplitudeBlue) == 5 & strcmp(runGroupsParametersTable.functionNameBlue,'[onCellPattern, outsideCellPattern] = activationPatternForPdCO(3,1,0,200,20)'));
selectedRunGroupsAnalysisTableBlueRed = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroupsBlueRed),:);

dataRed = vertcat(selectedRunGroupsAnalysisTableRed.(queryRunFeature){:});
dataBlueRed = vertcat(selectedRunGroupsAnalysisTableBlueRed.(queryRunFeature){:});

categoricalRed = zeros(numel(dataRed),1);
categoricalBlueRed = ones(numel(dataBlueRed),1);

allData = [dataRed; dataBlueRed]';
allCategoricalData = [categoricalRed; categoricalBlueRed]';

meanDataRed = nanmean(dataRed,1);
meanDataBlueRed = nanmean(dataBlueRed,1);
stdDataRed = nanstd(dataRed,1);
stdDataBlueRed = nanstd(dataBlueRed,1);

allMeansData = [meanDataRed; meanDataBlueRed]';
allStdsData = [stdDataRed; stdDataBlueRed]';
allCategoricalMeansStd = [0; 1]';

figure;

scatter(allCategoricalData, allData,'k','filled')
hold on;
bar(allCategoricalMeansStd, allMeansData,'FaceColor',[0 0.4470 0.7410],'BarWidth', 0.5);
errorbar(allCategoricalMeansStd, allMeansData, allStdsData,'.','Color',[0 0 0],'LineWidth', 1.2, 'CapSize', 14)
[h,p,ci,stats] = ttest2(dataRed', dataBlueRed');
text(1.2, 100, ['p-value: ', sprintf(['%0.', num2str(5), 'f'], p)], 'FontSize', 12, 'BackgroundColor', 'none', 'FontName', 'Arial');

hScatter = findobj(gca, 'Type', 'scatter');
uistack(hScatter, 'top');

xticks(allCategoricalMeansStd)
xticklabels({'ChRimson', 'ChRimson + PdCO'})
xlabel(['Optogenetic activation']);
ylabel(queryRunFeature);
box on;
title([typePattern ' - ' queryRunFeature]);
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');
ylim([-100,600])

plotName = [experimentName '_' queryRunFeature '_' typePattern];
plotFullPath = fullfile(savePlotPath, plotName);

if savePlots == 1
    saveas(gcf, plotFullPath, 'pdf')
    saveas(gcf, plotFullPath, 'fig')
end

%% fullField - heightPulsePeak

typePattern = 'fullField';
queryRunFeature = 'heightPulsePeak';

selectedCellsIDs = [3,4];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

selectedRunGroupsRed1 = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao2') & cell2mat(runGroupsParametersTable.nPulsesRed) == 1 & cell2mat(runGroupsParametersTable.pulseWidthRed) == 100 & cell2mat(runGroupsParametersTable.amplitudeRed) == 5 & strcmp(runGroupsParametersTable.functionNameRed,'dmdOutputRect(1, 1, 1, 1, 1, 1)'));
selectedRunGroupsRed2 = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao2') & cell2mat(runGroupsParametersTable.nPulsesRed) == 1 & cell2mat(runGroupsParametersTable.pulseWidthRed) == 100 & cell2mat(runGroupsParametersTable.amplitudeRed) == 5 & strcmp(runGroupsParametersTable.functionNameRed,'dmdOutputRect(1, 1, 1, 1, 2, 1)'));
selectedRunGroupsRed = [selectedRunGroupsRed1; selectedRunGroupsRed2];
selectedRunGroupsAnalysisTableRed = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroupsRed),:);

selectedRunGroupsBlueRed1 = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1, ao2') & cell2mat(runGroupsParametersTable.nPulsesRed) == 1 & cell2mat(runGroupsParametersTable.pulseWidthRed) == 100 & cell2mat(runGroupsParametersTable.delayPulseBlue) == 29500 & cell2mat(runGroupsParametersTable.amplitudeBlue) == 5 & strcmp(runGroupsParametersTable.functionNameBlue,'dmdOutputRect(1, 1, 1, 1, 1, 1)'));
selectedRunGroupsBlueRed2 = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1, ao2') & cell2mat(runGroupsParametersTable.nPulsesRed) == 1 & cell2mat(runGroupsParametersTable.pulseWidthRed) == 100 & cell2mat(runGroupsParametersTable.delayPulseBlue) == 29500 & cell2mat(runGroupsParametersTable.amplitudeBlue) == 5 & strcmp(runGroupsParametersTable.functionNameBlue,'dmdOutputRect(1, 1, 1, 1, 2, 1)'));
selectedRunGroupsBlueRed = [selectedRunGroupsBlueRed1; selectedRunGroupsBlueRed2];
selectedRunGroupsAnalysisTableBlueRed = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroupsBlueRed),:);

dataRed = vertcat(selectedRunGroupsAnalysisTableRed.(queryRunFeature){:});
dataBlueRed = vertcat(selectedRunGroupsAnalysisTableBlueRed.(queryRunFeature){:});

categoricalRed = zeros(numel(dataRed),1);
categoricalBlueRed = ones(numel(dataBlueRed),1);

allData = [dataRed; dataBlueRed]';
allCategoricalData = [categoricalRed; categoricalBlueRed]';

meanDataRed = nanmean(dataRed,1);
meanDataBlueRed = nanmean(dataBlueRed,1);
stdDataRed = nanstd(dataRed,1);
stdDataBlueRed = nanstd(dataBlueRed,1);

allMeansData = [meanDataRed; meanDataBlueRed]';
allStdsData = [stdDataRed; stdDataBlueRed]';
allCategoricalMeansStd = [0; 1]';

figure;

scatter(allCategoricalData, allData,'k','filled')
hold on;
bar(allCategoricalMeansStd, allMeansData,'FaceColor',[0 0.4470 0.7410],'BarWidth', 0.5);
errorbar(allCategoricalMeansStd, allMeansData, allStdsData,'.','Color',[0 0 0],'LineWidth', 1.2, 'CapSize', 14)
[h,p,ci,stats] = ttest2(dataRed', dataBlueRed');
text(1.2, 100, ['p-value: ', sprintf(['%0.', num2str(5), 'f'], p)], 'FontSize', 12, 'BackgroundColor', 'none', 'FontName', 'Arial');


hScatter = findobj(gca, 'Type', 'scatter');
uistack(hScatter, 'top');

xticks(allCategoricalMeansStd)
xticklabels({'ChRimson', 'ChRimson + PdCO'})
xlabel(['Optogenetic activation']);
ylabel(queryRunFeature);
box on;
title([typePattern ' - ' queryRunFeature]);
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');
ylim([-100,600])

plotName = [experimentName '_' queryRunFeature '_' typePattern];
plotFullPath = fullfile(savePlotPath, plotName);

if savePlots == 1
    saveas(gcf, plotFullPath, 'pdf')
    saveas(gcf, plotFullPath, 'fig')
end

%% spotOnCell - areaPulse

typePattern = 'spotOnCell';
queryRunFeature = 'areaPulse';

selectedCellsIDs = [3,4];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

selectedRunGroupsRed = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao2') & cell2mat(runGroupsParametersTable.pulseWidthRed) == 20 & cell2mat(runGroupsParametersTable.amplitudeRed) == 5 & strcmp(runGroupsParametersTable.functionNameRed,'[onCellPattern, outsideCellPattern] = activationPatternForPdCO(3,1,0,200,20)'));
selectedRunGroupsAnalysisTableRed = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroupsRed),:);

selectedRunGroupsBlueRed = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1, ao2') & cell2mat(runGroupsParametersTable.pulseWidthRed) == 20 & cell2mat(runGroupsParametersTable.delayPulseBlue) == 28500 & cell2mat(runGroupsParametersTable.amplitudeBlue) == 5 & strcmp(runGroupsParametersTable.functionNameBlue,'[onCellPattern, outsideCellPattern] = activationPatternForPdCO(3,1,0,200,20)'));
selectedRunGroupsAnalysisTableBlueRed = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroupsBlueRed),:);

dataRed = vertcat(selectedRunGroupsAnalysisTableRed.(queryRunFeature){:});
dataBlueRed = vertcat(selectedRunGroupsAnalysisTableBlueRed.(queryRunFeature){:});

categoricalRed = zeros(numel(dataRed),1);
categoricalBlueRed = ones(numel(dataBlueRed),1);

allData = [dataRed; dataBlueRed]';
allCategoricalData = [categoricalRed; categoricalBlueRed]';

meanDataRed = nanmean(dataRed,1);
meanDataBlueRed = nanmean(dataBlueRed,1);
stdDataRed = nanstd(dataRed,1);
stdDataBlueRed = nanstd(dataBlueRed,1);

allMeansData = [meanDataRed; meanDataBlueRed]';
allStdsData = [stdDataRed; stdDataBlueRed]';
allCategoricalMeansStd = [0; 1]';

figure;

scatter(allCategoricalData, allData,'k','filled')
hold on;
bar(allCategoricalMeansStd, allMeansData,'FaceColor',[0 0.4470 0.7410],'BarWidth', 0.5);
errorbar(allCategoricalMeansStd, allMeansData, allStdsData,'.','Color',[0 0 0],'LineWidth', 1.2, 'CapSize', 14)
[h,p,ci,stats] = ttest2(dataRed', dataBlueRed');
text(1.2, 2, ['p-value: ', sprintf(['%0.', num2str(5), 'f'], p)], 'FontSize', 12, 'BackgroundColor', 'none', 'FontName', 'Arial');

hScatter = findobj(gca, 'Type', 'scatter');
uistack(hScatter, 'top');

xticks(allCategoricalMeansStd)
xticklabels({'ChRimson', 'ChRimson + PdCO'})
xlabel(['Optogenetic activation']);
ylabel(queryRunFeature);
box on;
title([typePattern ' - ' queryRunFeature]);
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');
ylim([-1,3])

plotName = [experimentName '_' queryRunFeature '_' typePattern];
plotFullPath = fullfile(savePlotPath, plotName);

if savePlots == 1
    saveas(gcf, plotFullPath, 'pdf')
    saveas(gcf, plotFullPath, 'fig')
end

%% fullField - areaPulse

typePattern = 'fullField';
queryRunFeature = 'areaPulse';

selectedCellsIDs = [3];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

selectedRunGroupsRed1 = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao2') & cell2mat(runGroupsParametersTable.nPulsesRed) == 1 & cell2mat(runGroupsParametersTable.pulseWidthRed) == 100 & cell2mat(runGroupsParametersTable.amplitudeRed) == 5 & strcmp(runGroupsParametersTable.functionNameRed,'dmdOutputRect(1, 1, 1, 1, 1, 1)'));
selectedRunGroupsRed2 = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao2') & cell2mat(runGroupsParametersTable.nPulsesRed) == 1 & cell2mat(runGroupsParametersTable.pulseWidthRed) == 100 & cell2mat(runGroupsParametersTable.amplitudeRed) == 5 & strcmp(runGroupsParametersTable.functionNameRed,'dmdOutputRect(1, 1, 1, 1, 2, 1)'));
selectedRunGroupsRed = [selectedRunGroupsRed1; selectedRunGroupsRed2];
selectedRunGroupsAnalysisTableRed = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroupsRed),:);
notExcludedEpochs = find(cell2mat(selectedRunGroupsAnalysisTableRed.epoch) ~= 37 & cell2mat(selectedRunGroupsAnalysisTableRed.epoch) ~= 39 & cell2mat(selectedRunGroupsAnalysisTableRed.epoch) ~= 8);
selectedRunGroupsAnalysisTableRed = selectedRunGroupsAnalysisTableRed(notExcludedEpochs,:);

selectedRunGroupsBlueRed1 = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1, ao2') & cell2mat(runGroupsParametersTable.nPulsesRed) == 1 & cell2mat(runGroupsParametersTable.pulseWidthRed) == 100 & cell2mat(runGroupsParametersTable.delayPulseBlue) == 29500 & cell2mat(runGroupsParametersTable.amplitudeBlue) == 5 & strcmp(runGroupsParametersTable.functionNameBlue,'dmdOutputRect(1, 1, 1, 1, 1, 1)'));
selectedRunGroupsBlueRed2 = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1, ao2') & cell2mat(runGroupsParametersTable.nPulsesRed) == 1 & cell2mat(runGroupsParametersTable.pulseWidthRed) == 100 & cell2mat(runGroupsParametersTable.delayPulseBlue) == 29500 & cell2mat(runGroupsParametersTable.amplitudeBlue) == 5 & strcmp(runGroupsParametersTable.functionNameBlue,'dmdOutputRect(1, 1, 1, 1, 2, 1)'));
selectedRunGroupsBlueRed = [selectedRunGroupsBlueRed1; selectedRunGroupsBlueRed2];
selectedRunGroupsAnalysisTableBlueRed = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroupsBlueRed),:);
notExcludedEpochs = find(cell2mat(selectedRunGroupsAnalysisTableBlueRed.epoch) ~= 37 & cell2mat(selectedRunGroupsAnalysisTableBlueRed.epoch) ~= 39 & cell2mat(selectedRunGroupsAnalysisTableBlueRed.epoch) ~= 8);
selectedRunGroupsAnalysisTableBlueRed = selectedRunGroupsAnalysisTableBlueRed(notExcludedEpochs,:);

dataRed = vertcat(selectedRunGroupsAnalysisTableRed.(queryRunFeature){:});
dataBlueRed = vertcat(selectedRunGroupsAnalysisTableBlueRed.(queryRunFeature){:});

categoricalRed = zeros(numel(dataRed),1);
categoricalBlueRed = ones(numel(dataBlueRed),1);

allData = [dataRed; dataBlueRed]';
allCategoricalData = [categoricalRed; categoricalBlueRed]';

meanDataRed = nanmean(dataRed,1);
meanDataBlueRed = nanmean(dataBlueRed,1);
stdDataRed = nanstd(dataRed,1);
stdDataBlueRed = nanstd(dataBlueRed,1);

allMeansData = [meanDataRed; meanDataBlueRed]';
allStdsData = [stdDataRed; stdDataBlueRed]';
allCategoricalMeansStd = [0; 1]';

figure;

scatter(allCategoricalData, allData,'k','filled')
hold on;
bar(allCategoricalMeansStd, allMeansData,'FaceColor',[0 0.4470 0.7410],'BarWidth', 0.5);
errorbar(allCategoricalMeansStd, allMeansData, allStdsData,'.','Color',[0 0 0],'LineWidth', 1.2, 'CapSize', 14)
[h,p,ci,stats] = ttest2(dataRed', dataBlueRed');
text(1.2, 2, ['p-value: ', sprintf(['%0.', num2str(6), 'f'], p)], 'FontSize', 12, 'BackgroundColor', 'none', 'FontName', 'Arial');

hScatter = findobj(gca, 'Type', 'scatter');
uistack(hScatter, 'top');

xticks(allCategoricalMeansStd)
xticklabels({'ChRimson', 'ChRimson + PdCO'})
xlabel(['Optogenetic activation']);
ylabel(queryRunFeature);
box on;
title([typePattern ' - ' queryRunFeature]);
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');
ylim([-4,14])

plotName = [experimentName '_' queryRunFeature '_' typePattern];
plotFullPath = fullfile(savePlotPath, plotName);

if savePlots == 1
    saveas(gcf, plotFullPath, 'pdf')
    saveas(gcf, plotFullPath, 'fig')
end

%% fullField - with PdCO only - heightPulsePeak

typePattern = 'fullField';
queryRunFeature = 'heightPulsePeak';

selectedCellsIDs = [4];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

selectedRunGroupsRed1 = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao2') & cell2mat(runGroupsParametersTable.nPulsesRed) == 1 & cell2mat(runGroupsParametersTable.pulseWidthRed) == 100 & cell2mat(runGroupsParametersTable.amplitudeRed) == 5 & strcmp(runGroupsParametersTable.functionNameRed,'dmdOutputRect(1, 1, 1, 1, 1, 1)'));
selectedRunGroupsRed2 = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao2') & cell2mat(runGroupsParametersTable.nPulsesRed) == 1 & cell2mat(runGroupsParametersTable.pulseWidthRed) == 100 & cell2mat(runGroupsParametersTable.amplitudeRed) == 5 & strcmp(runGroupsParametersTable.functionNameRed,'dmdOutputRect(1, 1, 1, 1, 2, 1)'));
selectedRunGroupsRed = [selectedRunGroupsRed1; selectedRunGroupsRed2];
selectedRunGroupsAnalysisTableRed = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroupsRed),:);
notExcludedEpochs = find(cell2mat(selectedRunGroupsAnalysisTableRed.epoch) ~= 37 & cell2mat(selectedRunGroupsAnalysisTableRed.epoch) ~= 39 & cell2mat(selectedRunGroupsAnalysisTableRed.epoch) ~= 8);
selectedRunGroupsAnalysisTableRed = selectedRunGroupsAnalysisTableRed(notExcludedEpochs,:);

selectedRunGroupsBlueRed1 = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1, ao2') & cell2mat(runGroupsParametersTable.nPulsesRed) == 1 & cell2mat(runGroupsParametersTable.pulseWidthRed) == 100 & cell2mat(runGroupsParametersTable.delayPulseBlue) == 29500 & cell2mat(runGroupsParametersTable.amplitudeBlue) == 5 & strcmp(runGroupsParametersTable.functionNameBlue,'dmdOutputRect(1, 1, 1, 1, 1, 1)'));
selectedRunGroupsBlueRed2 = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1, ao2') & cell2mat(runGroupsParametersTable.nPulsesRed) == 1 & cell2mat(runGroupsParametersTable.pulseWidthRed) == 100 & cell2mat(runGroupsParametersTable.delayPulseBlue) == 29500 & cell2mat(runGroupsParametersTable.amplitudeBlue) == 5 & strcmp(runGroupsParametersTable.functionNameBlue,'dmdOutputRect(1, 1, 1, 1, 2, 1)'));
selectedRunGroupsBlueRed = [selectedRunGroupsBlueRed1; selectedRunGroupsBlueRed2];
selectedRunGroupsAnalysisTableBlueRed = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroupsBlueRed),:);
notExcludedEpochs = find(cell2mat(selectedRunGroupsAnalysisTableBlueRed.epoch) ~= 37 & cell2mat(selectedRunGroupsAnalysisTableBlueRed.epoch) ~= 39 & cell2mat(selectedRunGroupsAnalysisTableBlueRed.epoch) ~= 8);
selectedRunGroupsAnalysisTableBlueRed = selectedRunGroupsAnalysisTableBlueRed(notExcludedEpochs,:);

selectedCellsIDs = [3];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

selectedRunGroupsBlue = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.pulseWidthBlue) == 100);
selectedRunGroupsAnalysisTableBlue = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroupsBlue),:);
notExcludedEpochs = find(cell2mat(selectedRunGroupsAnalysisTableBlue.epoch) ~= 37 & cell2mat(selectedRunGroupsAnalysisTableBlue.epoch) ~= 39 & cell2mat(selectedRunGroupsAnalysisTableBlue.epoch) ~= 8);
selectedRunGroupsAnalysisTableBlue = selectedRunGroupsAnalysisTableBlue(notExcludedEpochs,:);

dataRed = vertcat(selectedRunGroupsAnalysisTableRed.(queryRunFeature){:});
dataBlueRed = vertcat(selectedRunGroupsAnalysisTableBlueRed.(queryRunFeature){:});
dataBlue = vertcat(selectedRunGroupsAnalysisTableBlue.(queryRunFeature){:});

categoricalRed = zeros(numel(dataRed),1);
categoricalBlueRed = 2*ones(numel(dataBlueRed),1);
categoricalBlue = ones(numel(dataBlue),1);

allData = [dataRed; dataBlue; dataBlueRed]';
allCategoricalData = [categoricalRed; categoricalBlue; categoricalBlueRed]';

meanDataRed = nanmean(dataRed,1);
meanDataBlueRed = nanmean(dataBlueRed,1);
meanDataBlue = nanmean(dataBlue,1);
stdDataRed = nanstd(dataRed,1);
stdDataBlueRed = nanstd(dataBlueRed,1);
stdDataBlue = nanstd(dataBlue,1);

allMeansData = [meanDataRed; meanDataBlue; meanDataBlueRed]';
allStdsData = [stdDataRed; stdDataBlue; stdDataBlueRed]';
allCategoricalMeansStd = [0; 1; 2]';

figure;

scatter(allCategoricalData, allData,'k','filled')
hold on;
bar(allCategoricalMeansStd, allMeansData,'FaceColor',[0 0.4470 0.7410],'BarWidth', 0.5);
errorbar(allCategoricalMeansStd, allMeansData, allStdsData,'.','Color',[0 0 0],'LineWidth', 1.2, 'CapSize', 14)
[h,p,ci,stats] = ttest2(dataRed', dataBlueRed');
%text(1.2, 80, ['p-value: ', sprintf(['%0.', num2str(5), 'f'], p)], 'FontSize', 12, 'BackgroundColor', 'none', 'FontName', 'Arial');

hScatter = findobj(gca, 'Type', 'scatter');
uistack(hScatter, 'top');

xticks(allCategoricalMeansStd)
xticklabels({'ChRimson', 'PdCO', 'ChRimson + PdCO'})
xlabel(['Optogenetic activation']);
ylabel(queryRunFeature);
box on;
title([typePattern ' - ' queryRunFeature]);
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');
ylim([-20,100])

plotName = [experimentName '_' queryRunFeature '_' typePattern '_' 'withPdCOonly'];
plotFullPath = fullfile(savePlotPath, plotName);

if savePlots == 1
    saveas(gcf, plotFullPath, 'pdf')
    saveas(gcf, plotFullPath, 'fig')
end

%% fullField - with PdCO only - areaPulse

typePattern = 'fullField';
queryRunFeature = 'areaPulse';

selectedCellsIDs = [4];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

selectedRunGroupsRed1 = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao2') & cell2mat(runGroupsParametersTable.nPulsesRed) == 1 & cell2mat(runGroupsParametersTable.pulseWidthRed) == 100 & cell2mat(runGroupsParametersTable.amplitudeRed) == 5 & strcmp(runGroupsParametersTable.functionNameRed,'dmdOutputRect(1, 1, 1, 1, 1, 1)'));
selectedRunGroupsRed2 = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao2') & cell2mat(runGroupsParametersTable.nPulsesRed) == 1 & cell2mat(runGroupsParametersTable.pulseWidthRed) == 100 & cell2mat(runGroupsParametersTable.amplitudeRed) == 5 & strcmp(runGroupsParametersTable.functionNameRed,'dmdOutputRect(1, 1, 1, 1, 2, 1)'));
selectedRunGroupsRed = [selectedRunGroupsRed1; selectedRunGroupsRed2];
selectedRunGroupsAnalysisTableRed = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroupsRed),:);
notExcludedEpochs = find(cell2mat(selectedRunGroupsAnalysisTableRed.epoch) ~= 37 & cell2mat(selectedRunGroupsAnalysisTableRed.epoch) ~= 8);
selectedRunGroupsAnalysisTableRed = selectedRunGroupsAnalysisTableRed(notExcludedEpochs,:);

selectedRunGroupsBlueRed1 = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1, ao2') & cell2mat(runGroupsParametersTable.nPulsesRed) == 1 & cell2mat(runGroupsParametersTable.pulseWidthRed) == 100 & cell2mat(runGroupsParametersTable.delayPulseBlue) == 29500 & cell2mat(runGroupsParametersTable.amplitudeBlue) == 5 & strcmp(runGroupsParametersTable.functionNameBlue,'dmdOutputRect(1, 1, 1, 1, 1, 1)'));
selectedRunGroupsBlueRed2 = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1, ao2') & cell2mat(runGroupsParametersTable.nPulsesRed) == 1 & cell2mat(runGroupsParametersTable.pulseWidthRed) == 100 & cell2mat(runGroupsParametersTable.delayPulseBlue) == 29500 & cell2mat(runGroupsParametersTable.amplitudeBlue) == 5 & strcmp(runGroupsParametersTable.functionNameBlue,'dmdOutputRect(1, 1, 1, 1, 2, 1)'));
selectedRunGroupsBlueRed = [selectedRunGroupsBlueRed1; selectedRunGroupsBlueRed2];
selectedRunGroupsAnalysisTableBlueRed = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroupsBlueRed),:);
notExcludedEpochs = find(cell2mat(selectedRunGroupsAnalysisTableBlueRed.epoch) ~= 37 & cell2mat(selectedRunGroupsAnalysisTableBlueRed.epoch) ~= 8);
selectedRunGroupsAnalysisTableBlueRed = selectedRunGroupsAnalysisTableBlueRed(notExcludedEpochs,:);

selectedCellsIDs = [3];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

selectedRunGroupsBlue = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.pulseWidthBlue) == 100);
selectedRunGroupsAnalysisTableBlue = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroupsBlue),:);
notExcludedEpochs = find(cell2mat(selectedRunGroupsAnalysisTableBlue.epoch) ~= 37 & cell2mat(selectedRunGroupsAnalysisTableBlue.epoch) ~= 8);
selectedRunGroupsAnalysisTableBlue = selectedRunGroupsAnalysisTableBlue(notExcludedEpochs,:);


dataRed = vertcat(selectedRunGroupsAnalysisTableRed.(queryRunFeature){:});
dataBlueRed = vertcat(selectedRunGroupsAnalysisTableBlueRed.(queryRunFeature){:});
dataBlue = vertcat(selectedRunGroupsAnalysisTableBlue.(queryRunFeature){:});

categoricalRed = zeros(numel(dataRed),1);
categoricalBlueRed = 2*ones(numel(dataBlueRed),1);
categoricalBlue = ones(numel(dataBlue),1);

allData = [dataRed; dataBlue; dataBlueRed]';
allCategoricalData = [categoricalRed; categoricalBlue; categoricalBlueRed]';

meanDataRed = nanmean(dataRed,1);
meanDataBlueRed = nanmean(dataBlueRed,1);
meanDataBlue = nanmean(dataBlue,1);
stdDataRed = nanstd(dataRed,1);
stdDataBlueRed = nanstd(dataBlueRed,1);
stdDataBlue = nanstd(dataBlue,1);

allMeansData = [meanDataRed; meanDataBlue; meanDataBlueRed]';
allStdsData = [stdDataRed; stdDataBlue; stdDataBlueRed]';
allCategoricalMeansStd = [0; 1; 2]';

figure;

scatter(allCategoricalData, allData,'k','filled')
hold on;
bar(allCategoricalMeansStd, allMeansData,'FaceColor',[0 0.4470 0.7410],'BarWidth', 0.5);
errorbar(allCategoricalMeansStd, allMeansData, allStdsData,'.','Color',[0 0 0],'LineWidth', 1.2, 'CapSize', 14)
[h,p,ci,stats] = ttest2(dataRed', dataBlueRed');
text(1.2, 100, ['p-value: ', sprintf(['%0.', num2str(5), 'f'], p)], 'FontSize', 12, 'BackgroundColor', 'none', 'FontName', 'Arial');


hScatter = findobj(gca, 'Type', 'scatter');
uistack(hScatter, 'top');

xticks(allCategoricalMeansStd)
xticklabels({'ChRimson', 'PdCO', 'ChRimson + PdCO'})
xlabel(['Optogenetic activation']);
ylabel(queryRunFeature);
box on;
title([typePattern ' - ' queryRunFeature]);
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');
ylim([-1,3])

plotName = [experimentName '_' queryRunFeature '_' typePattern '_' 'withPdCOonly'];
plotFullPath = fullfile(savePlotPath, plotName);

if savePlots == 1
    saveas(gcf, plotFullPath, 'pdf')
    saveas(gcf, plotFullPath, 'fig')
end

%% spotOnCell - heightPulsePeak - before and after

typePattern = 'spotOnCell';
queryRunFeature = 'heightPulsePeak';

selectedCellsIDs = [4];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

selectedRunGroupsRed = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao2') & cell2mat(runGroupsParametersTable.pulseWidthRed) == 20 & cell2mat(runGroupsParametersTable.amplitudeRed) == 5 & strcmp(runGroupsParametersTable.functionNameRed,'[onCellPattern, outsideCellPattern] = activationPatternForPdCO(3,1,0,200,20)'));
selectedRunGroupsAnalysisTableRed = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroupsRed),:);
selectedCyclePositionBefore = find(cell2mat(selectedRunGroupsAnalysisTableRed.cyclePosition)==1);
selectedCyclePositionAfter = find(cell2mat(selectedRunGroupsAnalysisTableRed.cyclePosition)==3);
selectedRunGroupsAnalysisTableRedBefore = selectedRunGroupsAnalysisTableRed(selectedCyclePositionBefore,:);
selectedRunGroupsAnalysisTableRedAfter = selectedRunGroupsAnalysisTableRed(selectedCyclePositionAfter,:);

selectedRunGroupsBlueRed = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1, ao2') & cell2mat(runGroupsParametersTable.pulseWidthRed) == 20 & cell2mat(runGroupsParametersTable.delayPulseBlue) == 28500 & cell2mat(runGroupsParametersTable.amplitudeBlue) == 5 & strcmp(runGroupsParametersTable.functionNameBlue,'[onCellPattern, outsideCellPattern] = activationPatternForPdCO(3,1,0,200,20)'));
selectedRunGroupsAnalysisTableBlueRed = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroupsBlueRed),:);

dataRedBefore = vertcat(selectedRunGroupsAnalysisTableRedBefore.(queryRunFeature){:});
dataRedAfter = vertcat(selectedRunGroupsAnalysisTableRedAfter.(queryRunFeature){:});
dataBlueRed = vertcat(selectedRunGroupsAnalysisTableBlueRed.(queryRunFeature){:});

categoricalRedBefore = zeros(numel(dataRedBefore),1);
categoricalRedAfter = 2*ones(numel(dataRedAfter),1);
categoricalBlueRed = ones(numel(dataBlueRed),1);

allData = [dataRedBefore; dataBlueRed; dataRedAfter]';
allCategoricalData = [categoricalRedBefore; categoricalBlueRed; categoricalRedAfter]';

meanDataRedBefore = nanmean(dataRedBefore,1);
meanDataRedAfter = nanmean(dataRedAfter,1);
meanDataBlueRed = nanmean(dataBlueRed,1);
stdDataRedBefore = nanstd(dataRedBefore,1);
stdDataRedAfter = nanstd(dataRedAfter,1);
stdDataBlueRed = nanstd(dataBlueRed,1);

allMeansData = [meanDataRedBefore; meanDataBlueRed; meanDataRedAfter]';
allStdsData = [stdDataRedBefore; stdDataBlueRed; stdDataRedAfter]';
allCategoricalMeansStd = [0; 1; 2]';

figure;

for iPair = 1:numel(dataRedBefore)

    scatter([allCategoricalData(iPair) allCategoricalData(iPair+1*numel(dataRedBefore)) allCategoricalData(iPair+2*numel(dataRedBefore))], [allData(iPair) allData(iPair+1*numel(dataRedBefore)) allData(iPair+2*numel(dataRedBefore))],'k','filled')
    hold on;
    plot([allCategoricalData(iPair) allCategoricalData(iPair+1*numel(dataRedBefore)) allCategoricalData(iPair+2*numel(dataRedBefore))], [allData(iPair) allData(iPair+1*numel(dataRedBefore)) allData(iPair+2*numel(dataRedBefore))],'k')
    hold on;
    
end

bar(allCategoricalMeansStd, allMeansData,'FaceColor',[0 0.4470 0.7410],'BarWidth', 0.5);
errorbar(allCategoricalMeansStd, allMeansData, allStdsData,'.','Color',[0 0 0],'LineWidth', 1.2, 'CapSize', 14)
%[h,p,ci,stats] = ttest2(dataRed', dataBlueRed');
%text(1.2, 100, ['p-value: ', sprintf(['%0.', num2str(5), 'f'], p)], 'FontSize', 12, 'BackgroundColor', 'none', 'FontName', 'Arial');

hScatter = findobj(gca, 'Type', 'scatter');
hPlot = findobj(gca, 'Type', 'plot');
hBar = findobj(gca, 'Type', 'bar');
uistack(hScatter, 'top');
uistack(hPlot, 'top');
uistack(hBar, 'bottom');

xticks(allCategoricalMeansStd)
xticklabels({'ChRimson Before', 'ChRimson + PdCO', 'ChRimson After'})
xlabel(['Optogenetic activation']);
ylabel(queryRunFeature);
box on;
title([typePattern ' - ' queryRunFeature]);
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');
ylim([-20,100])

plotName = [experimentName '_' queryRunFeature '_' typePattern '_' 'beforeAfter'];
plotFullPath = fullfile(savePlotPath, plotName);

if savePlots == 1
    saveas(gcf, plotFullPath, 'pdf')
    saveas(gcf, plotFullPath, 'fig')
end

%% fullField - heightPulsePeak - before and after

typePattern = 'fullField';
queryRunFeature = 'heightPulsePeak';

selectedCellsIDs = [4];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

selectedRunGroupsRed1 = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao2') & cell2mat(runGroupsParametersTable.nPulsesRed) == 1 & cell2mat(runGroupsParametersTable.pulseWidthRed) == 100 & cell2mat(runGroupsParametersTable.amplitudeRed) == 5 & strcmp(runGroupsParametersTable.functionNameRed,'dmdOutputRect(1, 1, 1, 1, 1, 1)'));
selectedRunGroupsRed2 = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao2') & cell2mat(runGroupsParametersTable.nPulsesRed) == 1 & cell2mat(runGroupsParametersTable.pulseWidthRed) == 100 & cell2mat(runGroupsParametersTable.amplitudeRed) == 5 & strcmp(runGroupsParametersTable.functionNameRed,'dmdOutputRect(1, 1, 1, 1, 2, 1)'));
selectedRunGroupsRed = [selectedRunGroupsRed1; selectedRunGroupsRed2];
selectedRunGroupsAnalysisTableRed = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroupsRed),:);
selectedCyclePositionBefore = find(cell2mat(selectedRunGroupsAnalysisTableRed.cyclePosition)==1);
selectedCyclePositionAfter = find(cell2mat(selectedRunGroupsAnalysisTableRed.cyclePosition)==3);
selectedRunGroupsAnalysisTableRedBefore = selectedRunGroupsAnalysisTableRed(selectedCyclePositionBefore,:);
selectedRunGroupsAnalysisTableRedAfter = selectedRunGroupsAnalysisTableRed(selectedCyclePositionAfter,:);

selectedRunGroupsBlueRed1 = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1, ao2') & cell2mat(runGroupsParametersTable.nPulsesRed) == 1 & cell2mat(runGroupsParametersTable.pulseWidthRed) == 100 & cell2mat(runGroupsParametersTable.delayPulseBlue) == 29500 & cell2mat(runGroupsParametersTable.amplitudeBlue) == 5 & strcmp(runGroupsParametersTable.functionNameBlue,'dmdOutputRect(1, 1, 1, 1, 1, 1)'));
selectedRunGroupsBlueRed2 = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1, ao2') & cell2mat(runGroupsParametersTable.nPulsesRed) == 1 & cell2mat(runGroupsParametersTable.pulseWidthRed) == 100 & cell2mat(runGroupsParametersTable.delayPulseBlue) == 29500 & cell2mat(runGroupsParametersTable.amplitudeBlue) == 5 & strcmp(runGroupsParametersTable.functionNameBlue,'dmdOutputRect(1, 1, 1, 1, 2, 1)'));
selectedRunGroupsBlueRed = [selectedRunGroupsBlueRed1; selectedRunGroupsBlueRed2];
selectedRunGroupsAnalysisTableBlueRed = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroupsBlueRed),:);

dataRedBefore = vertcat(selectedRunGroupsAnalysisTableRedBefore.(queryRunFeature){:});
dataRedAfter = vertcat(selectedRunGroupsAnalysisTableRedAfter.(queryRunFeature){:});
dataBlueRed = vertcat(selectedRunGroupsAnalysisTableBlueRed.(queryRunFeature){:});

categoricalRedBefore = zeros(numel(dataRedBefore),1);
categoricalRedAfter = 2*ones(numel(dataRedAfter),1);
categoricalBlueRed = ones(numel(dataBlueRed),1);

allData = [dataRedBefore; dataBlueRed; dataRedAfter]';
allCategoricalData = [categoricalRedBefore; categoricalBlueRed; categoricalRedAfter]';

meanDataRedBefore = nanmean(dataRedBefore,1);
meanDataRedAfter = nanmean(dataRedAfter,1);
meanDataBlueRed = nanmean(dataBlueRed,1);
stdDataRedBefore = nanstd(dataRedBefore,1);
stdDataRedAfter = nanstd(dataRedAfter,1);
stdDataBlueRed = nanstd(dataBlueRed,1);

allMeansData = [meanDataRedBefore; meanDataBlueRed; meanDataRedAfter]';
allStdsData = [stdDataRedBefore; stdDataBlueRed; stdDataRedAfter]';
allCategoricalMeansStd = [0; 1; 2]';

figure;

for iPair = 1:numel(dataRedBefore)

    scatter([allCategoricalData(iPair) allCategoricalData(iPair+1*numel(dataRedBefore)) allCategoricalData(iPair+2*numel(dataRedBefore))], [allData(iPair) allData(iPair+1*numel(dataRedBefore)) allData(iPair+2*numel(dataRedBefore))],'k','filled')
    hold on;
    plot([allCategoricalData(iPair) allCategoricalData(iPair+1*numel(dataRedBefore)) allCategoricalData(iPair+2*numel(dataRedBefore))], [allData(iPair) allData(iPair+1*numel(dataRedBefore)) allData(iPair+2*numel(dataRedBefore))],'k')
    hold on;
    
end

bar(allCategoricalMeansStd, allMeansData,'FaceColor',[0 0.4470 0.7410],'BarWidth', 0.5);
errorbar(allCategoricalMeansStd, allMeansData, allStdsData,'.','Color',[0 0 0],'LineWidth', 1.2, 'CapSize', 14)
%[h,p,ci,stats] = ttest2(dataRed', dataBlueRed');
%text(1.2, 100, ['p-value: ', sprintf(['%0.', num2str(5), 'f'], p)], 'FontSize', 12, 'BackgroundColor', 'none', 'FontName', 'Arial');

hScatter = findobj(gca, 'Type', 'scatter');
hPlot = findobj(gca, 'Type', 'plot');
hBar = findobj(gca, 'Type', 'bar');
uistack(hScatter, 'top');
uistack(hPlot, 'top');
uistack(hBar, 'bottom');

xticks(allCategoricalMeansStd)
xticklabels({'ChRimson Before', 'ChRimson + PdCO', 'ChRimson After'})
xlabel(['Optogenetic activation']);
ylabel(queryRunFeature);
box on;
title([typePattern ' - ' queryRunFeature]);
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');
ylim([-20,100])

plotName = [experimentName '_' queryRunFeature '_' typePattern '_' 'beforeAfter'];
plotFullPath = fullfile(savePlotPath, plotName);

if savePlots == 1
    saveas(gcf, plotFullPath, 'pdf')
    saveas(gcf, plotFullPath, 'fig')
end

%% spotOnCell - areaPulse - before and after

typePattern = 'spotOnCell';
queryRunFeature = 'areaPulse';

selectedCellsIDs = [4];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

selectedRunGroupsRed = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao2') & cell2mat(runGroupsParametersTable.pulseWidthRed) == 20 & cell2mat(runGroupsParametersTable.amplitudeRed) == 5 & strcmp(runGroupsParametersTable.functionNameRed,'[onCellPattern, outsideCellPattern] = activationPatternForPdCO(3,1,0,200,20)'));
selectedRunGroupsAnalysisTableRed = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroupsRed),:);
selectedCyclePositionBefore = find(cell2mat(selectedRunGroupsAnalysisTableRed.cyclePosition)==1);
selectedCyclePositionAfter = find(cell2mat(selectedRunGroupsAnalysisTableRed.cyclePosition)==3);
selectedRunGroupsAnalysisTableRedBefore = selectedRunGroupsAnalysisTableRed(selectedCyclePositionBefore,:);
selectedRunGroupsAnalysisTableRedAfter = selectedRunGroupsAnalysisTableRed(selectedCyclePositionAfter,:);

selectedRunGroupsBlueRed = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1, ao2') & cell2mat(runGroupsParametersTable.pulseWidthRed) == 20 & cell2mat(runGroupsParametersTable.delayPulseBlue) == 28500 & cell2mat(runGroupsParametersTable.amplitudeBlue) == 5 & strcmp(runGroupsParametersTable.functionNameBlue,'[onCellPattern, outsideCellPattern] = activationPatternForPdCO(3,1,0,200,20)'));
selectedRunGroupsAnalysisTableBlueRed = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroupsBlueRed),:);

dataRedBefore = vertcat(selectedRunGroupsAnalysisTableRedBefore.(queryRunFeature){:});
dataRedAfter = vertcat(selectedRunGroupsAnalysisTableRedAfter.(queryRunFeature){:});
dataBlueRed = vertcat(selectedRunGroupsAnalysisTableBlueRed.(queryRunFeature){:});

categoricalRedBefore = zeros(numel(dataRedBefore),1);
categoricalRedAfter = 2*ones(numel(dataRedAfter),1);
categoricalBlueRed = ones(numel(dataBlueRed),1);

allData = [dataRedBefore; dataBlueRed; dataRedAfter]';
allCategoricalData = [categoricalRedBefore; categoricalBlueRed; categoricalRedAfter]';

meanDataRedBefore = nanmean(dataRedBefore,1);
meanDataRedAfter = nanmean(dataRedAfter,1);
meanDataBlueRed = nanmean(dataBlueRed,1);
stdDataRedBefore = nanstd(dataRedBefore,1);
stdDataRedAfter = nanstd(dataRedAfter,1);
stdDataBlueRed = nanstd(dataBlueRed,1);

allMeansData = [meanDataRedBefore; meanDataBlueRed; meanDataRedAfter]';
allStdsData = [stdDataRedBefore; stdDataBlueRed; stdDataRedAfter]';
allCategoricalMeansStd = [0; 1; 2]';

figure;

for iPair = 1:numel(dataRedBefore)

    scatter([allCategoricalData(iPair) allCategoricalData(iPair+1*numel(dataRedBefore)) allCategoricalData(iPair+2*numel(dataRedBefore))], [allData(iPair) allData(iPair+1*numel(dataRedBefore)) allData(iPair+2*numel(dataRedBefore))],'k','filled')
    hold on;
    plot([allCategoricalData(iPair) allCategoricalData(iPair+1*numel(dataRedBefore)) allCategoricalData(iPair+2*numel(dataRedBefore))], [allData(iPair) allData(iPair+1*numel(dataRedBefore)) allData(iPair+2*numel(dataRedBefore))],'k')
    hold on;
    
end

bar(allCategoricalMeansStd, allMeansData,'FaceColor',[0 0.4470 0.7410],'BarWidth', 0.5);
errorbar(allCategoricalMeansStd, allMeansData, allStdsData,'.','Color',[0 0 0],'LineWidth', 1.2, 'CapSize', 14)
%[h,p,ci,stats] = ttest2(dataRed', dataBlueRed');
%text(1.2, 100, ['p-value: ', sprintf(['%0.', num2str(5), 'f'], p)], 'FontSize', 12, 'BackgroundColor', 'none', 'FontName', 'Arial');

hScatter = findobj(gca, 'Type', 'scatter');
hPlot = findobj(gca, 'Type', 'plot');
hBar = findobj(gca, 'Type', 'bar');
uistack(hScatter, 'top');
uistack(hPlot, 'top');
uistack(hBar, 'bottom');

xticks(allCategoricalMeansStd)
xticklabels({'ChRimson Before', 'ChRimson + PdCO', 'ChRimson After'})
xlabel(['Optogenetic activation']);
ylabel(queryRunFeature);
box on;
title([typePattern ' - ' queryRunFeature]);
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');
ylim([-1,3])

plotName = [experimentName '_' queryRunFeature '_' typePattern '_' 'beforeAfter'];
plotFullPath = fullfile(savePlotPath, plotName);

if savePlots == 1
    saveas(gcf, plotFullPath, 'pdf')
    saveas(gcf, plotFullPath, 'fig')
end

%% fullField - areaPulse - before and after

typePattern = 'fullField';
queryRunFeature = 'areaPulse';

selectedCellsIDs = [4];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

selectedRunGroupsRed1 = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao2') & cell2mat(runGroupsParametersTable.nPulsesRed) == 1 & cell2mat(runGroupsParametersTable.pulseWidthRed) == 100 & cell2mat(runGroupsParametersTable.amplitudeRed) == 5 & strcmp(runGroupsParametersTable.functionNameRed,'dmdOutputRect(1, 1, 1, 1, 1, 1)'));
selectedRunGroupsRed2 = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao2') & cell2mat(runGroupsParametersTable.nPulsesRed) == 1 & cell2mat(runGroupsParametersTable.pulseWidthRed) == 100 & cell2mat(runGroupsParametersTable.amplitudeRed) == 5 & strcmp(runGroupsParametersTable.functionNameRed,'dmdOutputRect(1, 1, 1, 1, 2, 1)'));
selectedRunGroupsRed = [selectedRunGroupsRed1; selectedRunGroupsRed2];
selectedRunGroupsAnalysisTableRed = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroupsRed),:);
selectedCyclePositionBefore = find(cell2mat(selectedRunGroupsAnalysisTableRed.cyclePosition)==1);
selectedCyclePositionAfter = find(cell2mat(selectedRunGroupsAnalysisTableRed.cyclePosition)==3);
selectedRunGroupsAnalysisTableRedBefore = selectedRunGroupsAnalysisTableRed(selectedCyclePositionBefore,:);
selectedRunGroupsAnalysisTableRedAfter = selectedRunGroupsAnalysisTableRed(selectedCyclePositionAfter,:);

selectedRunGroupsBlueRed1 = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1, ao2') & cell2mat(runGroupsParametersTable.nPulsesRed) == 1 & cell2mat(runGroupsParametersTable.pulseWidthRed) == 100 & cell2mat(runGroupsParametersTable.delayPulseBlue) == 29500 & cell2mat(runGroupsParametersTable.amplitudeBlue) == 5 & strcmp(runGroupsParametersTable.functionNameBlue,'dmdOutputRect(1, 1, 1, 1, 1, 1)'));
selectedRunGroupsBlueRed2 = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1, ao2') & cell2mat(runGroupsParametersTable.nPulsesRed) == 1 & cell2mat(runGroupsParametersTable.pulseWidthRed) == 100 & cell2mat(runGroupsParametersTable.delayPulseBlue) == 29500 & cell2mat(runGroupsParametersTable.amplitudeBlue) == 5 & strcmp(runGroupsParametersTable.functionNameBlue,'dmdOutputRect(1, 1, 1, 1, 2, 1)'));
selectedRunGroupsBlueRed = [selectedRunGroupsBlueRed1; selectedRunGroupsBlueRed2];
selectedRunGroupsAnalysisTableBlueRed = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroupsBlueRed),:);

dataRedBefore = vertcat(selectedRunGroupsAnalysisTableRedBefore.(queryRunFeature){:});
dataRedAfter = vertcat(selectedRunGroupsAnalysisTableRedAfter.(queryRunFeature){:});
dataBlueRed = vertcat(selectedRunGroupsAnalysisTableBlueRed.(queryRunFeature){:});

categoricalRedBefore = zeros(numel(dataRedBefore),1);
categoricalRedAfter = 2*ones(numel(dataRedAfter),1);
categoricalBlueRed = ones(numel(dataBlueRed),1);

allData = [dataRedBefore; dataBlueRed; dataRedAfter]';
allCategoricalData = [categoricalRedBefore; categoricalBlueRed; categoricalRedAfter]';

meanDataRedBefore = nanmean(dataRedBefore,1);
meanDataRedAfter = nanmean(dataRedAfter,1);
meanDataBlueRed = nanmean(dataBlueRed,1);
stdDataRedBefore = nanstd(dataRedBefore,1);
stdDataRedAfter = nanstd(dataRedAfter,1);
stdDataBlueRed = nanstd(dataBlueRed,1);

allMeansData = [meanDataRedBefore; meanDataBlueRed; meanDataRedAfter]';
allStdsData = [stdDataRedBefore; stdDataBlueRed; stdDataRedAfter]';
allCategoricalMeansStd = [0; 1; 2]';

figure;

for iPair = 1:numel(dataRedBefore)

    scatter([allCategoricalData(iPair) allCategoricalData(iPair+1*numel(dataRedBefore)) allCategoricalData(iPair+2*numel(dataRedBefore))], [allData(iPair) allData(iPair+1*numel(dataRedBefore)) allData(iPair+2*numel(dataRedBefore))],'k','filled')
    hold on;
    plot([allCategoricalData(iPair) allCategoricalData(iPair+1*numel(dataRedBefore)) allCategoricalData(iPair+2*numel(dataRedBefore))], [allData(iPair) allData(iPair+1*numel(dataRedBefore)) allData(iPair+2*numel(dataRedBefore))],'k')
    hold on;
    
end

bar(allCategoricalMeansStd, allMeansData,'FaceColor',[0 0.4470 0.7410],'BarWidth', 0.5);
errorbar(allCategoricalMeansStd, allMeansData, allStdsData,'.','Color',[0 0 0],'LineWidth', 1.2, 'CapSize', 14)
%[h,p,ci,stats] = ttest2(dataRed', dataBlueRed');
%text(1.2, 100, ['p-value: ', sprintf(['%0.', num2str(5), 'f'], p)], 'FontSize', 12, 'BackgroundColor', 'none', 'FontName', 'Arial');

hScatter = findobj(gca, 'Type', 'scatter');
hPlot = findobj(gca, 'Type', 'plot');
hBar = findobj(gca, 'Type', 'bar');
uistack(hScatter, 'top');
uistack(hPlot, 'top');
uistack(hBar, 'bottom');

xticks(allCategoricalMeansStd)
xticklabels({'ChRimson Before', 'ChRimson + PdCO', 'ChRimson After'})
xlabel(['Optogenetic activation']);
ylabel(queryRunFeature);
box on;
title([typePattern ' - ' queryRunFeature]);
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');
ylim([-1,3])

plotName = [experimentName '_' queryRunFeature '_' typePattern '_' 'beforeAfter'];
plotFullPath = fullfile(savePlotPath, plotName);

if savePlots == 1
    saveas(gcf, plotFullPath, 'pdf')
    saveas(gcf, plotFullPath, 'fig')
end

epochTrace = 2;
traceFullFieldRedBefore = selectedRunGroupsAnalysisTableRedBefore.trace{epochTrace,1};
traceFullFieldBlueRed = selectedRunGroupsAnalysisTableBlueRed.trace{epochTrace,1};
traceFullFieldRedAfter = selectedRunGroupsAnalysisTableRedAfter.trace{epochTrace,1};

%% comparison traces fullField

typePattern = 'comparison_traces_fullField';
queryRunFeature = 'areaPulse';

selectedCellsIDs = [4];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

selectedRunGroupsRed1 = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao2') & cell2mat(runGroupsParametersTable.nPulsesRed) == 1 & cell2mat(runGroupsParametersTable.pulseWidthRed) == 100 & cell2mat(runGroupsParametersTable.amplitudeRed) == 5 & strcmp(runGroupsParametersTable.functionNameRed,'dmdOutputRect(1, 1, 1, 1, 1, 1)'));
selectedRunGroupsRed2 = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao2') & cell2mat(runGroupsParametersTable.nPulsesRed) == 1 & cell2mat(runGroupsParametersTable.pulseWidthRed) == 100 & cell2mat(runGroupsParametersTable.amplitudeRed) == 5 & strcmp(runGroupsParametersTable.functionNameRed,'dmdOutputRect(1, 1, 1, 1, 2, 1)'));
selectedRunGroupsRed = [selectedRunGroupsRed1; selectedRunGroupsRed2];
selectedRunGroupsAnalysisTableRed = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroupsRed),:);
selectedCyclePositionBefore = find(cell2mat(selectedRunGroupsAnalysisTableRed.cyclePosition)==1);
selectedCyclePositionAfter = find(cell2mat(selectedRunGroupsAnalysisTableRed.cyclePosition)==3);
selectedRunGroupsAnalysisTableRedBefore = selectedRunGroupsAnalysisTableRed(selectedCyclePositionBefore,:);
selectedRunGroupsAnalysisTableRedAfter = selectedRunGroupsAnalysisTableRed(selectedCyclePositionAfter,:);

selectedRunGroupsBlueRed1 = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1, ao2') & cell2mat(runGroupsParametersTable.nPulsesRed) == 1 & cell2mat(runGroupsParametersTable.pulseWidthRed) == 100 & cell2mat(runGroupsParametersTable.delayPulseBlue) == 29500 & cell2mat(runGroupsParametersTable.amplitudeBlue) == 5 & strcmp(runGroupsParametersTable.functionNameBlue,'dmdOutputRect(1, 1, 1, 1, 1, 1)'));
selectedRunGroupsBlueRed2 = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1, ao2') & cell2mat(runGroupsParametersTable.nPulsesRed) == 1 & cell2mat(runGroupsParametersTable.pulseWidthRed) == 100 & cell2mat(runGroupsParametersTable.delayPulseBlue) == 29500 & cell2mat(runGroupsParametersTable.amplitudeBlue) == 5 & strcmp(runGroupsParametersTable.functionNameBlue,'dmdOutputRect(1, 1, 1, 1, 2, 1)'));
selectedRunGroupsBlueRed = [selectedRunGroupsBlueRed1; selectedRunGroupsBlueRed2];
selectedRunGroupsAnalysisTableBlueRed = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroupsBlueRed),:);

traceFullFieldRedBefore = preprocessSignalVC(selectedRunGroupsAnalysisTableRedBefore.trace{3,1});
traceFullFieldBlueRed = preprocessSignalVC(selectedRunGroupsAnalysisTableBlueRed.trace{1,1});
traceFullFieldRedAfter = preprocessSignalVC(selectedRunGroupsAnalysisTableRedAfter.trace{3,1});

figure;
plot(traceFullFieldRedBefore,'LineWidth',2)
hold on;
plot(traceFullFieldBlueRed,'LineWidth',2)
hold on;
plot(traceFullFieldRedAfter,'LineWidth',2)

ylim([-10, 80]);
xLimMin = 3;
xLimMax = 7;
xTicks = linspace(xLimMin,xLimMax, 9);
xlim([xLimMin*1000, xLimMax*1000]);
set(gca, 'XTickLabel', num2str((xTicks'-xLimMin)))
xlabel('Time [s]');
ylabel('Current [pA]');
box on;
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');
legend('ChRimson Before','ChRimson + PdCO','ChRimson After','Location', 'northeast', 'FontSize',12)

plotName = [experimentName '_' typePattern];
plotFullPath = fullfile(savePlotPath, plotName);

if savePlots == 1
    saveas(gcf, plotFullPath, 'pdf')
    saveas(gcf, plotFullPath, 'fig')
end

%% comparison traces spotOnCell

typePattern = 'comparison_traces_spotOnCell';
queryRunFeature = 'areaPulse';

selectedCellsIDs = [4];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

selectedRunGroupsRed = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao2') & cell2mat(runGroupsParametersTable.pulseWidthRed) == 20 & cell2mat(runGroupsParametersTable.amplitudeRed) == 5 & strcmp(runGroupsParametersTable.functionNameRed,'[onCellPattern, outsideCellPattern] = activationPatternForPdCO(3,1,0,200,20)'));
selectedRunGroupsAnalysisTableRed = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroupsRed),:);
selectedCyclePositionBefore = find(cell2mat(selectedRunGroupsAnalysisTableRed.cyclePosition)==1);
selectedCyclePositionAfter = find(cell2mat(selectedRunGroupsAnalysisTableRed.cyclePosition)==3);
selectedRunGroupsAnalysisTableRedBefore = selectedRunGroupsAnalysisTableRed(selectedCyclePositionBefore,:);
selectedRunGroupsAnalysisTableRedAfter = selectedRunGroupsAnalysisTableRed(selectedCyclePositionAfter,:);

selectedRunGroupsBlueRed = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1, ao2') & cell2mat(runGroupsParametersTable.pulseWidthRed) == 20 & cell2mat(runGroupsParametersTable.delayPulseBlue) == 28500 & cell2mat(runGroupsParametersTable.amplitudeBlue) == 5 & strcmp(runGroupsParametersTable.functionNameBlue,'[onCellPattern, outsideCellPattern] = activationPatternForPdCO(3,1,0,200,20)'));
selectedRunGroupsAnalysisTableBlueRed = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroupsBlueRed),:);

epochTrace = 6;
traceSpotOnCellRedBefore = preprocessSignalVC(selectedRunGroupsAnalysisTableRedBefore.trace{epochTrace,1});
traceSpotOnCellBlueRed = preprocessSignalVC(selectedRunGroupsAnalysisTableBlueRed.trace{epochTrace,1});
traceSpotOnCellRedAfter = preprocessSignalVC(selectedRunGroupsAnalysisTableRedAfter.trace{epochTrace,1});

figure;
plot(traceSpotOnCellRedBefore,'LineWidth',2)
hold on;
plot(traceSpotOnCellBlueRed,'LineWidth',2)
hold on;
plot(traceSpotOnCellRedAfter,'LineWidth',2)

ylim([-10, 80]);
xLimMin = 3;
xLimMax = 7;
xTicks = linspace(xLimMin,xLimMax, 9);
xlim([xLimMin*1000, xLimMax*1000]);
set(gca, 'XTickLabel', num2str((xTicks'-xLimMin)))
xlabel('Time [s]');
ylabel('Current [pA]');
box on;
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');
legend('ChRimson Before','ChRimson + PdCO','ChRimson After','Location', 'northeast', 'FontSize',12)

plotName = [experimentName '_' typePattern];
plotFullPath = fullfile(savePlotPath, plotName);

if savePlots == 1
    saveas(gcf, plotFullPath, 'pdf')
    saveas(gcf, plotFullPath, 'fig')
end

%% comparison traces withPdCOonly

typePattern = 'comparison_traces_withPdCOonly';
queryRunFeature = 'areaPulse';

selectedCellsIDs = [4];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

selectedRunGroupsRed1 = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao2') & cell2mat(runGroupsParametersTable.nPulsesRed) == 1 & cell2mat(runGroupsParametersTable.pulseWidthRed) == 100 & cell2mat(runGroupsParametersTable.amplitudeRed) == 5 & strcmp(runGroupsParametersTable.functionNameRed,'dmdOutputRect(1, 1, 1, 1, 1, 1)'));
selectedRunGroupsRed2 = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao2') & cell2mat(runGroupsParametersTable.nPulsesRed) == 1 & cell2mat(runGroupsParametersTable.pulseWidthRed) == 100 & cell2mat(runGroupsParametersTable.amplitudeRed) == 5 & strcmp(runGroupsParametersTable.functionNameRed,'dmdOutputRect(1, 1, 1, 1, 2, 1)'));
selectedRunGroupsRed = [selectedRunGroupsRed1; selectedRunGroupsRed2];
selectedRunGroupsAnalysisTableRed = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroupsRed),:);
notExcludedEpochs = find(cell2mat(selectedRunGroupsAnalysisTableRed.epoch) ~= 37 & cell2mat(selectedRunGroupsAnalysisTableRed.epoch) ~= 8);
selectedRunGroupsAnalysisTableRed = selectedRunGroupsAnalysisTableRed(notExcludedEpochs,:);

selectedRunGroupsBlueRed1 = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1, ao2') & cell2mat(runGroupsParametersTable.nPulsesRed) == 1 & cell2mat(runGroupsParametersTable.pulseWidthRed) == 100 & cell2mat(runGroupsParametersTable.delayPulseBlue) == 29500 & cell2mat(runGroupsParametersTable.amplitudeBlue) == 5 & strcmp(runGroupsParametersTable.functionNameBlue,'dmdOutputRect(1, 1, 1, 1, 1, 1)'));
selectedRunGroupsBlueRed2 = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1, ao2') & cell2mat(runGroupsParametersTable.nPulsesRed) == 1 & cell2mat(runGroupsParametersTable.pulseWidthRed) == 100 & cell2mat(runGroupsParametersTable.delayPulseBlue) == 29500 & cell2mat(runGroupsParametersTable.amplitudeBlue) == 5 & strcmp(runGroupsParametersTable.functionNameBlue,'dmdOutputRect(1, 1, 1, 1, 2, 1)'));
selectedRunGroupsBlueRed = [selectedRunGroupsBlueRed1; selectedRunGroupsBlueRed2];
selectedRunGroupsAnalysisTableBlueRed = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroupsBlueRed),:);
notExcludedEpochs = find(cell2mat(selectedRunGroupsAnalysisTableBlueRed.epoch) ~= 37 & cell2mat(selectedRunGroupsAnalysisTableBlueRed.epoch) ~= 8);
selectedRunGroupsAnalysisTableBlueRed = selectedRunGroupsAnalysisTableBlueRed(notExcludedEpochs,:);

selectedCellsIDs = [3];
selectedCellsRows = find(ismember(cell2mat(mouseAnalysisTable.cellName), selectedCellsIDs));
mouseSelectedCellsAnalysisTable = mouseAnalysisTable(selectedCellsRows,:);

selectedRunGroupsBlue = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.pulseWidthBlue) == 100);
selectedRunGroupsAnalysisTableBlue = mouseSelectedCellsAnalysisTable(ismember(mouseSelectedCellsAnalysisTable.runGroup,selectedRunGroupsBlue),:);
notExcludedEpochs = find(cell2mat(selectedRunGroupsAnalysisTableBlue.epoch) ~= 37 & cell2mat(selectedRunGroupsAnalysisTableBlue.epoch) ~= 8);
selectedRunGroupsAnalysisTableBlue = selectedRunGroupsAnalysisTableBlue(notExcludedEpochs,:);

traceRed = preprocessSignalVC(selectedRunGroupsAnalysisTableRed.trace{3,1});
traceBlue = preprocessSignalVC(selectedRunGroupsAnalysisTableBlue.trace{1,1});
traceBlueRed = preprocessSignalVC(selectedRunGroupsAnalysisTableBlueRed.trace{1,1});

figure;
plot(traceRed,'LineWidth',2)
hold on;
plot(traceBlue,'LineWidth',2)
hold on;
plot(traceBlueRed,'LineWidth',2)

ylim([-10, 80]);
xLimMin = 3;
xLimMax = 7;
xTicks = linspace(xLimMin,xLimMax, 9);
xlim([xLimMin*1000, xLimMax*1000]);
set(gca, 'XTickLabel', num2str((xTicks'-xLimMin)))
xlabel('Time [s]');
ylabel('Current [pA]');
box on;
set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');
legend('ChRimson','PdCO','ChRimson + PdCO','Location', 'northeast', 'FontSize',12)

plotName = [experimentName '_' typePattern];
plotFullPath = fullfile(savePlotPath, plotName);

if savePlots == 1
    saveas(gcf, plotFullPath, 'pdf')
    saveas(gcf, plotFullPath, 'fig')
end