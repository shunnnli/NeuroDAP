% Clearing workspace
clear all;
close all;
clc;

% User environment (automatic)
user = getenv('USER'); 

% Indicate main directory and add it to path
matlabDirectory = ['/Users/' user '/Desktop/EPFL/Thesis/Matlab/'];
mainDirectory = ['/Users/' user '/Desktop/EPFL/Thesis/Matlab/Experiments/']; % including all Matlab toolboxes and necessary paths
addpath(genpath(matlabDirectory));

% Indicate experiment to analyze
experimentName = 'DCN_ZI_Ephys_1_2';

% Set path
experimentDirectory = [mainDirectory experimentName];
mouseAnalysisTablePath = [experimentDirectory filesep 'mouseAnalysis' filesep 'MouseAnalysisTable.mat'];
mouseParametersTablePath = [experimentDirectory filesep 'mouseAnalysis' filesep 'MouseParametersTable.mat'];

% Load tables
mouseAnalysisTable = load(mouseAnalysisTablePath);
mouseAnalysisTable = mouseAnalysisTable.mouseAnalysisTable;
mouseParametersTable = load(mouseParametersTablePath);
runGroupsParametersTable = mouseParametersTable.runGroupsParametersTable;

% Decide saving options
savePlotPath = [experimentDirectory filesep 'mouseAnalysis'];
savePlots = 1;

%% amplitudeBlue vs heightPulsePeak

cat1_queryOptoParameter =  'amplitudeBlue';
queryRunFeature = 'heightPulsePeak';

plotBaseName = [experimentName '_' 'amplitudeBlue' '_' queryRunFeature]; 
queryOptoParameterValuesOverall = [5];

cat1_selectedCellsIDs = [4,7,8,9];

cat1_desiredEpochs = [];

groupByCell = 1;

cat1_selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.amplitudeBlue) == 5 & cell2mat(runGroupsParametersTable.holdingVoltage) == -70 & strcmp(runGroupsParametersTable.functionNameBlue,'dmdOutputRect(1, 1, 1, 1, 1, 1)'));

% No need to modify
cat1_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat1_selectedCellsIDs);
cat1_selectedRunGroupsAnalysisTable = cat1_selectedMembersTable(ismember(cat1_selectedMembersTable.runGroup,cat1_selectedRunGroups),:);
cat1_selectedRunGroupsAnalysisTable = extractEpochSubset(cat1_selectedRunGroupsAnalysisTable, cat1_desiredEpochs);

generateBarPlotVCwith1category(cat1_queryOptoParameter, queryRunFeature, queryOptoParameterValuesOverall, cat1_selectedRunGroupsAnalysisTable, cat1_selectedRunGroups, cat1_selectedCellsIDs, groupByCell, plotBaseName, savePlotPath, savePlots)

%% amplitudeBlue vs areaPulse

cat1_queryOptoParameter =  'amplitudeBlue';
queryRunFeature = 'areaPulse';

plotBaseName = [experimentName '_' 'amplitudeBlue' '_' queryRunFeature]; 
queryOptoParameterValuesOverall = [5];

cat1_selectedCellsIDs = [4,7,8,9];

cat1_desiredEpochs = [];

groupByCell = 1;

cat1_selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.amplitudeBlue) == 5 & cell2mat(runGroupsParametersTable.holdingVoltage) == -70 & strcmp(runGroupsParametersTable.functionNameBlue,'dmdOutputRect(1, 1, 1, 1, 1, 1)'));

% No need to modify
cat1_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat1_selectedCellsIDs);
cat1_selectedRunGroupsAnalysisTable = cat1_selectedMembersTable(ismember(cat1_selectedMembersTable.runGroup,cat1_selectedRunGroups),:);
cat1_selectedRunGroupsAnalysisTable = extractEpochSubset(cat1_selectedRunGroupsAnalysisTable, cat1_desiredEpochs);

generateBarPlotVCwith1category(cat1_queryOptoParameter, queryRunFeature, queryOptoParameterValuesOverall, cat1_selectedRunGroupsAnalysisTable, cat1_selectedRunGroups, cat1_selectedCellsIDs, groupByCell, plotBaseName, savePlotPath, savePlots)

%% amplitude vs areaPulse - connection strength

cat1_queryOptoParameter =  'amplitudeBlue';
cat2_queryOptoParameter =  'amplitudeBlue';
queryRunFeature = 'areaPulse';

plotBaseName = [experimentName '_' 'connectionStrength' '_' queryRunFeature]; 
queryOptoParameterValuesOverall = [5];

cat1_selectedCellsIDs = [4,7,8,9];
cat2_selectedCellsIDs = [1,2];

cat1_desiredEpochs = [];
cat2_desiredEpochs = [];

groupByCell = 1;

cat1_selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.holdingVoltage) == -70 & cell2mat(runGroupsParametersTable.amplitudeBlue) == 5 & strcmp(runGroupsParametersTable.functionNameBlue,'dmdOutputRect(1, 1, 1, 1, 1, 1)'));
cat2_selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.holdingVoltage) == -70 & cell2mat(runGroupsParametersTable.amplitudeBlue) == 5 & strcmp(runGroupsParametersTable.functionNameBlue,'dmdOutputRect(1, 1, 1, 1, 1, 1)'));

% No need to modify
cat1_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat1_selectedCellsIDs);
cat1_selectedRunGroupsAnalysisTable = cat1_selectedMembersTable(ismember(cat1_selectedMembersTable.runGroup,cat1_selectedRunGroups),:);
cat1_selectedRunGroupsAnalysisTable = extractEpochSubset(cat1_selectedRunGroupsAnalysisTable, cat1_desiredEpochs);

cat2_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat2_selectedCellsIDs);
cat2_selectedRunGroupsAnalysisTable = cat2_selectedMembersTable(ismember(cat2_selectedMembersTable.runGroup,cat2_selectedRunGroups),:);
cat2_selectedRunGroupsAnalysisTable = extractEpochSubset(cat2_selectedRunGroupsAnalysisTable, cat2_desiredEpochs);

generateBarPlotVCwith2categoriesMouse(cat1_queryOptoParameter, cat2_queryOptoParameter, queryRunFeature, queryOptoParameterValuesOverall, cat1_selectedRunGroupsAnalysisTable, cat2_selectedRunGroupsAnalysisTable, cat1_selectedRunGroups, cat2_selectedRunGroups, cat1_selectedCellsIDs, cat2_selectedCellsIDs, groupByCell, plotBaseName, savePlotPath, savePlots)

%% amplitude vs heightPulsePeak - connection strength

cat1_queryOptoParameter =  'amplitudeBlue';
cat2_queryOptoParameter =  'amplitudeBlue';
queryRunFeature = 'heightPulsePeak';

plotBaseName = [experimentName '_' 'connectionStrength' '_' queryRunFeature]; 
queryOptoParameterValuesOverall = [5];

cat1_selectedCellsIDs = [4,7,8,9];
cat2_selectedCellsIDs = [1,2];

cat1_desiredEpochs = [];
cat2_desiredEpochs = [];

groupByCell = 1;

cat1_selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.holdingVoltage) == -70 & cell2mat(runGroupsParametersTable.amplitudeBlue) == 5 & strcmp(runGroupsParametersTable.functionNameBlue,'dmdOutputRect(1, 1, 1, 1, 1, 1)'));
cat2_selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.holdingVoltage) == -70 & cell2mat(runGroupsParametersTable.amplitudeBlue) == 5 & strcmp(runGroupsParametersTable.functionNameBlue,'dmdOutputRect(1, 1, 1, 1, 1, 1)'));

% No need to modify
cat1_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat1_selectedCellsIDs);
cat1_selectedRunGroupsAnalysisTable = cat1_selectedMembersTable(ismember(cat1_selectedMembersTable.runGroup,cat1_selectedRunGroups),:);
cat1_selectedRunGroupsAnalysisTable = extractEpochSubset(cat1_selectedRunGroupsAnalysisTable, cat1_desiredEpochs);

cat2_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat2_selectedCellsIDs);
cat2_selectedRunGroupsAnalysisTable = cat2_selectedMembersTable(ismember(cat2_selectedMembersTable.runGroup,cat2_selectedRunGroups),:);
cat2_selectedRunGroupsAnalysisTable = extractEpochSubset(cat2_selectedRunGroupsAnalysisTable, cat2_desiredEpochs);

generateBarPlotVCwith2categoriesMouse(cat1_queryOptoParameter, cat2_queryOptoParameter, queryRunFeature, queryOptoParameterValuesOverall, cat1_selectedRunGroupsAnalysisTable, cat2_selectedRunGroupsAnalysisTable, cat1_selectedRunGroups, cat2_selectedRunGroups, cat1_selectedCellsIDs, cat2_selectedCellsIDs, groupByCell, plotBaseName, savePlotPath, savePlots)

%% amplitudeRed vs heightPulsePeak

cat1_queryOptoParameter =  'amplitudeRed';
queryRunFeature = 'heightPulsePeak';

plotBaseName = [experimentName '_' 'amplitudeRed' '_' queryRunFeature]; 
queryOptoParameterValuesOverall = [5];

cat1_selectedCellsIDs = [6];%1,2,4
cat1_desiredEpochs = [];

groupByCell = 0;

cat1_selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao2') & cell2mat(runGroupsParametersTable.amplitudeRed) == 5 & cell2mat(runGroupsParametersTable.holdingVoltage) == 0 & strcmp(runGroupsParametersTable.functionNameRed,'dmdOutputRect(1, 1, 1, 1, 2, 1)'));

% No need to modify
cat1_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat1_selectedCellsIDs);
cat1_selectedRunGroupsAnalysisTable = cat1_selectedMembersTable(ismember(cat1_selectedMembersTable.runGroup,cat1_selectedRunGroups),:);
cat1_selectedRunGroupsAnalysisTable = extractEpochSubset(cat1_selectedRunGroupsAnalysisTable, cat1_desiredEpochs);

generateBarPlotVCwith1category(cat1_queryOptoParameter, queryRunFeature, queryOptoParameterValuesOverall, cat1_selectedRunGroupsAnalysisTable, cat1_selectedRunGroups, cat1_selectedCellsIDs, groupByCell, plotBaseName, savePlotPath, savePlots)

%% amplitudeRed vs areaPulse

cat1_queryOptoParameter =  'amplitudeRed';
queryRunFeature = 'areaPulse';

plotBaseName = [experimentName '_' 'amplitudeRed' '_' queryRunFeature]; 
queryOptoParameterValuesOverall = [5];

cat1_selectedCellsIDs = [6];
cat1_desiredEpochs = [];

groupByCell = 0;

cat1_selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao2') & cell2mat(runGroupsParametersTable.amplitudeRed) == 5 & cell2mat(runGroupsParametersTable.holdingVoltage) == 0 & strcmp(runGroupsParametersTable.functionNameRed,'dmdOutputRect(1, 1, 1, 1, 2, 1)'));

% No need to modify
cat1_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat1_selectedCellsIDs);
cat1_selectedRunGroupsAnalysisTable = cat1_selectedMembersTable(ismember(cat1_selectedMembersTable.runGroup,cat1_selectedRunGroups),:);
cat1_selectedRunGroupsAnalysisTable = extractEpochSubset(cat1_selectedRunGroupsAnalysisTable, cat1_desiredEpochs);

generateBarPlotVCwith1category(cat1_queryOptoParameter, queryRunFeature, queryOptoParameterValuesOverall, cat1_selectedRunGroupsAnalysisTable, cat1_selectedRunGroups, cat1_selectedCellsIDs, groupByCell, plotBaseName, savePlotPath, savePlots)

%% amplitude vs heightPulsePeak - noDrugsNBQX

cat1_queryOptoParameter =  'amplitudeBlue';
cat2_queryOptoParameter =  'amplitudeBlue';
queryRunFeature = 'heightPulsePeak';

plotBaseName = [experimentName '_' 'noDrugsNBQX' '_' queryRunFeature]; 
queryOptoParameterValuesOverall = [5];

cat1_selectedCellsIDs = [9];
cat2_selectedCellsIDs = [9];

cat1_desiredEpochs = [1:9];
cat2_desiredEpochs = [11:14];

groupByCell = 0;

cat1_selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.holdingVoltage) == -70 & cell2mat(runGroupsParametersTable.amplitudeBlue) == 5 & strcmp(runGroupsParametersTable.functionNameBlue,'dmdOutputRect(1, 1, 1, 1, 1, 1)'));
cat2_selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.holdingVoltage) == -70 & cell2mat(runGroupsParametersTable.amplitudeBlue) == 5 & strcmp(runGroupsParametersTable.functionNameBlue,'dmdOutputRect(1, 1, 1, 1, 1, 1)'));

% No need to modify
cat1_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat1_selectedCellsIDs);
cat1_selectedRunGroupsAnalysisTable = cat1_selectedMembersTable(ismember(cat1_selectedMembersTable.runGroup,cat1_selectedRunGroups),:);
cat1_selectedRunGroupsAnalysisTable = extractEpochSubset(cat1_selectedRunGroupsAnalysisTable, cat1_desiredEpochs);

cat2_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat2_selectedCellsIDs);
cat2_selectedRunGroupsAnalysisTable = cat2_selectedMembersTable(ismember(cat2_selectedMembersTable.runGroup,cat2_selectedRunGroups),:);
cat2_selectedRunGroupsAnalysisTable = extractEpochSubset(cat2_selectedRunGroupsAnalysisTable, cat2_desiredEpochs);

generateBarPlotVCwith2categoriesMouse(cat1_queryOptoParameter, cat2_queryOptoParameter, queryRunFeature, queryOptoParameterValuesOverall, cat1_selectedRunGroupsAnalysisTable, cat2_selectedRunGroupsAnalysisTable, cat1_selectedRunGroups, cat2_selectedRunGroups, cat1_selectedCellsIDs, cat2_selectedCellsIDs, groupByCell, plotBaseName, savePlotPath, savePlots)

%% amplitude vs areaPulse - noDrugsNBQX

cat1_queryOptoParameter =  'amplitudeBlue';
cat2_queryOptoParameter =  'amplitudeBlue';
queryRunFeature = 'areaPulse';

plotBaseName = [experimentName '_' 'noDrugsNBQX' '_' queryRunFeature]; 
queryOptoParameterValuesOverall = [5];

cat1_selectedCellsIDs = [9];
cat2_selectedCellsIDs = [9];

cat1_desiredEpochs = [1:9];
cat2_desiredEpochs = [11:14];

groupByCell = 0;

cat1_selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.holdingVoltage) == -70 & cell2mat(runGroupsParametersTable.amplitudeBlue) == 5 & strcmp(runGroupsParametersTable.functionNameBlue,'dmdOutputRect(1, 1, 1, 1, 1, 1)'));
cat2_selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.holdingVoltage) == -70 & cell2mat(runGroupsParametersTable.amplitudeBlue) == 5 & strcmp(runGroupsParametersTable.functionNameBlue,'dmdOutputRect(1, 1, 1, 1, 1, 1)'));

% No need to modify
cat1_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat1_selectedCellsIDs);
cat1_selectedRunGroupsAnalysisTable = cat1_selectedMembersTable(ismember(cat1_selectedMembersTable.runGroup,cat1_selectedRunGroups),:);
cat1_selectedRunGroupsAnalysisTable = extractEpochSubset(cat1_selectedRunGroupsAnalysisTable, cat1_desiredEpochs);

cat2_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat2_selectedCellsIDs);
cat2_selectedRunGroupsAnalysisTable = cat2_selectedMembersTable(ismember(cat2_selectedMembersTable.runGroup,cat2_selectedRunGroups),:);
cat2_selectedRunGroupsAnalysisTable = extractEpochSubset(cat2_selectedRunGroupsAnalysisTable, cat2_desiredEpochs);

generateBarPlotVCwith2categoriesMouse(cat1_queryOptoParameter, cat2_queryOptoParameter, queryRunFeature, queryOptoParameterValuesOverall, cat1_selectedRunGroupsAnalysisTable, cat2_selectedRunGroupsAnalysisTable, cat1_selectedRunGroups, cat2_selectedRunGroups, cat1_selectedCellsIDs, cat2_selectedCellsIDs, groupByCell, plotBaseName, savePlotPath, savePlots)

%% amplitude vs areaPulse - noDrugsCharbachol

cat1_queryOptoParameter =  'amplitudeBlue';
cat2_queryOptoParameter =  'amplitudeBlue';
queryRunFeature = 'areaPulse';

plotBaseName = [experimentName '_' 'noDrugsCharbachol' '_' queryRunFeature]; 
queryOptoParameterValuesOverall = [5];

cat1_selectedCellsIDs = [9];
cat2_selectedCellsIDs = [9];

cat1_desiredEpochs = [1:4];
cat2_desiredEpochs = [5:10];

groupByCell = 0;

cat1_selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.holdingVoltage) == -70 & cell2mat(runGroupsParametersTable.amplitudeBlue) == 5 & strcmp(runGroupsParametersTable.functionNameBlue,'dmdOutputRect(1, 1, 1, 1, 1, 1)'));
cat2_selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.holdingVoltage) == -70 & cell2mat(runGroupsParametersTable.amplitudeBlue) == 5 & strcmp(runGroupsParametersTable.functionNameBlue,'dmdOutputRect(1, 1, 1, 1, 1, 1)'));

% No need to modify
cat1_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat1_selectedCellsIDs);
cat1_selectedRunGroupsAnalysisTable = cat1_selectedMembersTable(ismember(cat1_selectedMembersTable.runGroup,cat1_selectedRunGroups),:);
cat1_selectedRunGroupsAnalysisTable = extractEpochSubset(cat1_selectedRunGroupsAnalysisTable, cat1_desiredEpochs);

cat2_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat2_selectedCellsIDs);
cat2_selectedRunGroupsAnalysisTable = cat2_selectedMembersTable(ismember(cat2_selectedMembersTable.runGroup,cat2_selectedRunGroups),:);
cat2_selectedRunGroupsAnalysisTable = extractEpochSubset(cat2_selectedRunGroupsAnalysisTable, cat2_desiredEpochs);

generateBarPlotVCwith2categoriesMouse(cat1_queryOptoParameter, cat2_queryOptoParameter, queryRunFeature, queryOptoParameterValuesOverall, cat1_selectedRunGroupsAnalysisTable, cat2_selectedRunGroupsAnalysisTable, cat1_selectedRunGroups, cat2_selectedRunGroups, cat1_selectedCellsIDs, cat2_selectedCellsIDs, groupByCell, plotBaseName, savePlotPath, savePlots)

%% amplitude vs heightPulsePeak - noDrugsCharbachol

cat1_queryOptoParameter =  'amplitudeBlue';
cat2_queryOptoParameter =  'amplitudeBlue';
queryRunFeature = 'heightPulsePeak';

plotBaseName = [experimentName '_' 'noDrugsCharbachol' '_' queryRunFeature]; 
queryOptoParameterValuesOverall = [5];

cat1_selectedCellsIDs = [9];
cat2_selectedCellsIDs = [9];

cat1_desiredEpochs = [1:4];
cat2_desiredEpochs = [5:10];

groupByCell = 0;

cat1_selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.holdingVoltage) == -70 & cell2mat(runGroupsParametersTable.amplitudeBlue) == 5 & strcmp(runGroupsParametersTable.functionNameBlue,'dmdOutputRect(1, 1, 1, 1, 1, 1)'));
cat2_selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.holdingVoltage) == -70 & cell2mat(runGroupsParametersTable.amplitudeBlue) == 5 & strcmp(runGroupsParametersTable.functionNameBlue,'dmdOutputRect(1, 1, 1, 1, 1, 1)'));

% No need to modify
cat1_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat1_selectedCellsIDs);
cat1_selectedRunGroupsAnalysisTable = cat1_selectedMembersTable(ismember(cat1_selectedMembersTable.runGroup,cat1_selectedRunGroups),:);
cat1_selectedRunGroupsAnalysisTable = extractEpochSubset(cat1_selectedRunGroupsAnalysisTable, cat1_desiredEpochs);

cat2_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat2_selectedCellsIDs);
cat2_selectedRunGroupsAnalysisTable = cat2_selectedMembersTable(ismember(cat2_selectedMembersTable.runGroup,cat2_selectedRunGroups),:);
cat2_selectedRunGroupsAnalysisTable = extractEpochSubset(cat2_selectedRunGroupsAnalysisTable, cat2_desiredEpochs);

generateBarPlotVCwith2categoriesMouse(cat1_queryOptoParameter, cat2_queryOptoParameter, queryRunFeature, queryOptoParameterValuesOverall, cat1_selectedRunGroupsAnalysisTable, cat2_selectedRunGroupsAnalysisTable, cat1_selectedRunGroups, cat2_selectedRunGroups, cat1_selectedCellsIDs, cat2_selectedCellsIDs, groupByCell, plotBaseName, savePlotPath, savePlots)

%% amplitude vs timePulsePeak - noDrugsCharbachol

cat1_queryOptoParameter =  'amplitudeBlue';
cat2_queryOptoParameter =  'amplitudeBlue';
queryRunFeature = 'timePulsePeak';

plotBaseName = [experimentName '_' 'noDrugsCharbachol' '_' queryRunFeature]; 
queryOptoParameterValuesOverall = [5];

cat1_selectedCellsIDs = [9];
cat2_selectedCellsIDs = [9];

cat1_desiredEpochs = [1:4];
cat2_desiredEpochs = [5:10];

groupByCell = 0;

cat1_selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.holdingVoltage) == -70 & cell2mat(runGroupsParametersTable.amplitudeBlue) == 5 & strcmp(runGroupsParametersTable.functionNameBlue,'dmdOutputRect(1, 1, 1, 1, 1, 1)'));
cat2_selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.holdingVoltage) == -70 & cell2mat(runGroupsParametersTable.amplitudeBlue) == 5 & strcmp(runGroupsParametersTable.functionNameBlue,'dmdOutputRect(1, 1, 1, 1, 1, 1)'));

% No need to modify
cat1_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat1_selectedCellsIDs);
cat1_selectedRunGroupsAnalysisTable = cat1_selectedMembersTable(ismember(cat1_selectedMembersTable.runGroup,cat1_selectedRunGroups),:);
cat1_selectedRunGroupsAnalysisTable = extractEpochSubset(cat1_selectedRunGroupsAnalysisTable, cat1_desiredEpochs);

cat2_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat2_selectedCellsIDs);
cat2_selectedRunGroupsAnalysisTable = cat2_selectedMembersTable(ismember(cat2_selectedMembersTable.runGroup,cat2_selectedRunGroups),:);
cat2_selectedRunGroupsAnalysisTable = extractEpochSubset(cat2_selectedRunGroupsAnalysisTable, cat2_desiredEpochs);

generateBarPlotVCwith2categoriesMouse(cat1_queryOptoParameter, cat2_queryOptoParameter, queryRunFeature, queryOptoParameterValuesOverall, cat1_selectedRunGroupsAnalysisTable, cat2_selectedRunGroupsAnalysisTable, cat1_selectedRunGroups, cat2_selectedRunGroups, cat1_selectedCellsIDs, cat2_selectedCellsIDs, groupByCell, plotBaseName, savePlotPath, savePlots)
