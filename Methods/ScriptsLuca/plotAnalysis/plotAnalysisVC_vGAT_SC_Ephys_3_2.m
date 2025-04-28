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
experimentName = 'vGAT_SC_Ephys_3_2';

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

%% amplitudeBlue vs areaPulse - different power 

cat1_queryOptoParameter =  'amplitudeBlue';
queryRunFeature = 'areaPulse';

plotBaseName = [experimentName '_' 'amplitudeBlue' '_' 'diffPower' '_' queryRunFeature]; 
queryOptoParameterValuesOverall = [0.1,0.3,0.5,1,3,5];

cat1_selectedCellsIDs = [3];

cat1_desiredEpochs = [];

groupByCell = 0;

cat1_selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & strcmp(runGroupsParametersTable.functionNameBlue,'dmdOutputRect(1, 1, 1, 1, 1, 1)'));

% No need to modify
cat1_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat1_selectedCellsIDs);
cat1_selectedRunGroupsAnalysisTable = cat1_selectedMembersTable(ismember(cat1_selectedMembersTable.runGroup,cat1_selectedRunGroups),:);
cat1_selectedRunGroupsAnalysisTable = extractEpochSubset(cat1_selectedRunGroupsAnalysisTable, cat1_desiredEpochs);

generateBarPlotVCwith1category(cat1_queryOptoParameter, queryRunFeature, queryOptoParameterValuesOverall, cat1_selectedRunGroupsAnalysisTable, cat1_selectedRunGroups, cat1_selectedCellsIDs, groupByCell, plotBaseName, savePlotPath, savePlots)

%% amplitudeBlue vs heightPulsePeak - different power 

cat1_queryOptoParameter =  'amplitudeBlue';
queryRunFeature = 'heightPulsePeak';

plotBaseName = [experimentName '_' 'amplitudeBlue' '_' 'diffPower' '_' queryRunFeature]; 
queryOptoParameterValuesOverall = [0.1,0.3,0.5,1,3,5];

cat1_selectedCellsIDs = [3];

cat1_desiredEpochs = [];

groupByCell = 0;

cat1_selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & strcmp(runGroupsParametersTable.functionNameBlue,'dmdOutputRect(1, 1, 1, 1, 1, 1)'));

% No need to modify
cat1_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat1_selectedCellsIDs);
cat1_selectedRunGroupsAnalysisTable = cat1_selectedMembersTable(ismember(cat1_selectedMembersTable.runGroup,cat1_selectedRunGroups),:);
cat1_selectedRunGroupsAnalysisTable = extractEpochSubset(cat1_selectedRunGroupsAnalysisTable, cat1_desiredEpochs);

generateBarPlotVCwith1category(cat1_queryOptoParameter, queryRunFeature, queryOptoParameterValuesOverall, cat1_selectedRunGroupsAnalysisTable, cat1_selectedRunGroups, cat1_selectedCellsIDs, groupByCell, plotBaseName, savePlotPath, savePlots)

%% amplitudeBlue vs heightPulsePeak

cat1_queryOptoParameter =  'amplitudeBlue';
queryRunFeature = 'heightPulsePeak';

plotBaseName = [experimentName '_' 'amplitudeBlue' '_' queryRunFeature]; 
queryOptoParameterValuesOverall = [5];

cat1_selectedCellsIDs = [5];

cat1_desiredEpochs = [];

groupByCell = 0;

cat1_selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.amplitudeBlue) == 5 & cell2mat(runGroupsParametersTable.holdingVoltage) == 0 & strcmp(runGroupsParametersTable.functionNameBlue,'dmdOutputRect(1, 1, 1, 1, 1, 1)'));

% No need to modify
cat1_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat1_selectedCellsIDs);
cat1_selectedRunGroupsAnalysisTable = cat1_selectedMembersTable(ismember(cat1_selectedMembersTable.runGroup,cat1_selectedRunGroups),:);
cat1_selectedRunGroupsAnalysisTable = extractEpochSubset(cat1_selectedRunGroupsAnalysisTable, cat1_desiredEpochs);

generateBarPlotVCwith1category(cat1_queryOptoParameter, queryRunFeature, queryOptoParameterValuesOverall, cat1_selectedRunGroupsAnalysisTable, cat1_selectedRunGroups, cat1_selectedCellsIDs, groupByCell, plotBaseName, savePlotPath, savePlots)

%% amplitude vs areaPulse - connection strength

cat1_queryOptoParameter =  'amplitudeBlue';
cat2_queryOptoParameter =  'amplitudeBlue';
queryRunFeature = 'areaPulse';

plotBaseName = [experimentName '_' 'connectionStrength_ipsiContra' '_' queryRunFeature]; 
queryOptoParameterValuesOverall = [5];

cat1_selectedCellsIDs = [3,4,5];
cat2_selectedCellsIDs = [1,2];

cat1_desiredEpochs = [];
cat2_desiredEpochs = [];

groupByCell = 1;

cat1_selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.holdingVoltage) == 0 & cell2mat(runGroupsParametersTable.amplitudeBlue) == 5 & strcmp(runGroupsParametersTable.functionNameBlue,'dmdOutputRect(1, 1, 1, 1, 1, 1)'));
cat2_selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.holdingVoltage) == 0 & cell2mat(runGroupsParametersTable.amplitudeBlue) == 5 & strcmp(runGroupsParametersTable.functionNameBlue,'dmdOutputRect(1, 1, 1, 1, 1, 1)'));

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

plotBaseName = [experimentName '_' 'connectionStrength_ipsiContra' '_' queryRunFeature]; 
queryOptoParameterValuesOverall = [5];

cat1_selectedCellsIDs = [3,4,5];
cat2_selectedCellsIDs = [1,2];

cat1_desiredEpochs = [];
cat2_desiredEpochs = [];

groupByCell = 1;

cat1_selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.holdingVoltage) == 0 & cell2mat(runGroupsParametersTable.amplitudeBlue) == 5 & strcmp(runGroupsParametersTable.functionNameBlue,'dmdOutputRect(1, 1, 1, 1, 1, 1)'));
cat2_selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.holdingVoltage) == 0 & cell2mat(runGroupsParametersTable.amplitudeBlue) == 5 & strcmp(runGroupsParametersTable.functionNameBlue,'dmdOutputRect(1, 1, 1, 1, 1, 1)'));

% No need to modify
cat1_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat1_selectedCellsIDs);
cat1_selectedRunGroupsAnalysisTable = cat1_selectedMembersTable(ismember(cat1_selectedMembersTable.runGroup,cat1_selectedRunGroups),:);
cat1_selectedRunGroupsAnalysisTable = extractEpochSubset(cat1_selectedRunGroupsAnalysisTable, cat1_desiredEpochs);

cat2_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat2_selectedCellsIDs);
cat2_selectedRunGroupsAnalysisTable = cat2_selectedMembersTable(ismember(cat2_selectedMembersTable.runGroup,cat2_selectedRunGroups),:);
cat2_selectedRunGroupsAnalysisTable = extractEpochSubset(cat2_selectedRunGroupsAnalysisTable, cat2_desiredEpochs);

generateBarPlotVCwith2categoriesMouse(cat1_queryOptoParameter, cat2_queryOptoParameter, queryRunFeature, queryOptoParameterValuesOverall, cat1_selectedRunGroupsAnalysisTable, cat2_selectedRunGroupsAnalysisTable, cat1_selectedRunGroups, cat2_selectedRunGroups, cat1_selectedCellsIDs, cat2_selectedCellsIDs, groupByCell, plotBaseName, savePlotPath, savePlots)

%% amplitude vs timePulsePeak - connection strength

cat1_queryOptoParameter =  'amplitudeBlue';
cat2_queryOptoParameter =  'amplitudeBlue';
queryRunFeature = 'timePulsePeak';

plotBaseName = [experimentName '_' 'connectionStrength_ipsiContra' '_' queryRunFeature]; 
queryOptoParameterValuesOverall = [5];

cat1_selectedCellsIDs = [3,4,5];
cat2_selectedCellsIDs = [1,2];

cat1_desiredEpochs = [];
cat2_desiredEpochs = [];

groupByCell = 1;

cat1_selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.holdingVoltage) == 0 & cell2mat(runGroupsParametersTable.amplitudeBlue) == 5 & strcmp(runGroupsParametersTable.functionNameBlue,'dmdOutputRect(1, 1, 1, 1, 1, 1)'));
cat2_selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTable.holdingVoltage) == 0 & cell2mat(runGroupsParametersTable.amplitudeBlue) == 5 & strcmp(runGroupsParametersTable.functionNameBlue,'dmdOutputRect(1, 1, 1, 1, 1, 1)'));

% No need to modify
cat1_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat1_selectedCellsIDs);
cat1_selectedRunGroupsAnalysisTable = cat1_selectedMembersTable(ismember(cat1_selectedMembersTable.runGroup,cat1_selectedRunGroups),:);
cat1_selectedRunGroupsAnalysisTable = extractEpochSubset(cat1_selectedRunGroupsAnalysisTable, cat1_desiredEpochs);

cat2_selectedMembersTable = extractSelectedMembersFromTable(mouseAnalysisTable, 'cellName', cat2_selectedCellsIDs);
cat2_selectedRunGroupsAnalysisTable = cat2_selectedMembersTable(ismember(cat2_selectedMembersTable.runGroup,cat2_selectedRunGroups),:);
cat2_selectedRunGroupsAnalysisTable = extractEpochSubset(cat2_selectedRunGroupsAnalysisTable, cat2_desiredEpochs);

generateBarPlotVCwith2categoriesMouse(cat1_queryOptoParameter, cat2_queryOptoParameter, queryRunFeature, queryOptoParameterValuesOverall, cat1_selectedRunGroupsAnalysisTable, cat2_selectedRunGroupsAnalysisTable, cat1_selectedRunGroups, cat2_selectedRunGroups, cat1_selectedCellsIDs, cat2_selectedCellsIDs, groupByCell, plotBaseName, savePlotPath, savePlots)
