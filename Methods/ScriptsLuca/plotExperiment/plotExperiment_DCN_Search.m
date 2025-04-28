%% Plot experiment

%% Matlab and data directories

clear all;
clc;

% User environment (automatic)
% user = getenv('USER'); 

% Indicate main directory and add it to path
matlabDirectory = ['P:\research.files.med.harvard.edu\neurobio\MICROSCOPE\Shun\Analysis\NeuroDAP\'];
dataDirectory = ['P:\MICROSCOPE\Paolo\InVitro_Ephys\'];
mainDirectory = ['C:\Users\paolo\Desktop\EphysAnalysisMaps\']; % including all Matlab toolboxes and necessary paths
addpath(genpath(matlabDirectory));

%% Indicate individual experiments to group

overallExperimentName = 'DCN_Search';

mouseNames = {};
mouseNames{1} = 'DCN_ZI_Ephys_3_2';
mouseNames{2} = 'DCN_ZI_Ephys_4_2';
mouseNames{3} = 'DCN_SNr_Ephys_3_2';
mouseNames{4} = 'DCN_ZI_Ephys_3_1';
mouseNames{5} = 'DCN_SNr_Ephys_2_1';
mouseNames{6} = 'DCN_ZI_Ephys_2_2';
mouseNames{7} = 'DCN_SNr_Ephys_2_2';
mouseNames{8} = 'DCN_ZI_Ephys_4_1';
mouseNames{9} = 'DCN_SNr_Ephys_3_1';
mouseNames{10} = 'DCN_SNr_Ephys_1_1';
mouseNames{11} = 'DCN_ZI_Ephys_2_1';
mouseNames{12} = 'DCN_ZI_Ephys_5_1';

mouseCells = {};
mouseCells{1} = [2,4];
mouseCells{2} = [5,10];
mouseCells{3} = [7,11];
mouseCells{4} = [2,5];
mouseCells{5} = [1,9];
mouseCells{6} = [1,6];
mouseCells{7} = [3,5,10,11];
mouseCells{8} = [4,8,9];
mouseCells{9} = [7];
mouseCells{10} = [1,5,7];
mouseCells{11} = [6];
mouseCells{12} = [4,5,7,8,9];

mouseEpochs = {};
mouseEpochs{1} = {[5],[5]};
mouseEpochs{2} = {[5],[5]};
mouseEpochs{3} = {[10],[8]};
mouseEpochs{4} = {[15],[5]};
mouseEpochs{5} = {[10],[5]};
mouseEpochs{6} = {[5],[5]};
mouseEpochs{7} = {[5],[9],[5],[9]};
mouseEpochs{8} = {[5],[5],[10]};
mouseEpochs{9} = {[11]};
mouseEpochs{10} = {[11],[11],[10]};
mouseEpochs{11} = {[5]};
mouseEpochs{12} = {[5],[6],[8],[5],[5]};

selectedSearches = {};
selectedSearches.mouseNames = mouseNames;
selectedSearches.mouseCells = mouseCells;
selectedSearches.mouseEpochs = mouseEpochs;

responseType = 'excitatory';
recreateExperimentTable = 1;
savePlots = 1;

%% Collect searches for the selected mice

overallExperimentDirectory = [mainDirectory 'GroupedExperiments' filesep overallExperimentName];
if ~exist(overallExperimentDirectory, 'dir'); mkdir(overallExperimentDirectory); addpath(genpath(overallExperimentDirectory)); end

experimentSearchFeaturesTable = organizeOverallExperimentSearchesFolder(mainDirectory,dataDirectory,overallExperimentDirectory,mouseNames,recreateExperimentTable);
selectedSearchesTable = extractSelectedSearches(experimentSearchFeaturesTable,selectedSearches);

%% Plots

scriptPlotSearchSummaryExperiment;