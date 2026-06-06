%% Plot experiment

%% Matlab and data directories

clear all;
clc;

% User environment (automatic)
% user = getenv('USER'); 

% Indicate main directory and add it to path
matlabDirectory = ['N:\research.files.med.harvard.edu\neurobio\MICROSCOPE\Shun\Analysis\NeuroDAP\'];
dataDirectory = ['N:\MICROSCOPE\Paolo\InVitro_Ephys\'];
mainDirectory = ['C:\Users\paolo\Desktop\EphysAnalysisMaps\']; % including all Matlab toolboxes and necessary paths
addpath(genpath(matlabDirectory));

%% Indicate individual experiments to group

overallExperimentName = 'ZI_Search';

mouseNames = {};
mouseNames{1} = 'DCN_ZI_Ephys_4_2';
mouseNames{2} = 'DCN_ZI_Ephys_3_2';
mouseNames{3} = 'DCN_ZI_Ephys_3_1';
mouseNames{4} = 'DCN_ZI_Ephys_2_2';
mouseNames{5} = 'DCN_ZI_Ephys_5_3';
mouseNames{6} = 'CTX_ZI_Ephys_2';

mouseCells = {};
mouseCells{1} = [11];
mouseCells{2} = [2];
mouseCells{3} = [1,3];
mouseCells{4} = [6];
mouseCells{5} = [3];
mouseCells{6} = [2,3,5,6];

mouseEpochs = {};
mouseEpochs{1} = {[9]};
mouseEpochs{2} = {[16]};
mouseEpochs{3} = {[8],[8]};
mouseEpochs{4} = {[12]};
mouseEpochs{5} = {[9]};
mouseEpochs{6} = {[8],[9],[5],[5]};

selectedSearches = {};
selectedSearches.mouseNames = mouseNames;
selectedSearches.mouseCells = mouseCells;
selectedSearches.mouseEpochs = mouseEpochs;

responseType = 'inhibitory';
recreateExperimentTable = 1;
savePlots = 1;

%% Collect searches for the selected mice

overallExperimentDirectory = [mainDirectory 'GroupedExperiments' filesep overallExperimentName];
if ~exist(overallExperimentDirectory, 'dir'); mkdir(overallExperimentDirectory); addpath(genpath(overallExperimentDirectory)); end

experimentSearchFeaturesTable = organizeOverallExperimentSearchesFolder(mainDirectory,dataDirectory,overallExperimentDirectory,mouseNames,recreateExperimentTable);
selectedSearchesTable = extractSelectedSearches(experimentSearchFeaturesTable,selectedSearches);

%% Plots

scriptPlotSearchSummaryExperiment;