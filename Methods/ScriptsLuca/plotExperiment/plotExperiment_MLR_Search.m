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

overallExperimentName = 'MLR_Search';

mouseNames = {};
mouseNames{1} = 'MLR_SNr_Ephys_1';
mouseNames{2} = 'MLR_ZI_Ephys_2_3';

mouseCells = {};
mouseCells{1} = [1];
mouseCells{2} = [2,3];

mouseEpochs = {};
mouseEpochs{1} = {[11]};
mouseEpochs{2} = {[5],[8]};

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