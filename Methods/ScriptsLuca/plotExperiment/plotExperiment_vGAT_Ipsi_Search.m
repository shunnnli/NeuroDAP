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

overallExperimentName = 'vGAT_Ipsi_Search';

mouseNames = {};
mouseNames{1} = 'vGAT_SC_Ephys_3_2';
mouseNames{2} = 'vGAT_IC_Ephys_3';

mouseCells = {};
mouseCells{1} = [1,3,5];
mouseCells{2} = [7,8,10,11];

mouseEpochs = {};
mouseEpochs{1} = {[8],[13],[13]};
mouseEpochs{2} = {[6],[5],[5],[5]};

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
