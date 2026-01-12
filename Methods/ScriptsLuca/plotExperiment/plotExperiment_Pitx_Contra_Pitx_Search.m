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

overallExperimentName = 'Pitx_Contra_Pitx_Search';

mouseNames = {};
mouseNames{1} = 'Pitx_H2BGFP_vGAT_IC_Ephys_1_2';
mouseNames{2} = 'Pitx_IC_Ephys_50_3';
mouseNames{3} = 'PitxH2BGFP_vGAT__Ephys_3_4';
mouseNames{4} = 'Pitx_IC_Ephys_60_1';

mouseCells = {};
mouseCells{1} = [2];
mouseCells{2} = [1,6];
mouseCells{3} = [2];
mouseCells{4} = [1];

mouseEpochs = {};
mouseEpochs{1} = {[9]};
mouseEpochs{2} = {[10],[8]};
mouseEpochs{3} = {[8]};
mouseEpochs{4} = {[5]};

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