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

overallExperimentName = 'vGAT_Contra_notvGAT_Search';

mouseNames = {};
mouseNames{1} = 'vGAT_IC_Ephys_1';
mouseNames{2} = 'vGAT_IC_Ephys_30_1';
mouseNames{3} = 'vGAT_IC_Ephys_30_3';
mouseNames{4} = 'vGAT_IC_Ephys_40_1';
mouseNames{5} = 'vGAT_IC_Ephys_50_1';
mouseNames{6} = 'vGAT_IC_Ephys_50_2';
mouseNames{7} = 'vGAT_IC_Ephys_50_3';

mouseCells = {};
mouseCells{1} = [5];
mouseCells{2} = [1,2,5];
mouseCells{3} = [1,3,4,5];
mouseCells{4} = [1,3,4];
mouseCells{5} = [5];
mouseCells{6} = [1,2,4];
mouseCells{7} = [2,6];

mouseEpochs = {};
mouseEpochs{1} = {[10]};
mouseEpochs{2} = {[7],[5],[5]};
mouseEpochs{3} = {[5],[10],[5],[8]};
mouseEpochs{4} = {[7],[5],[5]};
mouseEpochs{5} = {[6]};
mouseEpochs{6} = {[5],[6],[5]};
mouseEpochs{7} = {[5],[5]};

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