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

overallExperimentName = 'SNr_Search';

mouseNames = {};
mouseNames{1} = 'DCN_SNr_Ephys_1_1';
mouseNames{2} = 'DCN_SNr_Ephys_3_2';
mouseNames{3} = 'DCN_SNr_Ephys_2_1';
mouseNames{4} = 'DCN_SNr_Ephys_2_2';
mouseNames{5} = 'SNrtoSC_Ephys_1_1';

mouseCells = {};
mouseCells{1} = [7];
mouseCells{2} = [2,3,10,11];
mouseCells{3} = [5];
mouseCells{4} = [1,10];
mouseCells{5} = [3,4];

mouseEpochs = {};
mouseEpochs{1} = {[12]};
mouseEpochs{2} = {[9],[10],[9],[15]}; %check c2, 9
mouseEpochs{3} = {[8]};
mouseEpochs{4} = {[8],[17]};
mouseEpochs{5} = {[5],[5]};

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