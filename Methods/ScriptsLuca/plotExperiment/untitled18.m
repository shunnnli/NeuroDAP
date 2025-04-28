%% Plot experiment

%% Matlab and data directories

clear all;
clc;

% User environment (automatic)
user = getenv('USER'); 

% Indicate main directory and add it to path
matlabDirectory = ['/Users/' user '/Desktop/EPFL/Thesis/Matlab/'];
mainDirectory = ['/Users/' user '/Desktop/EPFL/Thesis/Matlab/Experiments/']; % including all Matlab toolboxes and necessary paths
addpath(genpath(matlabDirectory));

%% Options

savePlots = 1;

%% Indicate individual experiments to group

overallExperimentName = 'Test_Inputs_Ephys';

mouseExperimentNames = {};
%mouseExperimentNames{1} = 'DCN_ZI_Ephys_2_1';
mouseExperimentNames{1} = 'DCN_ZI_Ephys_2_2';
mouseExperimentNames{2} = 'DCN_ZI_Ephys_3_1';
%mouseExperimentNames{4} = 'DCN_ZI_Ephys_3_2';
mouseExperimentNames{3} = 'DCN_SNr_Ephys_2_1';

channel = 'red';

%% Housekeeping - Define and organize overall experiment directory

experimentName = overallExperimentName;
experimentFolderPath = [matlabDirectory 'GroupedExperiments' filesep overallExperimentName];

if ~exist(experimentFolderPath, 'dir'); mkdir(experimentFolderPath); addpath(genpath(experimentFolderPath)); end
[experimentAnalysisTable, runGroupsParametersTable, experimentCellsConnectionInfo] = organizeOverallExperimentFolder(matlabDirectory,mainDirectory,overallExperimentName,mouseExperimentNames);

%% Plots

scriptPlotOverallExperiment;