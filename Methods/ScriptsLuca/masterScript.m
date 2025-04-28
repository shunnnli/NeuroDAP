%% Master Script

%% Matlab and data directories

clear all;
clc;

% User environment (automatic)
user = getenv('USER'); 

% Indicate main directory and add it to path
matlabDirectory = ['/Users/' user '/Desktop/EPFL/Thesis/Matlab/'];
mainDirectory = ['/Users/' user '/Desktop/EPFL/Thesis/Matlab/Experiments/']; % including all Matlab toolboxes and necessary paths
addpath(genpath(matlabDirectory));

%% Enter experiment and cells to analyze

% Experiment name (e.g., 'Pitx_H2BGFP_vGAT_IC_Ephys_1_1')
experimentName = 'Pitx_H2BGFP_vGAT_Ephys_2_1';

% Response expected in each channel (e.g.,'excitatory','inhibitory','none').
responseTypeBlue = 'excitatory';
responseTypeRed = 'inhibitory';

% Selected cells (e.g., [1,2,3,4,5,6,7,8,9]). If empty, all the existing cells are processed.
selectedCellsNumberID = [];

% Experiment type (e.g., 'VC' or 'CC')
experimentType = 'VC';

% Options
overwriteProcessing = 1;
overwriteAnalysis = 1;
overwriteConnectionInfo = 0;
savePlots = 1;

%% Housekeeping - Set paths and selected cell folders

% Set experiment path
experimentDirectory = [mainDirectory experimentName];

if isempty(selectedCellsNumberID); selectedCellsNumberID = findCellsToProcess(experimentDirectory); end

folderContent = dir([experimentDirectory]);
folderContent = folderContent(startsWith({folderContent.name}, 'cell'));
nCells = sum([folderContent.isdir]);
startDir = 1;

warning('off','MATLAB:unknownObjectNowStruct')
warning('off','MATLAB:listener:callbackError')
selectedCells = arrayfun(@(x) sprintf('cell%d', x), selectedCellsNumberID, 'UniformOutput', false);
selectedCells = string(selectedCells);

%% Call scripts

if strcmp(experimentType,'VC')
    
    % %% General scripts %%
    % scriptOrganizeAcquisitionsVCwithRandomSearch;
    % scriptCollectFeaturesVC;
    
    % %% Standard plotting scripts %%
    % scriptPlotAnalysisCell;
    % scriptPlotAnalysisSummary;
    % scriptPlotHeatmapVerification;
    % scriptPlotPlasticityProtocol;

    % %% Other scripts %%
    % scriptCollectFeaturesRandomSearchNew;
    % scriptPlotAnalysisCellInhibitoryOpsinTesting;
    
elseif strcmp(experimentType,'CC')
    
    % %% General scripts %%
    % scriptOrganizeAcquisitionsCC;
    % scriptCollectFeaturesCC;

end
