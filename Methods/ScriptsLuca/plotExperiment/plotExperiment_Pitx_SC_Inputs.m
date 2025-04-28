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

%% Indicate individual experiments to group

overallExperimentName = 'Pitx_SC_Inputs';

mouseExperimentNames = {};
mouseExperimentNames{1} = 'Pitx_SC_Ephys_1';
mouseExperimentNames{2} = 'Pitx_SC_Ephys_2';
mouseExperimentNames{3} = 'Pitx_H2BGFP_vGAT_IC_Ephys_1_1';
mouseExperimentNames{4} = 'Pitx_H2BGFP_vGAT_IC_Ephys_1_2';
mouseExperimentNames{5} = 'Pitx_IC_Ephys_20_1';
mouseExperimentNames{6} = 'Pitx_IC_Ephys_10_1';

channel = 'blue|red';
responseType = 'excitatory';

queryParametersToFilter = {};
queryParametersToFilter{1} = {'cellType','NOTPITX'};
queryParametersToFilter{2} = {'cellType','PITX'};
queryParametersToFilter{3} = {'cellLocation','LEFT'};
queryParametersToFilter{4} = {'cellLocation','RIGHT'};
queryParametersToFilter{5} = {'cellLocation','LEFT'; 'cellType','NOTPITX'};
queryParametersToFilter{6} = {'cellLocation','LEFT'; 'cellType','PITX'};
queryParametersToFilter{7} = {'cellLocation','RIGHT'; 'cellType','NOTPITX'};
queryParametersToFilter{8} = {'cellLocation','RIGHT'; 'cellType','PITX'};

recreateExperimentTable = 1;
savePlots = 1;

%% Housekeeping - Define and organize overall experiment directory

experimentName = overallExperimentName;
experimentFolderPath = [matlabDirectory 'GroupedExperiments' filesep overallExperimentName];

if ~exist(experimentFolderPath, 'dir'); mkdir(experimentFolderPath); addpath(genpath(experimentFolderPath)); end
[experimentAnalysisTable, runGroupsParametersTable, experimentCellsConnectionInfo] = organizeOverallExperimentFolder(matlabDirectory,mainDirectory,overallExperimentName,mouseExperimentNames,recreateExperimentTable);

%% Plots

for iFilter = 1:numel(queryParametersToFilter)+1
    
    queryParameters = [];
    postTextSuptitle = {};
    
    if iFilter ~= 1
        
        for iSubfilter = 1:size(queryParametersToFilter{iFilter-1},1)
        
            queryParameters.(queryParametersToFilter{iFilter-1}{iSubfilter,1}) = queryParametersToFilter{iFilter-1}{iSubfilter,2};
            postTextSuptitle{end+1} = {'_',queryParametersToFilter{iFilter-1}{iSubfilter,1},'_',queryParametersToFilter{iFilter-1}{iSubfilter,2}};
        
        end
        
    end
    
    stringArray = cellfun(@(x) strjoin(x,''), postTextSuptitle, 'UniformOutput', false);
    postTextSuptitle = strjoin(stringArray,'');

    scriptPlotOverallExperiment;

end