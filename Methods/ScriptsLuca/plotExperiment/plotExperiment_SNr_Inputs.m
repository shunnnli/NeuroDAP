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

overallExperimentName = 'SNr_Inputs';

mouseExperimentNames = {};
mouseExperimentNames{1} = 'MLR_SNr_Ephys_1';
mouseExperimentNames{2} = 'DCN_SNr_Ephys_1_1';
mouseExperimentNames{3} = 'DCN_SNr_Ephys_2_1';
mouseExperimentNames{4} = 'DCN_SNr_Ephys_2_2';
mouseExperimentNames{5} = 'DCN_SNr_Ephys_3_1';
mouseExperimentNames{6} = 'DCN_SNr_Ephys_3_2';

channel = 'red';
responseType = 'inhibitory';

queryParametersToFilter = {};
queryParametersToFilter{1} = {'cellType','IRT'};

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