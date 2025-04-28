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

overallExperimentName = 'vGAT_SC_Inputs';

mouseExperimentNames = {};
mouseExperimentNames{1} = 'vGAT_IC_Ephys_1';
mouseExperimentNames{2} = 'vGAT_IC_Ephys_2';
mouseExperimentNames{3} = 'vGAT_IC_Ephys_3';
mouseExperimentNames{4} = 'vGAT_SC_Ephys_3';
mouseExperimentNames{5} = 'vGAT_SC_Ephys_2_1';
mouseExperimentNames{6} = 'vGAT_SC_Ephys_3_2';
mouseExperimentNames{7} = 'Pitx_H2BGFP_vGAT_IC_Ephys_1_1';
mouseExperimentNames{8} = 'Pitx_H2BGFP_vGAT_IC_Ephys_1_2';

channel = 'blue|red';
responseType = 'inhibitory';

queryParametersToFilter = {};
queryParametersToFilter{1} = {'cellType','NONE'};
queryParametersToFilter{2} = {'cellType','NOTVGAT'};
queryParametersToFilter{3} = {'cellType','PITX'};
queryParametersToFilter{4} = {'cellType','NOTPITX'};
queryParametersToFilter{5} = {'cellLocation','LEFT'};
queryParametersToFilter{6} = {'cellLocation','RIGHT'};
queryParametersToFilter{7} = {'cellLocation','LEFT'; 'cellType','NOTPITX'};
queryParametersToFilter{8} = {'cellLocation','LEFT'; 'cellType','PITX'};
queryParametersToFilter{9} = {'cellLocation','LEFT'; 'cellType','NOTVGAT'};
queryParametersToFilter{10} = {'cellLocation','LEFT'; 'cellType','VGAT'};
queryParametersToFilter{11} = {'cellLocation','RIGHT'; 'cellType','NOTPITX'};
queryParametersToFilter{12} = {'cellLocation','RIGHT'; 'cellType','PITX'};
queryParametersToFilter{13} = {'cellLocation','RIGHT'; 'cellType','NOTVGAT'};
queryParametersToFilter{14} = {'cellLocation','RIGHT'; 'cellType','VGAT'};

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