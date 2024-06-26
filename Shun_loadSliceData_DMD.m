% Shun_loadSliceData_DMD

% Load and analyze DMD slice data (single session)

%% Load sessions
clear; close all;
addpath(genpath(osPathSwitch('/Volumes/MICROSCOPE/Shun/Analysis/NeuroDAP/Methods')));

% Select sessions for analysis
% parentPath = osPathSwitch('/Volumes/MICROSCOPE/wengang/Exp_withShun/');
parentPath = osPathSwitch('/Volumes/MICROSCOPE/Shun/Project valence/Patch/');
expPath = uipickfiles('FilterSpec',parentPath,'Prompt','Select experiment folders');
rawDataPath = 'default'; % strcat(parentPath,filesep,'20231221_ally');

% Set comman params
[sessionParams,canceled] = inputSessionParams_singleSlice(expPath,...
                                paradigm=2,redStim=true,...
                                reload=false,calculateQC=false,...
                                timeRange='[-20,100]',nArtifactSamples='10');
taskOptions = {'random','reward pairing','punish pairing'};
task = taskOptions{sessionParams.Paradigm};
timeRange = eval(sessionParams.timeRange);
nArtifactSamples = str2double(sessionParams.nArtifactSamples);

if ~isscalar(expPath); error('Multiple sessions were selected, for multi-session analysis see analyzeSlice pipeline!'); end
expPath = expPath{1};

%%

cellList = dir(fullfile(expPath,'cell*'));
nCells = length(cellList);
startDir = 1;

loadCells = struct;
selectedCells = ["cell1","cell2","cell3","cell4"];

warning('off','MATLAB:unknownObjectNowStruct')

% Call scripts
scriptOrganizeAcquisitionsVCwithRandomSearch; % find corresponding sweeps for each stim
scriptCollectFeaturesVC; % calculate features of full field VC recordings
scriptCollectFeaturesRandomSearch; % same as above for random search

% separately run plotAnalysisRandomSearch