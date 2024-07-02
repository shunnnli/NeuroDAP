% Shun_analyzeSlice_OptoPair
% 06/14/24

%% Load sessions
clear; close all;
addpath(genpath(osPathSwitch('/Volumes/MICROSCOPE/Shun/Analysis/NeuroDAP/Methods')));

% Select slice sessions for analysis
% parentPath = osPathSwitch('/Volumes/MICROSCOPE/wengang/Exp_withShun/');
parentPath = osPathSwitch('/Volumes/MICROSCOPE/Shun/Project valence/Patch/');
expPaths = uipickfiles('FilterSpec',parentPath,'Prompt','Select experiment folders');
rawDataPath = 'default'; % strcat(parentPath,filesep,'20231221_ally');

% Select corresponding photometry folder
parentPath_photometry = osPathSwitch('/Volumes/MICROSCOPE/Shun/Project valence/Recordings/');

% Set comman params
timeRange = [-20,100]; % in ms
nArtifactSamples = 10;

%% Concatenate epochs.mat from selected sessions

% Initialize empty allEpochs
allEpochs = [];

% Concatenate epochs
for i = 1:length(expPaths)
    dirsplit = split(expPaths{i},filesep); expName = dirsplit{end};
    sessionEpoch = loadSlices(expPaths{i},reload=false);
    allEpochs = [allEpochs; sessionEpoch];
    disp(strcat('Finished concatenating session: ',expName, ...
        ' (',num2str(i),'/',num2str(length(expPaths)),')'));
end

%% Find corresponding photometry sessions and concatenate those

animalList = unique(allEpochs{:,"Animal"});
allStimTraces = [];

for i = 1:length(animalList)
    % Find corresponding physiology session
    animalRows = allEpochs.Animal == animalList(i);
    sessionPath = unique(allEpochs{animalRows,'Session'});
    if length(sessionPath) ~= 1
        if isempty(sessionPath)
            error(['Cannot find corresponding physiology session for animal',animalList(i)]); 
        else
            warning([animalList(i),': multiple physiology session detected, choose the first one by default.']);
            sessionPath = sessionPath(1); 
        end
    end
    dirsplit = split(sessionPath,filesep); sessionName = dirsplit{end}; sessionDate = sessionName(1:end-5);
    disp(strcat(animalList(i),": found physiology session ",sessionName));

    % Find corresponding photometry session
    photometryName = char(strcat(sessionDate,'-',animalList(i)));
    photometryPath = dir(fullfile(parentPath_photometry,'*',[photometryName,'*']));
    if length(photometryPath) ~= 1
        if isempty(photometryPath)
            error(['Cannot find corresponding photometry session for animal', animalList(i)]);
        else
            warning([animalList(i),': multiple photometry session detected, choose the first one by default.']);
            photometryPath = photometryPath(1); 
        end
    end

    % Open analysis.mat to extract stim response
    load(strcat(photometryPath.folder,filesep,photometryPath.name,filesep,'analysis_',photometryPath.name,'.mat'));
    allStimTraces = [allStimTraces, analysis(contains({analysis.event},'stim',IgnoreCase=true))];
    disp(strcat(animalList(i),": found photometry session ",photometryPath.name));
end

%% Extract EPSC/IPSC ratio


%% Plot EPSC/IPSC ratio against DA response



%% Plot multiple sessions