%% Shun_loadSliceData
% Modified from Shun_analyzeSlice

% 09/13/23
% Separated from Shun_analyzeSlice, the idea is to plot individual and
% average trace from each epoch without referencing Excel data

% 09/14/23
% Package loading part and anlaysis part into separate function

%% Define data path
clear; close all;
addpath(genpath(osPathSwitch('/Volumes/MICROSCOPE/Shun/Analysis/NeuroDAP/Methods')));

% Select sessions for analysis
parentPath = osPathSwitch('/Volumes/MICROSCOPE/wengang/Exp_withShun/');
expPaths = uipickfiles('FilterSpec',parentPath,'Prompt','Select experiment folders');

% Set comman params
timeRange = [-20,50]; % in ms
nArtifactSamples = 10;
reload = false;

if length(expPaths) == 1; multipleSessions = false;
else; multipleSessions = true; end

%% Load epoch for single session

if ~multipleSessions
    epochs = loadSlices(expPaths{1},reload=reload,...
                filterSignal=false,filterSweeps=true,...
                nArtifactSamples=nArtifactSamples);
end

%% Save modified epochs if neccessary

if ~multipleSessions
    disp(expPaths{1});
    
    epochs = sortrows(epochs, [3 8]);
    expName = erase(expPaths{1},osPathSwitch(parentPath));
    save(strcat(expPaths{1},filesep,'epochs_',expName),'epochs','-v7.3');
    disp(strcat("Saved: ",expName));
    
    plotSliceEpochs(expPaths{1},...
                timeRange=timeRange,...
                nArtifactSamples=nArtifactSamples);
end

%% Plot multiple sessions

if multipleSessions
    groups = [1,2,2,2,2,1,1,0,0];
    
    plotSliceEpochs(expPaths,groups=groups,...
                timeRange=timeRange,nArtifactSamples=nArtifactSamples,...
                fitScatter=false);
end

return

%% Load and plot epoch for each sessions

if multipleSessions
    % Start analysis
    for i = 1:length(expPaths)
        close all
        try
            disp(strcat('********** ',erase(expPaths{i},parentPath),'**********'));
            epochs = loadSlices(expPaths{i},reload=reload,...
                filterSignal=false,filterSweeps=true,...
                nArtifactSamples=nArtifactSamples);
            plotSliceEpochs(expPaths{i},...
                timeRange=timeRange,...
                nArtifactSamples=nArtifactSamples);
        catch ME
            disp(getReport(ME));
            warning(['Session ', erase(expPaths{i},parentPath), ' have an error, skipped for now!!!!']);
            continue
        end
    end
end
