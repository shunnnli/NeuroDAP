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
% parentPath = osPathSwitch('/Volumes/MICROSCOPE/wengang/Exp_withShun/');
parentPath = osPathSwitch('/Volumes/MICROSCOPE/Shun/Project valence/Patch/');
expPaths = uipickfiles('FilterSpec',parentPath,'Prompt','Select experiment folders');
rawDataPath = 'default'; % strcat(parentPath,filesep,'20231221_ally');

% Set comman params
timeRange = [-20,50]; % in ms
nArtifactSamples = 10;
reload = true;
qc = true;

if isscalar(expPaths); multipleSessions = false;
else; multipleSessions = true; end

%% Load epoch for single session

if ~multipleSessions
    dirsplit = split(expPaths{1},filesep); expName = dirsplit{end};

    epochs = loadSlices(expPaths{1},reload=reload,...
                filterSignal=false,filterSweeps=true,...
                nArtifactSamples=nArtifactSamples,...
                rawDataPath=rawDataPath);

    plotSliceEpochs(expPaths{1},...
                timeRange=timeRange,...
                nArtifactSamples=nArtifactSamples);

    save(strcat(expPaths{1},filesep,'PreQC',filesep,'epochs_',expName),'epochs','-v7.3');
    disp(strcat("Saved: ",expName," in PreQC folder"));
    close all
end

% Message for editing quality checks
f = msgbox("Edit the epochs table and run following block after finished",...
    "Quality check","help");
return

%% Quality checks: save modified epoch file

% Message for editing quality checks
answer = questdlg('Confirm and save quality check results?', ...
    'Quality check confirmation','Yes','Not yet','Not yet');
switch answer
    case 'Yes'; qc = true;
    case 'Cake'; qc = false;
end

if ~multipleSessions && qc
    dirsplit = split(expPaths{1},filesep); expName = dirsplit{end};
    if qc
        save(strcat(expPaths{1},filesep,'epochs_',expName),'epochs','qc','-v7.3');
        disp(strcat("Saved: ",expName));
    end

    plotSliceEpochs(expPaths{1},...
                timeRange=timeRange,nArtifactSamples=nArtifactSamples,...
                resultPath='PostQC',plotAll=false);
    close all
end

%% Plot multiple sessions

% if multipleSessions
%     groups = [1,2,2,2,2,1,1,0,0];
% 
%     plotSliceEpochs(expPaths,groups=groups,...
%                 timeRange=timeRange,nArtifactSamples=nArtifactSamples,...
%                 fitScatter=false);
% end

return

%% Load and plot epoch for each sessions

if multipleSessions
    % Start analysis
    for i = 1:length(expPaths)
        try
            disp(strcat('********** ',erase(expPaths{i},parentPath),'**********'));
            epochs = loadSlices(expPaths{i},reload=reload,...
                filterSignal=false,filterSweeps=true,...
                nArtifactSamples=nArtifactSamples);
            plotSliceEpochs(expPaths{i},...
                timeRange=timeRange,...
                nArtifactSamples=nArtifactSamples,plotAll=false);
        catch ME
            disp(getReport(ME));
            warning(['Session ', erase(expPaths{i},parentPath), ' have an error, skipped for now!!!!']);
            continue
        end
        close all
    end
end

%% Useful code to plot raw sweeps

row = 19;
plotAll = false;

% Find event window
timeRange = [-20,50];
timeRangeStartSample = 10000 + 10000*timeRange(1)/1000;
timeRangeEndSample = 10000 + 10000*timeRange(2)/1000;
plotWindow = timeRangeStartSample : timeRangeEndSample;
timeRangeInms = (plotWindow-1*10000) ./ (10000/1000);
analysisWindow = (10000+nArtifactSamples)-timeRangeStartSample : length(plotWindow);

% Plot traces
if plotAll; included = 1:size(epochs{row,'Raw sweeps'}{1},1);
else; included = epochs{row,'Included'}{1}; end
traces = epochs{row,'Raw sweeps'}{1}(included==1,plotWindow);
plotSEM(timeRangeInms,traces,[0.343, 0.75, 0.232],...
        meanOnly=true,plotIndividual=true);
xlabel('Time (ms)');
ylabel('Current (pA)');
yMin = min(traces(:,analysisWindow),[],"all");
yMax = max(traces(:,analysisWindow),[],"all");
yPad = abs(yMax-yMin)*0.1;
ylim([yMin-yPad,yMax+yPad]);
title(strcat('Epochs #',num2str(epochs{row,'Epoch'})));