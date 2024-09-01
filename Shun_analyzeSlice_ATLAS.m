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
parentPath = osPathSwitch('/Volumes/MICROSCOPE/wengang/Exp_withAlly/');
expPaths = uipickfiles('FilterSpec',parentPath,'Prompt','Select experiment folders');

% Set comman params
timeRange = [-20,50]; % in ms
nArtifactSamples = 10;
reload = true;

if isscalar(expPaths); multipleSessions = false;
else; multipleSessions = true; end
[~,~,~,~,~,~,bluePurpleRed] = loadColors;

%% Load epoch for single session

if ~multipleSessions
    epochs = loadSlices(expPaths{1},reload=reload,...
                filterSignal=false,filterSweeps=true,...
                nArtifactSamples=nArtifactSamples,getCellTable=false,...
                defaultStimOnset=1000);
else
    combined = [];
    for i = 1:length(expPaths)
        strsplit = split(expPaths{i}); session = strsplit{end};
        disp(['Loading session: ',session]);
        epochs = loadSlices(expPaths{i},reload=reload,...
                filterSignal=false,filterSweeps=true,...
                nArtifactSamples=nArtifactSamples,getCellTable=false,...
                defaultStimOnset=1000);
        combined = [combined; epochs];
    end
end

% rawDataPath=strcat(parentPath,filesep,'20231221_ally')

%% Test: plot responses for selected epoch

close all;
row = [15 16];

timeRange = [-20,50]; % in ms
eventSample = 10000; outputFs = 10000;
timeRangeStartSample = eventSample + outputFs*timeRange(1)/1000;
timeRangeEndSample = eventSample + outputFs*timeRange(2)/1000;
plotWindow = timeRangeStartSample : timeRangeEndSample;


initializeFig(0.5,0.5); tiledlayout(1,2);
nexttile;
for i = 1:length(row)
    rowIdx = row(i);
    included = combined{rowIdx,'Included'}{1};
    data = combined{rowIdx,'Raw sweeps'}{1}(included==1,plotWindow);
    plotSEM(plotWindow,data,bluePurpleRed(i*100,:),plotIndividual=true,label=num2str(rowIdx));
end

nexttile;
for i = 1:length(row)
    rowIdx = row(i);
    included = combined{rowIdx,'Included'}{1};
    data = combined{rowIdx,'Processed sweeps'}{1}(included==1,plotWindow);
    plotSEM(plotWindow,data,bluePurpleRed(i*100,:),plotIndividual=true,label=num2str(rowIdx));
end
legend();

%% Optional: save modified epoch file
% 
% dirsplit = split(expPaths{1},filesep); expName = dirsplit{end};
% save(strcat(expPaths{1},filesep,'epochs_',expName),'epochs','-v7.3');
% disp(strcat("Saved: ",expName));

%% Load combined.mat

clear; close all;
addpath(genpath(osPathSwitch('/Volumes/MICROSCOPE/Shun/Analysis/NeuroDAP/Methods')));

% Select sessions for analysis
parentPath = osPathSwitch('/Volumes/MICROSCOPE/wengang/Exp_withAlly/Compiled_Data');
expPaths = uipickfiles('FilterSpec',parentPath,'Prompt','Select experiment folders');
resultsPath = '/Volumes/MICROSCOPE/wengang/Exp_withAlly/Compiled_Data/ATLAS-Cre';

if isscalar(expPaths); expPaths = expPaths{1}; end
load(expPaths);

%% Set comman params
timeRange = [-20,50]; % in ms
nArtifactSamples = 10;
reload = false;

% Set up
eventSample = 10000; outputFs = 10000;

timeRangeStartSample = eventSample + outputFs*timeRange(1)/1000;
timeRangeEndSample = eventSample + outputFs*timeRange(2)/1000;
plotWindow = timeRangeStartSample : timeRangeEndSample;
timeRangeInms = (plotWindow-1*outputFs) ./ (outputFs/1000);
analysisWindow = eventSample : timeRangeEndSample;

[~,~,~,~,~,~,bluePurpleRed] = loadColors;

% Extract traces
% Included
redIncluded = cell2mat(combined{combined.Red == 1,'Included'});
nonRedIncluded = cell2mat(combined{combined.Red == 0,'Included'});
% Traces
redTraces = cell2mat(combined{combined.Red == 1,'Processed sweeps'});
nonRedTraces = cell2mat(combined{combined.Red == 0,'Processed sweeps'});
% Remove not included traces
redTraces = redTraces(redIncluded,:);
nonRedTraces = nonRedTraces(nonRedIncluded,:);

% Calculate statistics
pairList = unique(combined.Pair);
pairPeaks = nan(length(pairList),2);
pairAUCs = nan(length(pairList),2);

for p = 1:length(pairList)
    % Find pair
    curPeaks = cellfun(@mean,combined{combined.Pair == pairList(p),'Peaks'});
    curAUCs = cellfun(@mean,combined{combined.Pair == pairList(p),'AUCs'});
    pairRedness = logical(combined{combined.Pair == pairList(p),'Red'});

    pairPeaks(p,1) = curPeaks(pairRedness);
    pairPeaks(p,2) = curPeaks(~pairRedness);
    pairAUCs(p,1) = curAUCs(pairRedness);
    pairAUCs(p,2) = curAUCs(~pairRedness);
end

pairAUCs = pairAUCs / 1e5; % change AUC unit to pC

% Calculate mean and sem of peak and AUCs
% Mean
avgRedPeaks = mean(pairPeaks(:,1));
avgNonRedPeaks = mean(pairPeaks(:,2));
avgRedAUCs = mean(pairAUCs(:,1));
avgNonRedAUCs = mean(pairAUCs(:,2));

% SEM
semRedPeaks = getSEM(pairPeaks(:,1));
semNonRedPeaks = getSEM(pairPeaks(:,2));
semRedAUCs = getSEM(pairAUCs(:,1));
semNonRedAUCs = getSEM(pairAUCs(:,2));

%% Optional: Calculate peaks

for i = 1:size(combined,1)
    included = combined{i,'Included'}{1};
    sweeps = combined{i,'Processed sweeps'}{1}(included==1,analysisWindow);

    combined{i,'Peaks'} = {min(sweeps,[],2)};
    combined{i,'AUCs'} = {sum(sweeps,2)};
end

%% Plot 1: average red vs non red

ylimit = [-400,50];

initializeFig(.5,.5); tiledlayout('flow');
nexttile;
% Plot red cells
plotTraces(redTraces(:,plotWindow),timeRangeInms,color=bluePurpleRed(end,:),...
        plotPatch=true,plotIndividual=true,...
        xlabel='Time (ms)', ylabel='Amplitude (pA)');
ylim(ylimit);

nexttile;
% Plot non red cells
plotTraces(nonRedTraces(:,plotWindow),timeRangeInms,color=[.2, .2, .2],...
        plotPatch=true,plotIndividual=true,...
        xlabel='Time (ms)', ylabel='Amplitude (pA)');
ylim(ylimit);

saveFigures(gcf,'Summary_redVsNonRed',resultsPath,savePDF=true,saveFIG=true);

%% Plot 2: mean of red vs non-red

ylimit = [-150,10];

initializeFig(.5,.5); tiledlayout('flow');

nexttile;
plotTraces(redTraces(:,plotWindow),timeRangeInms,plotPatch=true,color=bluePurpleRed(end,:),...
    xlabel='Time (ms)', ylabel='Amplitude (pA)');
plotTraces(nonRedTraces(:,plotWindow),timeRangeInms,plotPatch=true,color=[0.2,0.2,0.2],...
    xlabel='Time (ms)', ylabel='Amplitude (pA)');
ylim(ylimit);

nexttile;
plotTraces(redTraces(:,plotWindow),timeRangeInms,plotPatch=false,color=bluePurpleRed(end,:),...
    xlabel='Time (ms)', ylabel='Amplitude (pA)');
plotTraces(nonRedTraces(:,plotWindow),timeRangeInms,plotPatch=false,color=[0.2,0.2,0.2],...
    xlabel='Time (ms)', ylabel='Amplitude (pA)');
ylim(ylimit);

saveFigures(gcf,'Summary_redVsNonRed_mean',resultsPath,savePDF=true,saveFIG=true);

%% Plot 3: all pairs

ylimit = [-350,50];

initializeFig(1,1); tiledlayout('flow');

pairList = unique(combined.Pair);
for p = 1:length(pairList)
    % Find pair
    pairTraces = combined{combined.Pair == pairList(p),'Processed sweeps'};
    pairRedness = logical(combined{combined.Pair == pairList(p),'Red'});
    pairIncluded = combined{combined.Pair == pairList(p),'Included'};
    redTrace = pairTraces{pairRedness};
    nonRedTrace = pairTraces{~pairRedness};
    redIncluded = pairIncluded{pairRedness};
    nonRedIncluded = pairIncluded{~pairRedness};

    % Plot
    nexttile;
    plotTraces(redTrace(redIncluded==1,plotWindow),timeRangeInms,color=bluePurpleRed(end,:),...
        plotPatch=true,plotIndividual=true,...
        xlabel='Time (ms)', ylabel='Amplitude (pA)');
    plotTraces(nonRedTrace(nonRedIncluded==1,plotWindow),timeRangeInms,color=[0.2,0.2,0.2],...
        plotPatch=true,plotIndividual=true,...
        xlabel='Time (ms)', ylabel='Amplitude (pA)');
    % ylim(ylimit);
end

saveFigures(gcf,'Summary_redVsNonRed_allCell',resultsPath,savePDF=true,saveFIG=true);

%% Plot 4: bar graph

conditionLabels = {'Red','Non-red'};
initializeFig(.5,.5); tiledlayout('flow');

% Plot peaks
nexttile;
plot([1,2],pairPeaks,color=[.75,.75,.75]); hold on;
plotScatterBar(pairPeaks(:,1),1,color=bluePurpleRed(end,:),XJitterWidth=0.01,dotSize=100);
plotScatterBar(pairPeaks(:,2),2,color=[0.2,0.2,0.2],XJitterWidth=0.01,dotSize=100);

% [~,p,~] = kstest2(pairPeaks(:,1),pairPeaks(:,2));
p = signrank(pairPeaks(:,1),pairPeaks(:,2));
plotSignificance(p,[1 2],0.95);
xticks([1 2]); xticklabels(conditionLabels);
ylabel('Amplitude (pA)');

% Plot AUCs
nexttile;
plot([1,2],pairAUCs,color=[.75,.75,.75]); hold on;
plotScatterBar(pairAUCs(:,1),1,color=bluePurpleRed(end,:),XJitterWidth=0.01,dotSize=100);
plotScatterBar(pairAUCs(:,2),2,color=[0.2,0.2,0.2],XJitterWidth=0.01,dotSize=100);

% [~,p,~] = kstest2(pairAUCs(:,1),pairAUCs(:,2));
p = signrank(pairAUCs(:,1),pairAUCs(:,2));
plotSignificance(p,[1 2],0.95);
xticks([1 2]); xticklabels(conditionLabels);
ylabel('Total charge (pC)');

saveFigures(gcf,'Summary_amplitude_signrank',resultsPath,savePDF=true,saveFIG=true);