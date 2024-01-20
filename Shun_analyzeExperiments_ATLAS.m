%% Shun_analyzeExperiments_ATLAS.m

% 2024/01/18

%% Load combined.mat

clear; close all;
addpath(genpath(osPathSwitch('/Volumes/MICROSCOPE/Shun/Analysis/NeuroDAP/Methods')));

% Select sessions for analysis
parentPath = osPathSwitch('/Volumes/MICROSCOPE/wengang/Exp_withAlly/');
expPaths = uipickfiles('FilterSpec',parentPath,'Prompt','Select experiment folders');
resultsPath = '/Volumes/MICROSCOPE/wengang/Exp_withAlly/Compiled_Data';

% Set comman params
timeRange = [-20,50]; % in ms
nArtifactSamples = 10;
reload = false;

if length(expPaths) == 1; expPaths = expPaths{1}; end
load('/Volumes/MICROSCOPE/wengang/Exp_withAlly/Compiled_Data/combined.mat');

% Set up
eventSample = 10000; outputFs = 10000;

timeRangeStartSample = eventSample + outputFs*timeRange(1)/1000;
timeRangeEndSample = eventSample + outputFs*timeRange(2)/1000;
plotWindow = timeRangeStartSample : timeRangeEndSample;
timeRangeInms = (plotWindow-1*outputFs) ./ (outputFs/1000);
analysisWindow = eventSample-timeRangeStartSample : length(plotWindow);

[~,~,~,~,~,~,bluePurpleRed] = loadColors;

%%
red = [];
for i = 1:height(combined)
    red = [red; combined{i,'Red'}{1}(1)];
end

%% Plot 1: average red vs non red

ylimit = [-350,50];

initializeFig(0.5,0.5); tiledlayout('flow');
nexttile;
% Plot red cells
traces = combined{combined.red == 1,'Processed sweeps'};
for i = 1:length(traces)
    plotTraces(traces{i}(:,plotWindow),timeRangeInms,meanOnly=true,color=bluePurpleRed(end,:),...
        xlabel='Time (ms)', ylabel='Amplitude (pA)');
end
ylim(ylimit);

nexttile;
% Plot non red cells
traces = combined{combined.red == 0,'Processed sweeps'};
for i = 1:length(traces)
    plotTraces(traces{i}(:,plotWindow),timeRangeInms,meanOnly=true,color=[.2, .2, .2],...
        xlabel='Time (ms)', ylabel='Amplitude (pA)');
end
ylim(ylimit);

saveFigures(gcf,'Summary_redVsNonRed',resultsPath,savePDF=true,saveFIG=true);

%% Plot 2: mean of red vs non-red

ylimit = [-200,50];

initializeFig(0.5,0.5); tiledlayout('flow');

redTraces = combined{combined.red == 1,'Processed sweeps'};
redMeanTraces = cell2mat(cellfun(@mean,redTraces,'UniformOutput',false));
nonRedTraces = combined{combined.red == 0,'Processed sweeps'};
nonRedMeanTraces = cell2mat(cellfun(@mean,nonRedTraces,'UniformOutput',false));

nexttile;
plotTraces(redMeanTraces(:,plotWindow),timeRangeInms,meanOnly=false,color=bluePurpleRed(end,:),...
    xlabel='Time (ms)', ylabel='Amplitude (pA)');
plotTraces(nonRedMeanTraces(:,plotWindow),timeRangeInms,meanOnly=false,color=[0.2,0.2,0.2],...
    xlabel='Time (ms)', ylabel='Amplitude (pA)');
ylim(ylimit);

nexttile;
plotTraces(redMeanTraces(:,plotWindow),timeRangeInms,meanOnly=true,color=bluePurpleRed(end,:),...
    xlabel='Time (ms)', ylabel='Amplitude (pA)');
plotTraces(nonRedMeanTraces(:,plotWindow),timeRangeInms,meanOnly=true,color=[0.2,0.2,0.2],...
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
    pairRedness = combined{combined.Pair == pairList(p),'red'};
    redTrace = pairTraces{find(pairRedness)};
    nonRedTrace = pairTraces{find(~pairRedness)};

    % Plot
    nexttile;
    plotTraces(redTrace(:,plotWindow),timeRangeInms,meanOnly=true,color=bluePurpleRed(end,:),...
        xlabel='Time (ms)', ylabel='Amplitude (pA)');
    plotTraces(nonRedTrace(:,plotWindow),timeRangeInms,meanOnly=true,color=[0.2,0.2,0.2],...
        xlabel='Time (ms)', ylabel='Amplitude (pA)');
    % ylim(ylimit);
end

saveFigures(gcf,'Summary_redVsNonRed_allCell',resultsPath,savePDF=true,saveFIG=true);

%% Plot 4: bar grapth

pairList = unique(combined.Pair);
pairPeaks = nan(length(pairList),2);
pairAUCs = nan(length(pairList),2);

for p = 1:length(pairList)
    % Find pair
    curPeaks = cellfun(@mean,combined{combined.Pair == pairList(p),'Peaks'});
    curAUCs = cellfun(@mean,combined{combined.Pair == pairList(p),'AUCs'});
    pairRedness = combined{combined.Pair == pairList(p),'red'};

    pairPeaks(p,1) = curPeaks(find(pairRedness));
    pairPeaks(p,2) = curPeaks(find(~pairRedness));
    pairAUCs(p,1) = curAUCs(find(pairRedness));
    pairAUCs(p,2) = curAUCs(find(~pairRedness));
end

%% Plot bar plot

conditionLabels = {'Red','Non-red'};
initializeFig(0.5,0.5); tiledlayout('flow');

% plot peaks
nexttile; 
boxchart(ones(length(pairList),1),pairPeaks(:,1),'BoxFaceColor',bluePurpleRed(end,:)); hold on
boxchart(ones(length(pairList),1)*2,pairPeaks(:,2),'BoxFaceColor',[0.2,0.2,0.2]);

plot([1,2],pairPeaks,color=[.75,.75,.75]);
swarmchart(ones(length(pairList),1),pairPeaks(:,1),150,bluePurpleRed(end,:),'filled',...
                'MarkerFaceAlpha',0.8,'XJitter','density','XJitterWidth',0.01); hold on
swarmchart(ones(length(pairList),1)*2,pairPeaks(:,2),150,[0.2,0.2,0.2],'filled',...
                'MarkerFaceAlpha',0.8,'XJitter','density','XJitterWidth',0.01); hold on

% [~,p,~] = kstest2(pairPeaks(:,1),pairPeaks(:,2));
p = signrank(pairPeaks(:,1),pairPeaks(:,2));
plotSignificance(p,[1 2],0.95);
xticks([1 2]); xticklabels(conditionLabels);
ylabel('Amplitude (pA)');


% plot AUCs
nexttile; 
boxchart(ones(length(pairList),1),pairAUCs(:,1),'BoxFaceColor',bluePurpleRed(end,:)); hold on
boxchart(ones(length(pairList),1)*2,pairAUCs(:,2),'BoxFaceColor',[0.2,0.2,0.2]);

plot([1,2],pairAUCs,color=[.75,.75,.75]);
swarmchart(ones(length(pairList),1),pairAUCs(:,1),150,bluePurpleRed(end,:),'filled',...
                'MarkerFaceAlpha',0.8,'XJitter','density','XJitterWidth',0.01); hold on
swarmchart(ones(length(pairList),1)*2,pairAUCs(:,2),150,[0.2,0.2,0.2],'filled',...
                'MarkerFaceAlpha',0.8,'XJitter','density','XJitterWidth',0.01); hold on

% [~,p,~] = kstest2(pairAUCs(:,1),pairAUCs(:,2)); 
p = signrank(pairAUCs(:,1),pairAUCs(:,2)); 
plotSignificance(p,[1 2],0.95);
xticks([1 2]); xticklabels(conditionLabels);
ylabel('AUC (A.U.)');

saveFigures(gcf,'Summary_amplitude_signrank',resultsPath,savePDF=true,saveFIG=true);
