% Shun_figure_WiChR.m

%% Setup

clear; close all;
addpath(genpath(osPathSwitch('/Volumes/MICROSCOPE/Shun/Analysis/NeuroDAP/Methods')));
[~,~,~,~,~,~,bluePurpleRed] = loadColors;

% Define result directory
resultspath = osPathSwitch('/Volumes/MICROSCOPE/Shun/Project valence/Results');

fileList = uipickfiles('FilterSpec',osPathSwitch('/Volumes/MICROSCOPE/Shun/Project valence/Results'));
groupSessions = false;
% Update resultspath
dirsplit = strsplit(fileList{1},filesep); projectName = dirsplit{end-1}; 
resultspath = strcat(resultspath,filesep,projectName);
% Load selected files
for file = 1:length(fileList)
    dirsplit = strsplit(fileList{file},filesep);
    disp(['Ongoing: loading ',dirsplit{end}]);
    load(fileList{file});
    disp(['Finished: loaded ',dirsplit{end}]);
end

%% Baseline dLight: plot water, airpuff, stim, tone

timeRange = [-0.5,3];
eventRange = {'WiChR'};
animalRange = 'All';
taskRange = 'Random';
trialRange = 'All'; % range of trials in each session
totalTrialRange = 'All';
signalRange = {'dLight'};

colorList = {[.213 .543 .324]};
eventDuration = [.5,.5];

for s = 1:length(signalRange)
    initializeFig(.5,.5);
    legendEntries = {};
    for i = 1:length(eventRange)
        combined = combineTraces(animals,timeRange=timeRange,...
                                    eventRange=eventRange{i},...
                                    animalRange=animalRange,...
                                    taskRange=taskRange,...
                                    totalTrialRange=totalTrialRange,...
                                    trialRange=trialRange,...
                                    signalRange=signalRange{s});
        plotTraces(combined.data{1},combined.timestamp,color=colorList{i},plotShuffled=false);
        xlabel('Time (s)'); ylabel([signalRange{s},' z-score']); ylim([-0.5,1]);
        plotEvent(eventRange{i},eventDuration(i),color=colorList{i});
        legendEntries{end+1} = [eventRange{i},' (n=',num2str(size(combined.data{1},1)),')'];
    end
    legend(legendEntries, 'Location', 'northeast');
    % saveFigures(gcf,'Summary_random_dLight_WiChR',...
    %         strcat(resultspath),...
    %         saveFIG=true,savePDF=true);
end

%% dLight: Plot overall to show animal learned

timeRange = [-0.5,3];
eventRange = {'Stim','Pair','Tone','WiChR'};
animalRange = {'SL351','SL352','SL354','SL355','SL356'};
totalTrialRange = 'All';
trialRange = 'All';
signalRange = 'dLight';
trialConditions = 'trials.performing == 1';

colorList = {bluePurpleRed(500,:),bluePurpleRed(300,:),bluePurpleRed(100,:),[.213 .543 .324]};
groupSizeList = [30;30;10;20];
nGroupsList = [15;15;15;15];

taskRange = {'Reward','Punish'};
ylimList = [-1,4.5; -1,3.5];

for task = 1:length(taskRange)
    initializeFig(1,0.5); tiledlayout('flow');
    for event = 1:length(eventRange)
        nexttile;
        if event == 4
            combined = combineTraces(animals,timeRange=timeRange,...
                                        eventRange=eventRange{event},...
                                        animalRange=animalRange,...
                                        taskRange=taskRange{task},...
                                        totalTrialRange=totalTrialRange,...
                                        trialRange=trialRange,...
                                        signalRange=signalRange);
        else
            combined = combineTraces(animals,timeRange=timeRange,...
                                        eventRange=eventRange{event},...
                                        animalRange=animalRange,...
                                        taskRange=taskRange{task},...
                                        totalTrialRange=totalTrialRange,...
                                        trialRange=trialRange,...
                                        signalRange=signalRange,...
                                        trialConditions=trialConditions);
        end
        legendList = plotGroupTraces(combined.data{1},combined.timestamp,bluePurpleRed,...
                        groupSize=groupSizeList(event),nGroups=nGroupsList(event),...
                        groupby='sessions',startIdx=combined.options.startIdx,remaining='include');
        ylim(ylimList(task,:));
        plotEvent(eventRange{event},.5,color=colorList{event});
        xlabel('Time (s)'); ylabel([signalRange,' z-score']);
        legend(legendList,'Location','northeast');
    end
    % saveFigures(gcf,['Summary_dLight_WiChR_',taskRange{task}],...
    %     strcat(resultspath),...
    %     saveFIG=true,savePDF=true);
end


%% Plot grouped CS DA response (grouped across animal and 10 trials)

taskRange = {'Reward','Punish'};
statsTypes = {'stageAmp','stageAmp'}; ylabelList = {'Amp DA response during cue','Amp DA response during cue'};
ylimit = [-1,4; -1,5];

animalRange = {'SL351','SL352','SL354','SL355','SL356'};
conditionRange = 'All';
signalRange = 'dLight';
trialRange = 'All';
trialConditions = 'trials.performing == 1';

eventRange = {'WiChR','Baseline','Stim','Pair'};
colorList = {[.213 .543 .324],[0.8,0.8,0.8],bluePurpleRed(500,:),bluePurpleRed(300,:)};

combinedStats = getGroupedTrialStats(animals,statsTypes,...
                            eventRange=eventRange,...
                            animalRange=animalRange,...
                            taskRange=taskRange,...
                            trialRange=trialRange,...
                            totalTrialRange=conditionRange,...
                            signalRange=signalRange,...
                            trialConditions=trialConditions);

initializeFig(.7,.7); tiledlayout('flow');
results = plotGroupedTrialStats(combinedStats,ylabelList,groupSize=10,plotCommonTrials=false,...
        color=colorList,xlimIdx=3,xlim=[0,150],ylim=ylimit);

% saveFigures(gcf,'Summary_CSvsTrialsGrouped_dLight',...
%         strcat(resultspath),...
%         saveFIG=true,savePDF=true);

%% Plot bar plot of DA slopes
initializeFig(.7,.7); tiledlayout('flow');
for task = 1:length(results.stats)
    nexttile;
    cur_stats = results.stats{task};
    for event = 1:length(eventRange)
        slopes = cur_stats{event}(:,1);
        plotScatterBar(event,slopes,color=colorList{event},...
                       style='bar',dotSize=200,LineWidth=2,connectPairs=true);

        % Calculate significance
        if event < length(eventRange)
            for i = event+1:length(eventRange)
                plotStats(slopes,cur_stats{i}(:,1),[event i],testType='kstest');
            end
        end
    end

    xticks(1:length(eventRange)); xticklabels(eventRange);
    ylabel('Slope of DA amplitude during CS');
end

% saveFigures(gcf,'Summary_CSvsTrialsGrouped_dLight_slopeBar',...
%         strcat(resultspath),...
%         saveFIG=true,savePDF=true);