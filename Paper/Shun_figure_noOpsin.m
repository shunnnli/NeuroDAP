% Shun_figure_noOpsin.m

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


%% Baseline DA: plot water, airpuff, stim, tone

timeRange = [-0.5,3];
eventRange = {'Water','Airpuff','Stim','Tone'};
animalRange = 'All';
taskRange = 'Random';
trialRange = 'All'; % range of trials in each session
totalTrialRange = 'All';
signalRange = {'NAc'};

colorList = {bluePurpleRed(1,:),[.2,.2,.2],bluePurpleRed(500,:),bluePurpleRed(100,:)};
eventDuration = [0,.1,.5,.5];

for s = 1:length(signalRange)
    initializeFig(.5,.5); tiledlayout(2,2);
    for i = 1:length(eventRange)
        nexttile;
        combined = combineTraces(animals,timeRange=timeRange,...
                                    eventRange=eventRange{i},...
                                    animalRange=animalRange,...
                                    taskRange=taskRange,...
                                    totalTrialRange=totalTrialRange,...
                                    trialRange=trialRange,...
                                    signalRange=signalRange{s});
        plotTraces(combined.data{1},combined.timestamp,color=colorList{i});
        xlabel('Time (s)'); ylabel([signalRange{s},' z-score']); ylim([-1,4]);
        plotEvent(eventRange{i},eventDuration(i),color=colorList{i});
        legend([eventRange{i},' (n=',num2str(size(combined.data{1},1)),')'],...
                'Location','northeast');
    end
    % saveFigures(gcf,'Summary_random_dLight',...
    %         strcat(resultspath),...
    %         saveFIG=true,savePDF=true);
end

%% Plot overall to show animal learned

timeRange = [-0.5,3];
eventRange = {'Stim','Pair'};
animalRange = 'All';
totalTrialRange = 'All';
trialRange = 'All';
signalRange = 'NAc';

colorList = {bluePurpleRed(500,:),bluePurpleRed(300,:),bluePurpleRed(100,:)};
groupSizeList = [50;50;10];
nGroupsList = [15;15;15];
taskRange = {'Reward1','Punish1'};

for task = 1:length(taskRange)
    initializeFig(1,0.5); tiledlayout('flow');
    for event = 1:length(eventRange)
        nexttile;
        combined = combineTraces(animals,timeRange=timeRange,...
                                    eventRange=eventRange{event},...
                                    animalRange=animalRange,...
                                    taskRange=taskRange{task},...
                                    totalTrialRange=totalTrialRange,...
                                    trialRange=trialRange,...
                                    signalRange=signalRange);
        legendList = plotGroupTraces(combined.data{1},combined.timestamp,bluePurpleRed,...
                        groupSize=groupSizeList(event),nGroups=nGroupsList(event),...
                        groupby='trials',startIdx=combined.options.startIdx,remaining='include');
        % ylim(ylimList(task,:));
        plotEvent(eventRange{event},.5,color=colorList{event});
        xlabel('Time (s)'); ylabel([signalRange,' z-score']);
        legend(legendList,'Location','northeast');
    end
    % saveFigures(gcf,['Summary_pairing_',taskRange{task}],...
    %     strcat(resultspath),...
    %     saveFIG=true,savePDF=true);
end

%% Plot grouped CS DA response (grouped across animal and 10 trials)

% First part
taskRange = {'Reward1','Punish1'};
statsTypes = {'stageMax','stageMin'}; ylabelList = {'DA amplitude','DA amplitude'};
% ylim = [-1.5,3; -1.5,3];

animalRange = 'All';

eventRange = {'Baseline','Stim','Pair'};
colorList = {[0.8,0.8,0.8],bluePurpleRed(500,:),bluePurpleRed(300,:)};

combinedStats = getGroupedTrialStats(animals,statsTypes,...
                            eventRange=eventRange,...
                            animalRange=animalRange,...
                            taskRange=taskRange,...
                            signalRange='NAc');

initializeFig(.7,.7); tiledlayout('flow');
results = plotGroupedTrialStats(combinedStats,ylabelList,groupSize=10,color=colorList,...
                                xlimIdx=2,xlim=[1,150],plotCommonTrials=true);

% saveFigures(gcf,'Summary_CSvsTrialsGrouped_Pairing1_Amp',...
%         strcat(resultspath),...
%         saveFIG=true,savePDF=true);

%% Plot bar plot of DA slopes
initializeFig(.7,.7); tiledlayout('flow'); clearvars ylim
for task = 1:length(results.stats)
    nexttile;
    cur_stats = results.stats{task};
    for event = 1:length(eventRange)
        slopes = cur_stats{event}(:,1);
        plotScatterBar(event,slopes,color=colorList{event},...
                       style='bar',dotSize=200,LineWidth=2);

        % Calculate significance
        if event < length(eventRange)
            for i = event+1:length(eventRange)
                plotStats(slopes,cur_stats{i}(:,1),[event i],testType='kstest');
            end
        end
    end

    xticks(1:length(eventRange)); xticklabels(eventRange);
    ylabel('Slope of DA amplitude during CS');
    % ylim([-1,1]);
end

% saveFigures(gcf,'Summary_CSvsTrialsGrouped_Pairing1_slopeBox',...
%         strcat(resultspath),...
%         saveFIG=true,savePDF=true);