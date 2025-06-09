% Shun_figure_OptoPair.m

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


%% Plot overall to show animal learned

timeRange = [-0.5,3];
eventRange = 'Stim';
animalRange = 'SL043';
totalTrialRange = [1 200];
trialRange = 'All';
signalRange = 'dLight';
trialConditions = '';

colorList = bluePurpleRed(500,:);
groupSize = 30; nGroups = 5;

taskRange = {'Reward','Punish'};
ylimList = [-1.5,3; -2,1.5];

initializeFig(0.5,0.5); tiledlayout('flow');
for task = 1:length(taskRange)
    nexttile;
    combined = combineTraces(animals,timeRange=timeRange,...
                                eventRange=eventRange,...
                                animalRange=animalRange,...
                                taskRange=taskRange{task},...
                                totalTrialRange=totalTrialRange,...
                                trialRange=trialRange,...
                                signalRange=signalRange,...
                                trialConditions=trialConditions);
    legendList = plotGroupTraces(combined.data{1},combined.timestamp,bluePurpleRed,...
                    groupSize=groupSize,nGroups=nGroups,...
                    groupby='trials',startIdx=combined.options.startIdx,remaining='exclude');
    ylim(ylimList(task,:));
    plotEvent(eventRange,.5,color=colorList);
    xlabel('Time (s)'); ylabel([signalRange,' z-score']);
    legend(legendList,'Location','northeast');
    % saveFigures(gcf,['Summary_pairing_',taskRange{task}],...
    %     strcat(resultspath),...
    %     saveFIG=true,savePDF=true);
end


%% Plot grouped CS DA response (grouped across animal and 10 trials)

taskRange = {'Reward','Punish'};
statsTypes = {'stageAmp','stageAmp'}; ylabelList = {'Amp DA response during cue','Amp DA response during cue'};
ylimList = [-0.5,2.5; -0.5,2];

% animalRange = {'SL043','SL044','SL046','SL060','SL062','SL064','SL066','SL068'};%'All';
animalRange = 'all';
conditionRange = 'All';
signalRange = 'dLight';
trialConditions = 'trials.performing';

eventRange = {'Baseline','Stim','Pair'};
colorList = {[0.8,0.8,0.8],bluePurpleRed(500,:),bluePurpleRed(300,:)};

groupSize = 10; % numbers of trials to calculate average
combinedStats = getGroupedTrialStats(animals,statsTypes,...
                            eventRange=eventRange,...
                            animalRange=animalRange,...
                            taskRange=taskRange,...
                            totalTrialRange=conditionRange,...
                            signalRange=signalRange,...
                            trialConditions=trialConditions);

initializeFig(.7,.7); tiledlayout('flow');
results = plotGroupedTrialStats(combinedStats,ylabelList,groupSize=10,color=colorList,ylim=ylimList,...
                                xlimIdx=2,xlim=[1,150],plotIndividual=false);

% saveFigures(gcf,'Summary_CSvsTrialsGrouped_Amp',...
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
end
% 
% saveFigures(gcf,'Summary_CSvsTrialsGrouped_slopeBar',...
%         strcat(resultspath),...
%         saveFIG=true,savePDF=true);


%% Plot grouped anticipatory lick changes

taskRange = {'Reward','Punish'};
ylabelList = {'Anticipatory licks','Anticipatory licks'};
ylimList = [1.5,4.5; 0,4];

animalRange = 'all';
conditionRange = 'All';
signalRange = 'Lick';
trialConditions = 'trials.performing';

eventRange = {'Stim'};
colorList = {bluePurpleRed(500,:)};

groupSize = 10; % numbers of trials to calculate average
combinedStats = getGroupedTrialStats(animals,'nAnticipatoryLicks',...
                            eventRange=eventRange,...
                            animalRange=animalRange,...
                            taskRange=taskRange,...
                            totalTrialRange=conditionRange,...
                            signalRange=signalRange, ...
                            trialConditions=trialConditions);

initializeFig(.5,.5); tiledlayout('flow');
results = plotGroupedTrialStats(combinedStats,ylabelList,groupSize=10,color=colorList,...
                                xlim=[1,150],plotIndividual=false);

% saveFigures(gcf,'Summary_CSvsTrialsGrouped_Lick',...
%         strcat(resultspath),...
%         saveFIG=true,savePDF=true);


%% Plot grouped anticipatory lick changes

taskRange = {'Reward','Punish'};
ylabelList = {'Anticipatory licks','Anticipatory licks'};

animalRange = 'all';
conditionRange = 'All';
signalRange = 'Lick';
trialConditions = 'trials.performing';

eventRange = {'Baseline','Pair','Stim'};
colorList = {[.7 .7 .7],bluePurpleRed(300,:),bluePurpleRed(500,:)};

groupSize = 10; % numbers of trials to calculate average
combinedStats = getGroupedTrialStats(animals,'OutcomeTime',...
                            eventRange=eventRange,...
                            animalRange=animalRange,...
                            taskRange=taskRange,...
                            totalTrialRange=conditionRange,...
                            signalRange=signalRange, ...
                            trialConditions=trialConditions);

initializeFig(.5,.8); tiledlayout('flow');
results = plotGroupedTrialStats(combinedStats,ylabelList,groupSize=10,color=colorList,...
                                xlim=[1,150]);


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
    ylabel('Slope');
end
% saveFigures(gcf,'Summary_CSvsTrialsGrouped_Lick',...
%         strcat(resultspath),...
%         saveFIG=true,savePDF=true);


