%% Shun_figure_DAclamp

clear; close all;
addpath(genpath(osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Analysis/NeuroDAP/Methods')));
[~,~,~,~,~,~,bluePurpleRed] = loadColors;

% Define result directory
resultspath = osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Project clamping/');

% Load animals.mat files
load('path_to_file');

%% Set up animals

clampAnimals = {'SL287'};
unclampAnimals = setdiff(unique({animals.animal}),clampAnimals);

%% dLight: clampped vs non clampped response to water

timeRange = [-0.5,3];
eventRange = {'Rewarded Licks','Water'};
taskRange = 'Reward';
trialRange = 'All'; % range of trials in each session
totalTrialRange = 'All';
signalRange = {'dLight'};

ctrlColorList = {[.7 .7 .7],[.7 .7 .7],[.7 .7 .7],[.7 .7 .7]};
pidColorList = {bluePurpleRed(1,:),bluePurpleRed(100,:)};

eventDuration = [0,0,.5,.5];

for s = 1:length(signalRange)
    initializeFig(.5,.5); tiledlayout(1,2);
    for i = 1:length(eventRange)
        nexttile;
        combined_ctrl = combineTraces(animals,timeRange=timeRange,...
                                    eventRange=eventRange{i},...
                                    animalRange=clampAnimals,...
                                    taskRange=taskRange,...
                                    totalTrialRange=totalTrialRange,...
                                    trialRange=trialRange,...
                                    signalRange=signalRange{s});
        plotTraces(combined_ctrl.data{1},combined_ctrl.timestamp,color=pidColorList{i},plotShuffled=false);

        combined_pid = combineTraces(animals,timeRange=timeRange,...
                                    eventRange=eventRange{i},...
                                    animalRange=unclampAnimals,...
                                    taskRange=taskRange,...
                                    totalTrialRange=totalTrialRange,...
                                    trialRange=trialRange,...
                                    signalRange=signalRange{s});
        plotTraces(combined_pid.data{1},combined_pid.timestamp,color=ctrlColorList{i},plotShuffled=false);

        xlabel('Time (s)'); ylabel([signalRange{s},' z-score']); %ylim([-2,4]);
        plotEvent(eventRange{i},eventDuration(i),color=pidColorList{i});
        legend({['PID (n=',num2str(size(combined_pid.data{1},1)),')'],...
                ['Ctrl (n=',num2str(size(combined_ctrl.data{1},1)),')']},...
                'Location','northeast');
    end
    % saveFigures(gcf,'Summary_random_NAc',...
    %         strcat(resultspath),...
    %         saveFIG=true,savePDF=true);
end

%% dLight: Plot overall to show animal learned

timeRange = [-0.5,3];
eventRange = {'Tone'};
totalTrialRange = 'All';
trialRange = 'All';
signalRange = 'dLight';

colorList = {bluePurpleRed(500,:),bluePurpleRed(300,:),bluePurpleRed(100,:)};
groupSizeList = [50;50;10];
nGroupsList = [15;15;15];

taskRange = {'Reward'};
ylimList = [-1,2.5; -1,2.5; -1.5,1; -1.25,1.75];

for task = 1:length(taskRange)
    initializeFig(0.5,0.5); tiledlayout('flow');
    for event = 1:length(eventRange)
        nexttile;
        combined_pid = combineTraces(animals,timeRange=timeRange,...
                                    eventRange=eventRange{event},...
                                    animalRange=clampAnimals,...
                                    taskRange=taskRange{task},...
                                    totalTrialRange=totalTrialRange,...
                                    trialRange=trialRange,...
                                    signalRange=signalRange);
        legendList = plotGroupTraces(combined_pid.data{1},combined_pid.timestamp,bluePurpleRed,...
                        groupSize=groupSizeList(event),nGroups=nGroupsList(event),...
                        groupby='trials',startIdx=combined_pid.options.startIdx,remaining='include');
        ylim(ylimList(task,:));
        plotEvent(eventRange{event},.5,color=colorList{event});
        xlabel('Time (s)'); ylabel([signalRange,' z-score']);
        legend(legendList,'Location','northeast');

        nexttile;
        combined_ctrl = combineTraces(animals,timeRange=timeRange,...
                                    eventRange=eventRange{event},...
                                    animalRange=unclampAnimals,...
                                    taskRange=taskRange{task},...
                                    totalTrialRange=totalTrialRange,...
                                    trialRange=trialRange,...
                                    signalRange=signalRange);
        legendList = plotGroupTraces(combined_ctrl.data{1},combined_ctrl.timestamp,bluePurpleRed,...
                        groupSize=groupSizeList(event),nGroups=nGroupsList(event),...
                        groupby='trials',startIdx=combined_ctrl.options.startIdx,remaining='include');
        ylim(ylimList(task,:));
        plotEvent(eventRange{event},.5,color=colorList{event});
        xlabel('Time (s)'); ylabel([signalRange,' z-score']);
        legend(legendList,'Location','northeast');
    end
    % saveFigures(gcf,['Summary_pairing_',taskRange{task}],...
    %     strcat(resultspath),...
    %     saveFIG=true,savePDF=true);
end

%% Licking: 

timeRange = [-0.5,3];
eventRange = {'Tone'};
totalTrialRange = 'All';
trialRange = 'All';
signalRange = 'Lick';

colorList = {bluePurpleRed(500,:),bluePurpleRed(300,:),bluePurpleRed(100,:)};
groupSizeList = [30;30;30];
nGroupsList = [5;5;5];

% Second part
taskRange = {'Reward'}; 
ylimList = [0,3; 0,3; -1.5,1; -1.25,1.75]; % Second part

for task = 1:length(taskRange)
    initializeFig(.5,.5); tiledlayout('flow');
    for event = 1:length(eventRange)
        nexttile;
        combined = combineTraces(animals,timeRange=timeRange,...
                                    eventRange=eventRange{event},...
                                    animalRange=clampAnimals,...
                                    taskRange=taskRange{task},...
                                    totalTrialRange=totalTrialRange,...
                                    trialRange=trialRange,...
                                    signalRange=signalRange);
        legendList = plotGroupTraces(combined.data{1},combined.timestamp,bluePurpleRed,...
                        groupSize=groupSizeList(event),nGroups=nGroupsList(event),...
                        groupby='sessions',startIdx=combined.options.startIdx);
        ylim(ylimList(task,:));
        plotEvent(eventRange{event},.5,color=colorList{event});
        xlabel('Time (s)'); ylabel('Licks/s'); ylim([0 Inf]);
        legend(legendList,'Location','northeast');

        nexttile;
        combined = combineTraces(animals,timeRange=timeRange,...
                                    eventRange=eventRange{event},...
                                    animalRange=unclampAnimals,...
                                    taskRange=taskRange{task},...
                                    totalTrialRange=totalTrialRange,...
                                    trialRange=trialRange,...
                                    signalRange=signalRange);
        legendList = plotGroupTraces(combined.data{1},combined.timestamp,bluePurpleRed,...
                        groupSize=groupSizeList(event),nGroups=nGroupsList(event),...
                        groupby='sessions',startIdx=combined.options.startIdx);
        ylim(ylimList(task,:));
        plotEvent(eventRange{event},.5,color=colorList{event});
        xlabel('Time (s)'); ylabel('Licks/s'); ylim([0 Inf]);
        legend(legendList,'Location','northeast');
    end
    % saveFigures(gcf,['Summary_licking_',taskRange{task}],...
    %     strcat(resultspath),...
    %     saveFIG=true,savePDF=true);
end