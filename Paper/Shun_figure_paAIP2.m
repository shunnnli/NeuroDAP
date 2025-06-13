% Shun_figure_paAIP2.m

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
eventRange = {'Stim','Pair','Tone'};
animalRange = 'All';
totalTrialRange = 'All';
trialRange = 'All';
signalRange = 'dLight';

colorList = {bluePurpleRed(500,:),bluePurpleRed(300,:),bluePurpleRed(100,:)};
groupSizeList = [50;50;10];
nGroupsList = [15;15;15];

taskRange = {'Reward1 (paAIP2)','Punish1 (paAIP2)'}; % first part
ylimList = [-0.75,2.5; -0.75,1.5; -1,2.5; -1.25,1.75]; % first part

% taskRange = {'Reward2 (Ctrl)'}; % second part
% ylimList = [-1.5,1; -1.5,1; -1.5,1; -1.25,1.75]; % Second part

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
taskRange = {'Reward1 (paAIP2)','Punish1 (paAIP2)'};
statsTypes = {'stageAmp','stageAmp'}; ylabelList = {'Amp DA response during cue','Amp DA response during cue'};
ylim = [-0.5,2; -0.5,3];

animalRange = 'All';

eventRange = {'Baseline','Stim','Pair'};
colorList = {[0.8,0.8,0.8],bluePurpleRed(500,:),bluePurpleRed(300,:)};

combinedStats = getGroupedTrialStats(animals,statsTypes,...
                            eventRange=eventRange,...
                            animalRange=animalRange,...
                            taskRange=taskRange,...
                            signalRange='dLight');

initializeFig(.7,.7); tiledlayout('flow');
results = plotGroupedTrialStats(combinedStats,ylabelList,groupSize=10,color=colorList,ylim=ylim,...
                                xlimIdx=2,xlim=[1,150]);

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
    if task == 2; ylim([-0.3,0.2]); end
end

% saveFigures(gcf,'Summary_CSvsTrialsGrouped_Pairing1_slopeBox',...
%         strcat(resultspath),...
%         saveFIG=true,savePDF=true);

%%

























%% (Supp) Plot grouped CS DA response (grouped across animal and 10 trials)

% First part
taskRange = {'Reward1 (paAIP2)','Punish1 (paAIP2)'};
% taskRange = {'Reward2 (Ctrl)','Punish2 (paAIP2)'};
% statsTypes = {'stageMax','stageMin'}; ylabelList = {'Max DA response during cue','Max DA response during cue'};
statsTypes = {'stageAmp','stageAmp'}; ylabelList = {'Amp DA response during cue','Amp DA response during cue'};
ylim = [-0.5,2; -0.5,3];

% Second part
% taskRange = {'Punish (paAIP2)'};
% statsTypes = {'stageAmp','stageAmp'}; ylabelList = {'Amp DA response during cue','Amp DA response during cue'};
% ylim = [-1,2.1; -1.5,2.5];

animalRange = 'All';
conditionRange = 'All';
signalRange = 'dLight';

% eventRange = {'Baseline','Pair','Tone','Stim'};
% colorList = {[0.8,0.8,0.8],bluePurpleRed(300,:),bluePurpleRed(100,:),bluePurpleRed(500,:)};
eventRange = {'Baseline','Stim','Pair'};
colorList = {[0.8,0.8,0.8],bluePurpleRed(500,:),bluePurpleRed(300,:)};

combinedStats = getGroupedTrialStats(animals,statsTypes,...
                            eventRange=eventRange,...
                            animalRange=animalRange,...
                            taskRange=taskRange,...
                            signalRange=signalRange);

initializeFig(.7,.7); tiledlayout('flow');
results = plotGroupedTrialStats(combinedStats,ylabelList,groupSize=10,color=colorList,ylim=ylim,...
                                xlimIdx=2,xlim=[1,150],plotIndividual=false);

% saveFigures(gcf,'Summary_CSvsTrialsGrouped_Pairing1_Amp',...
%         strcat(resultspath),...
%         saveFIG=true,savePDF=true);