% Shun_figure_iGluSnFR

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

new_group = {'SL359','SL360','SL361','SL362','SL363'};
old_group = {'RCL','SL174','SL175','SL319','SL320','SL321','SL322','SL323'};

%% Baseline iGluSnFR: plot water, airpuff, stim, tone

timeRange = [-0.5,3];
eventRange = {'Airpuff','Stim'};
animalRange = old_group;
taskRange = 'Random';
trialRange = 'All'; % range of trials in each session
totalTrialRange = 'All';
signalRange = {'iGluSnFR'};

colorList = {[.2,.2,.2],bluePurpleRed(500,:),bluePurpleRed(100,:)};
eventDuration = [.1,.5,.5];

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
        xlabel('Time (s)'); ylabel([signalRange{s},' z-score']); ylim([-1,4]);
        plotEvent(eventRange{i},eventDuration(i),color=colorList{i});
        legendEntries{end+1} = [eventRange{i},' (n=',num2str(size(combined.data{1},1)),')'];
    end
    legend(legendEntries, 'Location', 'northeast');
    % saveFigures(gcf,'Summary_random_iGluSnFR',...
    %         strcat(resultspath),...
    %         saveFIG=true,savePDF=true);
end


ylabelList = {'Glu amplitude','Glu amplitude'};
eventRange = {'Stim','Airpuff'};
randomStats_Glu = getGroupedTrialStats(animals,'stageMax',...
                            eventRange=eventRange,...
                            animalRange='All',...
                            taskRange='Random',...
                            signalRange='iGluSnFR');
randomResults_Glu = plotGroupedTrialStats(randomStats_Glu,ylabelList,groupSize=1,plot=false);

colorList = {bluePurpleRed(100,:),bluePurpleRed(500,:)};
cur_randomTraces = randomResults_Glu.traces{1};
stimData = cur_randomTraces{1};
airpuffData = cur_randomTraces{2};

initializeFig(.5,.5); tiledlayout('flow');
nexttile;
stimData_animal = mean(stimData,2,'omitnan');
airpuffData_animal = mean(airpuffData,2,'omitnan');
plotScatterBar([1 2],[stimData_animal,airpuffData_animal],connectPairs=true,...
                LineWidth=5,dotSize=400,color=colorList,style='bar');
plotStats(stimData_animal,airpuffData_animal,[1,2],testType='kstest');
xticks(1:2); 
xticklabels({['Stim (n=',length(stimData_animal),')'],['Airpuff (',length(airpuffData_animal),')']});
ylabel('Amp iGluSnFR response during cue');

%% DA & iGluSnFR: Plot overall to show animal learned

timeRange = [-0.5,3];
eventRange = 'Stim';
animalRange = 'All';
totalTrialRange = 'All';
rewardTrialRange = [1,150];
punishTrialRange = [1,150];
trialConditions = '';

colorList = bluePurpleRed(500,:);
groupSizeList = 30;
nGroupsList = 15;

initializeFig(0.6,1); tiledlayout(2,2);
nexttile; 
signalRange = 'dLight';
taskRange = 'Reward';
combined = combineTraces(animals,timeRange=timeRange,...
                            eventRange=eventRange,...
                            animalRange=animalRange,...
                            taskRange=taskRange,...
                            totalTrialRange=totalTrialRange,...
                            trialRange=rewardTrialRange,...
                            signalRange=signalRange,...
                            trialConditions=trialConditions);
legendList = plotGroupTraces(combined.data{1},combined.timestamp,bluePurpleRed,...
                groupSize=groupSizeList,nGroups=nGroupsList,...
                groupby='trials',startIdx=combined.options.startIdx,remaining='include');
ylim([-1,3.5]);
plotEvent(eventRange,.5,color=colorList);
xlabel('Time (s)'); ylabel([signalRange,' z-score']);
legend(legendList,'Location','northeast');

nexttile; 
signalRange = 'iGluSnFR';
taskRange = 'Reward';
combined = combineTraces(animals,timeRange=timeRange,...
                            eventRange=eventRange,...
                            animalRange=animalRange,...
                            taskRange=taskRange,...
                            totalTrialRange=totalTrialRange,...
                            trialRange=rewardTrialRange,...
                            signalRange=signalRange,...
                            trialConditions=trialConditions);
legendList = plotGroupTraces(combined.data{1},combined.timestamp,bluePurpleRed,...
                groupSize=groupSizeList,nGroups=nGroupsList,...
                groupby='trials',startIdx=combined.options.startIdx,remaining='include');
ylim([-1,3]);
plotEvent(eventRange,.5,color=colorList);
xlabel('Time (s)'); ylabel([signalRange,' z-score']);
legend(legendList,'Location','northeast');

nexttile; 
signalRange = 'dLight';
taskRange = 'Punish';
combined = combineTraces(animals,timeRange=timeRange,...
                            eventRange=eventRange,...
                            animalRange=animalRange,...
                            taskRange=taskRange,...
                            totalTrialRange=totalTrialRange,...
                            trialRange=punishTrialRange,...
                            signalRange=signalRange,...
                            trialConditions=trialConditions);
legendList = plotGroupTraces(combined.data{1},combined.timestamp,bluePurpleRed,...
                groupSize=groupSizeList,nGroups=nGroupsList,...
                groupby='trials',startIdx=combined.options.startIdx,remaining='include');
ylim([-1,3]);
plotEvent(eventRange,.5,color=colorList);
xlabel('Time (s)'); ylabel([signalRange,' z-score']);
legend(legendList,'Location','northeast');

nexttile; 
signalRange = 'iGluSnFR';
taskRange = 'Punish';
combined = combineTraces(animals,timeRange=timeRange,...
                            eventRange=eventRange,...
                            animalRange=animalRange,...
                            taskRange=taskRange,...
                            totalTrialRange=totalTrialRange,...
                            trialRange=punishTrialRange,...
                            signalRange=signalRange,...
                            trialConditions=trialConditions);
legendList = plotGroupTraces(combined.data{1},combined.timestamp,bluePurpleRed,...
                groupSize=groupSizeList,nGroups=nGroupsList,...
                groupby='trials',startIdx=combined.options.startIdx,remaining='include');
ylim([-1,3]);
plotEvent(eventRange,.5,color=colorList);
xlabel('Time (s)'); ylabel([signalRange,' z-score']);
legend(legendList,'Location','northeast');


%% Plot grouped CS iGluSnFR response (grouped across animal and 10 trials)

taskRange = {'Reward','Punish'};
statsTypes = {'stageAmp','stageAmp'}; ylabelList = {'Glu amplitude','Glu amplitude'};
ylimit = [-1,3; -1,3];

animalRange = 'All';
conditionRange = 'All';
signalRange = 'iGluSnFR';
trialRange = [1,200];

% eventRange = {'Baseline','Pair','Tone','Stim'};
% colorList = {[0.8,0.8,0.8],bluePurpleRed(300,:),bluePurpleRed(100,:),bluePurpleRed(500,:)};
eventRange = {'Baseline','Stim','Pair'};
colorList = {[0.8,0.8,0.8],bluePurpleRed(500,:),bluePurpleRed(300,:)};

pairingStats_Glu = getGroupedTrialStats(animals,statsTypes,...
                            eventRange=eventRange,...
                            animalRange=animalRange,...
                            taskRange=taskRange,...
                            totalTrialRange=conditionRange,...
                            trialRange=trialRange,...
                            signalRange=signalRange);

initializeFig(.5,1); tiledlayout('flow');
randomResults_Glu = plotGroupedTrialStats(pairingStats_Glu,ylabelList,groupSize=10,color=colorList,...
                                          xlimIdx=2,xlim=[0,150],ylim=ylimit);


% Plot bar plot of iGluSnFR slopes
for task = 1:length(randomResults_Glu.stats)
    nexttile;
    cur_stats = randomResults_Glu.stats{task};
    for event = 1:length(eventRange)
        slopes = cur_stats{event}(:,1);
        plotScatterBar(event,slopes,color=colorList{event},...
                       style='bar',dotSize=200,LineWidth=2,connectPairs=true);
        ylim([-0.5,0.75]);

        % Calculate significance
        if event < length(eventRange)
            for i = event+1:length(eventRange)
                plotStats(slopes,cur_stats{i}(:,1),[event i],testType='kstest');
            end
        end
    end

    xticks(1:length(eventRange)); xticklabels(eventRange);
    ylabel('iGluSnFR slope');
end

% saveFigures(gcf,'Summary_CSvsTrialsGrouped_slopeBar',...
%         strcat(resultspath),...
%         saveFIG=true,savePDF=true);


%% Plot stim DA and iGluSnFR together


statsTypes = {'stageAmp','stageAmp'};
ylimit = [-0.5,3; -0.5,3];

animalRange = 'All';
conditionRange = 'All';

rewardStats_DA = getGroupedTrialStats(animals,statsTypes,...
                            eventRange={'Stim'},...
                            animalRange=animalRange,...
                            taskRange='Reward',...
                            totalTrialRange=conditionRange,...
                            signalRange='dLight');

rewardStats_Glu = getGroupedTrialStats(animals,statsTypes,...
                            eventRange={'Stim'},...
                            animalRange=animalRange,...
                            taskRange='Reward',...
                            totalTrialRange=conditionRange,...
                            signalRange='iGluSnFR');

punishStats_DA = getGroupedTrialStats(animals,statsTypes,...
                            eventRange={'Stim','Baseline'},...
                            animalRange=animalRange,...
                            taskRange='Punish',...
                            totalTrialRange=conditionRange,...
                            signalRange='dLight');

punishStats_Glu = getGroupedTrialStats(animals,statsTypes,...
                            eventRange={'Stim','Baseline'},...
                            animalRange=animalRange,...
                            taskRange='Punish',...
                            totalTrialRange=conditionRange,...
                            signalRange='iGluSnFR');


initializeFig(.5,.5); tiledlayout('flow');
rewardColor = {[0,158,115]./255,[0.7,0.7,0.7]};
punishColor = {[135,104,247]./255,[0.7,0.7,0.7]};

punishResults_dLight = plotGroupedTrialStats(punishStats_DA,'DA amplitude',groupSize=10,color=punishColor,xlim=[0,150],ylim=ylimit,plotNextTile=true);
rewardResults_dLight = plotGroupedTrialStats(rewardStats_DA,'DA amplitude',groupSize=10,color=rewardColor,xlim=[0,200],ylim=ylimit,plotNextTile=false);

punishResults_Glu = plotGroupedTrialStats(punishStats_Glu,'Glu amplitude',groupSize=10,color=punishColor,xlim=[0,150],ylim=ylimit,plotNextTile=true);
rewardResults_Glu = plotGroupedTrialStats(rewardStats_Glu,'Glu amplitude',groupSize=10,color=rewardColor,xlim=[0,200],ylim=ylimit,plotNextTile=false);

% saveFigures(gcf,'Summary_CSvsTrialsGrouped_stimOnlyCombined',...
%         strcat(resultspath),...
%         saveFIG=true,savePDF=true);


%% Plot iGluSnFR and DA slope together


initializeFig(.4,.5); tiledlayout('flow');

rewardSlopes_DA = rewardResults_dLight.stats{1}{1}(:,1);
rewardSlopes_Glu = rewardResults_Glu.stats{1}{1}(:,1);
punishSlopes_DA = punishResults_dLight.stats{1}{1}(:,1);
punishSlopes_Glu = punishResults_Glu.stats{1}{1}(:,1);
ctrlSlopes_DA = punishResults_dLight.stats{1}{2}(:,1);
ctrlSlopes_Glu = punishResults_Glu.stats{1}{2}(:,1);
rewardColor = [0,158,115]./255;
punishColor = [135,104,247]./255;
colorList = {[.7 .7 .7],rewardColor, punishColor};

% Plot DA
nexttile;
plotScatterBar(1,ctrlSlopes_DA,color=colorList{1},...
                   style='bar',dotSize=400,LineWidth=5,connectPairs=true);
plotScatterBar(2,rewardSlopes_DA,color=colorList{2},...
               style='bar',dotSize=400,LineWidth=5,connectPairs=true);
plotScatterBar(3,punishSlopes_DA,color=colorList{3},...
               style='bar',dotSize=400,LineWidth=5,connectPairs=true);
ylim([-0.5,0.5]);

plotStats(ctrlSlopes_DA,rewardSlopes_DA,[1 2],testType='kstest');
plotStats(ctrlSlopes_DA,punishSlopes_DA,[1 3],testType='kstest');
plotStats(rewardSlopes_DA,punishSlopes_DA,[2 3],testType='kstest');

xticks(1:3); xticklabels({'Ctrl','Glu','DA'});
ylabel('Slope');


% Plot reward DA vs Glu
nexttile;
plotScatterBar(1,ctrlSlopes_Glu,color=colorList{1},...
                   style='bar',dotSize=400,LineWidth=5,connectPairs=true);
plotScatterBar(2,rewardSlopes_Glu,color=colorList{2},...
               style='bar',dotSize=400,LineWidth=5,connectPairs=true);
plotScatterBar(3,punishSlopes_Glu,color=colorList{3},...
               style='bar',dotSize=400,LineWidth=5,connectPairs=true);

ylim([-0.5,0.75]);

plotStats(ctrlSlopes_Glu,rewardSlopes_Glu,[1 2],testType='kstest');
plotStats(ctrlSlopes_Glu,punishSlopes_Glu,[1 3],testType='kstest');
plotStats(rewardSlopes_Glu,punishSlopes_Glu,[2 3],testType='kstest');

xticks(1:3); xticklabels({'Ctrl','Reward','Punish'});
ylabel('Slope');