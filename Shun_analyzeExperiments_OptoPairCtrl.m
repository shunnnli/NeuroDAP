% Shun_analyzeExperiments_template
% 2023/12/04

% Template for multi-session analysis of an experiment

% To start, copy this matlab file and replace template with specific
% experiments. In theory, there will be a specific analyzeExperiments file
% for each individual experiments, as the specific needs for analysis
% varies between different experiments.

%% Analysis pipeline
% The pipeline in general is the following:

% 1. Select whether to load a previously animals struct (described below)
% or select individual session to combine.

% 2. After selecting ALL SESSIONS from an experiments, the pipeline will
% automatically concatenate analysis.mat for each recording sessions.
% Rename properties as needed in order to facilitate further analysis.

% 3. Run getAnimalStruct.m function to recreate animals struct from
% summary. This combines all sessions from the same animals together while
% cutoffs between individual sessions are also recorded.

% 4. Save animals struct if needed. Note: saving summary struct will take
% extremely long (>5hrs) so while saving animals struct is much shorter 
% (~2min). animals struct should contain information that satisfies MOST 
% plotting requirements so saving summary struct is not needed.

% 5. Data analysis and plotting. This part is designed to vary across
% experiments. Thus, following codes are just for demonstration of 
% essential functions.

%% Essential functions

% getAnimalStruct(summary)
% combine sessions of the same animal, from the same task,
% of the same event, recorded from the same signal (eg NAc, LHb, cam,
% Lick) together. As described above, animals struct will be the MOST
% IMPORTANT struct that stores information about the experiments.

% combineTraces(animals,options)
% combine traces and their relevant statstics of selected animals, 
% selected tasks, selected trialRange, selected totalTrialRange, 
% selected signals, and selected events together. 
% This is the MOST IMPORTANT and USED function in this script. 
% Important features are listed as follows:
    % 1. the function returns a structure with fields. data fields stores
    % the data (photometry, cam, lick rate traces) of selected sessions.
    % 2. Field stats stores stageAvg/Max/Min of each traces at selected stage
    % time (often determined when creating analysis.mat but can modify later).
    % 3. Field options contains following important variables:
        % options.empty: true if no session is found that fits the input criteria. 
        % Should skip during plotting or further analysis 
        % options.animalStartIdx: Records index (in field data) of the
        % first trace for each animals. Used in plotGroupTraces
        % options.sessionStartIdx: Records index (in field data) of the
        % first trace for each session.
    % 4. totalTrialRange and trialRange
        % totalTrialRange selects the ACTUAL trial number within each
        % session while trialRange selects the samples across selected
        % sessions. For example: I have 3 session where I inhibit CaMKII
        % activity for the first 60 trials of each session. Within these
        % first 60 trials, 30% of them are stim-only trials. If I want to
        % only plot the 50-100th stim-only trials with CaMKII inhibition 
        % across all sessions, I will set totalTrialRange=[1,60] and
        % trialRange=[50,100]. Detailed description and automatic handling
        % of edge cases is documented within the method.
    % 5. The function can take both animals and summary struct as inputs.

% plotGroupTraces(combined.data,combined.timestamp,options)
% While plotGroupTraces is also used in analyzeSessions.m; here, we can
% plot traces across all animal easily (see code below). Key options are as
% follows:
    % 1. groupSize and nGroups
        % You need to provide either groupSize or nGroups for the function to
        % run. If you provide both, plotGroupTraces will plot to the maximum
        % number of groups based on groupSize. Thus, for a input with 50
        % trials and groupSize = 10, the function will automatically plot 5
        % lines even when nGroups=10
    % 2. options.animalStartIdx
        % Use this to reorganize input data so that its plotted based on
        % animals. eg when I want to plot Trial 1-10, 11-20 for EACH animal
        % across all sessions
    % 3. options.remaining
        % There inevitably will be some traces that does not fully form a group
        % (eg 5 traces remaining for a groupSize of 50 traces). These traces,
        % if plotted separately, can induce lines with great variations and
        % error bars. To address this, one can either set remaining='include'
        % to include these traces to prev group; set to 'exclude' to not plot
        % these traces, or 'separate' if you really want to plot these traces
        % separately

%% Intro of sample data set

% The sample data set is recorded by Shun Li in 2023. It contains 4
% animals, with 1 animals with off-target expression ('SL137'). dLight
% signals in NAc, pupil/Eye area, and lick are simultaneously recorded for
% all sessions.

% There are 5 major phases:
    % 1. Random: water, airpuff, EP stim, and tone (75dB) are delivered
    % randomly
    % 2. Reward1/2: where EP stim and tone are paired with water
    % 3. Punish1/2: where EP stim and tone are paired with airpuff
    % 4. Timeline: Random (2 sessions) -> Reward1 (3 sessions) -> Punish1
    % (3 sessions) -> Reward2 (3 sessions) -> 1 week rest -> Punish2 (3 sessions but 3 animals)

%% Setup

clear; close all;
addpath(genpath(osPathSwitch('/Volumes/MICROSCOPE/Shun/Analysis/NeuroDAP/Methods')));
[~,~,~,~,~,~,bluePurpleRed] = loadColors;
today = char(datetime('today','Format','yyyyMMdd'));

% Define result directory
resultspath = osPathSwitch('/Volumes/MICROSCOPE/Shun/Project valence/Results');

% Building summary struct from selected sessions
answer = questdlg('Group sessions or load combined data?','Select load sources',...
                  'Group single sessions','Load combined data','Load sample data','Load combined data');

if strcmpi(answer,'Group single sessions')
    sessionList = uipickfiles('FilterSpec',osPathSwitch('/Volumes/MICROSCOPE/Shun/Project valence/Recordings'));
    groupSessions = true;
    % Update resultspath
    dirsplit = strsplit(sessionList{1},filesep); projectName = dirsplit{end-1}; 
    resultspath = strcat(resultspath,filesep,projectName);
    % Create resultspath if necessary
    if isempty(dir(resultspath)); mkdir(resultspath); end

elseif strcmpi(answer,'Load combined data')
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

elseif strcmpi(answer,'Load sample data')
    groupSessions = false;
    % Update resultspath
    dirsplit = strsplit(fileList{1},filesep); projectName = dirsplit{end-1}; 
    resultspath = strcat(osPathSwith('/Volumes/MICROSCOPE/Shun/Analysis/NeuroDAP/Tutorials/Sample data/Results'),filesep,projectName);
    % Load selected files
    for file = 1:length(fileList)
        dirsplit = strsplit(fileList{file},filesep);
        disp(['Ongoing: loading ',dirsplit{end}]);
        load(fileList{file});
        disp(['Finished: loaded ',dirsplit{end}]);
    end
end

%% Optional: Create summary struct (only need to do this for initial loading)

if groupSessions    
    summary = concatAnalysis(sessionList);
    trialTables = loadTrialTables(sessionList);
end

% Check summary format (should all be chars NOT strings)
stringColumnsLabels = {'animal','date','session','task','event','name','system'};
for i = 1:length(stringColumnsLabels)
    for row = 1:length(summary)
        if isstring(summary(row).(stringColumnsLabels{i})) 
            summary(row).(stringColumnsLabels{i}) = convertStringsToChars(summary(row).(stringColumnsLabels{i}));
        end
    end
end
disp('Finished: summary struct and trialtables loaded');

%% Optional: Make changes to summary for further analysis

% for i = 1:211; summary(i).task = 'Random'; end
% for i = 212:499; summary(i).task = 'Reward1'; end
% for i = 500:919; summary(i).task = 'Punish1'; end
% for i = 920:1279; summary(i).task = 'Reward2'; end
% for i = 1280:1594; summary(i).task = 'Punish2'; end

for i = 1:length(summary)
    cur_task = summary(i).task;
    if contains('random',cur_task,IgnoreCase=true)
        summary(i).task = 'Random';
    elseif contains('reward pairing',cur_task,IgnoreCase=true)
        summary(i).task = 'Reward1';
    elseif contains('punish pairing',cur_task,IgnoreCase=true)
        summary(i).task = 'Punish1';
    end
end

%% Create animals struct

if isempty(dir(fullfile(resultspath,'animals*.mat'))) || groupSessions
    animals = getAnimalsStruct(summary);
end

%% Save animals struct

% Save animals.mat
disp(['Ongoing: saving animals.mat (',char(datetime('now','Format','HH:mm:ss')),')']);
save(strcat(resultspath,filesep,'animals_',today),'animals','trialTables','sessionList','-v7.3');
disp(['Finished: saved animals.mat (',char(datetime('now','Format','HH:mm:ss')),')']);

% Save summary.mat (not recommend!! Will take forever!!)
% disp(['Ongoing: saving summary.mat (',char(datetime('now','Format','HH:mm:ss')),')']);
% save(strcat(resultspath,filesep,'summary_',today),'summary','sessionList','-v7.3');
% disp(['Finished: saved summary.mat (',char(datetime('now','Format','HH:mm:ss')),')']);

%% Test: Plot traces from summary/animals struct

initializeFig(0.5,0.5);
combined = combineTraces(animals,timeRange=[-0.5,3],...
                            eventRange='Water',...
                            animalRange="SL160",...
                            taskRange='Punish1',...
                            totalTrialRange='All',...
                            trialRange='All',...
                            signalRange='LHb',...
                            sessionRange='20240101-SL161-P4');
plotTraces(combined.data{1},combined.timestamp,color=bluePurpleRed(1,:));
plotEvent('Water',0,color=bluePurpleRed(1,:))
xlabel('Time (s)'); ylabel('z-score');

%% Optional (calculate LHb weird oscillation frequency)
endPoint = min(1e6, length(combined.data{1}(2,:)));
sData_fft = fft(normalize(trace(1:endPoint)));
P2 = abs(sData_fft/endPoint);
P1 = P2(1:endPoint/2+1);
P1(2:end-1) = 2*P1(2:end-1);

% make the frequency bins
fftFreq = 50 * (0:(endPoint/2))/endPoint;

figure; plot(fftFreq,P1);
title('FFT'); 
set(gca, 'YScale', 'log', 'XScale', 'log');

%% Baseline DA: plot water, airpuff, stim, tone

timeRange = [-0.5,3];
eventRange = {'Water','Airpuff','Stim','Tone'};
animalRange = {'SL155','SL156','SL157','SL158','SL159'};
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
        plotTraces(combined.data{1},combined.timestamp,color=colorList{i},plotShuffled=true);
        xlabel('Time (s)'); ylabel([signalRange{s},' z-score']); ylim([-1,4]);
        plotEvent(eventRange{i},eventDuration(i),color=colorList{i});
        legend({['Shuffled (n=',num2str(size(combined.data{1},1)),')'],...
                [eventRange{i},' (n=',num2str(size(combined.data{1},1)),')']},...
                'Location','northeast');
    end
    saveFigures(gcf,'Summary_random_dLight',...
            strcat(resultspath),...
            saveFIG=true,savePDF=true);
end

%% Baseline GCaMP: plot water, airpuff, stim, tone

timeRange = [-0.5,3];
eventRange = {'Water','Airpuff','Stim','Tone'};
animalRange = {'SL157','SL158','SL159'};
taskRange = 'Random';
trialRange = 'All'; % range of trials in each session
totalTrialRange = 'All';
signalRange = {'LHb'};

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
                                    signalRange=signalRange);
        plotTraces(combined.data{1},combined.timestamp,color=colorList{i},plotShuffled=true);
        xlabel('Time (s)'); ylabel([signalRange{s},' z-score']); ylim([-1,4]);
        plotEvent(eventRange{i},eventDuration(i),color=colorList{i});
        legend({['Shuffled (n=',num2str(size(combined.data{1},1)),')'],...
                [eventRange{i},' (n=',num2str(size(combined.data{1},1)),')']},...
                'Location','northeast');
    end
    saveFigures(gcf,'Summary_random_GCaMP',...
            strcat(resultspath),...
            saveFIG=true,savePDF=true);
end

%% Baseline iGluSnFR: plot water, airpuff, stim, tone

timeRange = [-0.5,3];
eventRange = {'Water','Airpuff','Stim','Tone'};
animalRange = {'SL155','SL156'};
taskRange = 'Punish1';
trialRange = 'All'; % range of trials in each session
totalTrialRange = 'All';
signalRange = {'LHb'};

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
                                    signalRange=signalRange);
        plotTraces(combined.data{1},combined.timestamp,color=colorList{i},plotShuffled=true);
        xlabel('Time (s)'); ylabel([signalRange{s},' z-score']); ylim([-1,4]);
        plotEvent(eventRange{i},eventDuration(i),color=colorList{i});
        legend({['Shuffled (n=',num2str(size(combined.data{1},1)),')'],...
                [eventRange{i},' (n=',num2str(size(combined.data{1},1)),')']},...
                'Location','northeast');
    end
    % saveFigures(gcf,'Summary_random_iGluSnFR',...
    %         strcat(resultspath),...
    %         saveFIG=true,savePDF=true);
end

%% Plot overall to show animal learned

timeRange = [-0.5,3];
eventRange = {'Stim','Pair','Tone'};
animalRange = 'All'; %{'SL155','SL156','SL157','SL158'};
taskRange = {'Reward1','Punish1'};
totalTrialRange = 'All';
trialRange = 'All';
signalRange = 'NAc';

colorList = {bluePurpleRed(500,:),bluePurpleRed(300,:),bluePurpleRed(100,:)};
groupSizeList = [50;50;50];
nGroupsList = [15;15;15];
ylimList = [-1.2,3.5; -1.2,1];

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
                        groupby='sessions',startIdx=combined.options.startIdx,remaining='include');
        ylim(ylimList(task,:));
        plotEvent(eventRange{event},.5,color=colorList{event});
        xlabel('Time (s)'); ylabel([signalRange,' z-score']);
        legend(legendList,'Location','northeast');
    end
    % saveFigures(gcf,['Summary_pairing_',taskRange{task}],...
    %     strcat(resultspath),...
    %     saveFIG=true,savePDF=true);
end
% autoArrangeFigures

%% Plot lick trace

timeRange = [-0.5,3];
eventRange = {'Stim','Pair','Tone'};
animalRange = 'All';
taskRange = {'Reward1','Punish1'};
totalTrialRange = 'All';
trialRange = 'All';
signalRange = 'Lick';

colorList = {bluePurpleRed(500,:),bluePurpleRed(300,:),bluePurpleRed(100,:)};
groupSizeList = [30;30;30];
nGroupsList = [5;5;5];

for task = 1:length(taskRange)
    initializeFig(1,.5); tiledlayout('flow');
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
                        groupby='sessions',startIdx=combined.options.startIdx);
        plotEvent(eventRange{event},.5,color=colorList{event});
        xlabel('Time (s)'); ylabel('Licks/s'); ylim([0 Inf]);
        legend(legendList,'Location','northeast');
    end
    % saveFigures(gcf,['Summary_licking_',taskRange{task},'-',eventRange{event}],...
    %     strcat(resultspath),...
    %     saveFIG=true,savePDF=true);
end
% autoArrangeFigures

%% Plot grouped CS DA response (grouped across animal and 10 trials)

animalRange = 'All';
taskRange = {'Reward1','Punish1'};
conditionRange = [1,200];
signalRange = 'NAc';

% eventRange = {'Stim','Pair','Tone','Baseline'};
% colorList = {bluePurpleRed(500,:),bluePurpleRed(300,:),bluePurpleRed(100,:),[0.75,0.75,0.75]};

eventRange = {'Baseline','Stim','Pair'};
colorList = {[0.75,0.75,0.75],bluePurpleRed(500,:),bluePurpleRed(300,:)};
stage = 2; % Plot CS only
statsTypes = {'stageMax','stageMin','stageAvg','stageAvg'};
ylabelList = {'Max DA response during cue','Min DA response during cue'};

groupSize = 10; % numbers of trials to calculate average

combinedStats = getGroupedTrialStats(animals,statsTypes,...
                            eventRange=eventRange,...
                            animalRange=animalRange,...
                            taskRange=taskRange,...
                            totalTrialRange=conditionRange,...
                            signalRange=signalRange);

initializeFig(.7,.7); tiledlayout('flow');
results = plotGroupedTrialStats(combinedStats,ylabelList,groupSize=10,color=colorList,xlimIdx=2);

% saveFigures(gcf,['Summary_CSvsTrialsGrouped'],...
%         strcat(resultspath),...
%         saveFIG=true,savePDF=true);

%% Plot bar plot of DA slopes
initializeFig(.7,.7); tiledlayout('flow');
for task = 1:length(results.stats)
    nexttile;
    cur_stats = results.stats{task};
    for event = 1:length(eventRange)
        slopes = cur_stats{event}(:,1);

        % plot bar plots
        scatter(event,slopes,200,colorList{event},'filled'); hold on
        bar(event,mean(slopes),EdgeColor=colorList{event},FaceColor=colorList{event},...
            FaceAlpha=0.5,LineWidth=3);
        errorbar(event,mean(slopes),getSEM(slopes),Color=colorList{event},LineWidth=2);

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

%% Plot bar plot of DA final values
initializeFig(.7,.7); tiledlayout('flow');
for task = 1:length(results.stats)
    nexttile;
    cur_traces = results.traces{task};
    for event = 1:length(eventRange)
        finalData = cur_traces{event}(:,end-3:end);
        animalData = mean(finalData,2);

        % plot bar plots
        scatter(event,animalData,200,colorList{event},'filled'); hold on
        bar(event,mean(animalData),EdgeColor=colorList{event},FaceColor=colorList{event},...
            FaceAlpha=0.5,LineWidth=3);
        errorbar(event,mean(animalData),getSEM(animalData),Color=colorList{event},LineWidth=2);

        % Calculate significancesea
        if event < length(eventRange)
            for i = event+1:length(eventRange)
                plotStats(animalData,mean(cur_traces{i}(:,end-3:end),2),[event i],testType='kstest');
            end
        end
    end
    xticks(1:length(eventRange)); xticklabels(eventRange);
    ylabel('Learned CS response for DA');
end

%% Plot grouped anticipatory lick changes

eventRange = {'Stim','Pair'};
animalRange = 'All';
taskRange = {'Reward1','Punish1'};
signalRange = 'NAc';
colorList = {bluePurpleRed(500,:),bluePurpleRed(300,:),bluePurpleRed(100,:)};
ylabelList = {'Anticipatory licks','Anticipatory licks'};

combinedStats = getGroupedTrialStats(animals,'nAnticipatoryLicks',...
                            eventRange=eventRange,...
                            animalRange=animalRange,...
                            taskRange=taskRange,...
                            totalTrialRange=conditionRange,...
                            signalRange=signalRange);

initializeFig(.7,.7); tiledlayout('flow');
plotGroupedTrialStats(combinedStats,ylabelList,groupSize=10,color=colorList,...
                       eventRange=eventRange);

% saveFigures(gcf,['Summary_AnticipatoryLicksvsTrials_',taskRange{task},'-',eventRange{event}],...
%         strcat(resultspath),...
%         saveFIG=true,savePDF=true);

%% Plot anticipatory lick changes

eventRange = {'Stim','Pair'};
animalRange = 'All';
taskRange = {'Reward1','Punish1'};
signalRange = 'NAc';
conditionColors = {bluePurpleRed(1,:),[.213 .543 .324]};

% Get subtrial stats
stats = getGroupedTrialStats(animals,'nAnticipatoryLicks',...
                            eventRange=eventRange,...
                            animalRange=animalRange,...
                            taskRange=taskRange,...
                            totalTrialRange=conditionRange,...
                            signalRange=signalRange,...
                            concatSessions=false);

% Plot scatter plot and best fit line
for task = 1:length(taskRange)
    initializeFig(.7,.7); tiledlayout('flow');
    stats_combined = stats{task};
    for event = 1:length(eventRange)
        for animal = 1:length(animalList)
            nexttile;
            data = stats_combined{animal,event};
            lastTrial = 0;
            for row = 1:size(data,1)
                for col = 1:size(data,2)
                    plotData = data{row,col};
                    x = (1:size(plotData,1))+lastTrial;
                    y = plotData;
                    scatter(x,y,100,conditionColors{col},"filled",'MarkerFaceAlpha',0.5,HandleVisibility='off'); hold on
                    lastTrial = size(plotData,1)+lastTrial;
    
                    % Calc best fit line of trial vs CS response
                    p = polyfit(x,y',1);
                    plot(x,polyval(p,x),Color=conditionColors{col},lineWidth=5);
                end
            end
            xlabel('Trials'); ylabel([signalRange,' CS response']);
            title([taskRange{task}, ': ',animalList{animal},' -> ',eventRange{event}]);
        end
    end
    % saveFigures(gcf,['Summary_AnticipatoryLicksvsTrials_',taskRange{task},'-',eventRange{event},'_',signalRange],...
    %     strcat(resultspath),...
    %     saveFIG=true,savePDF=true);
end
% autoArrangeFigures

%% Plot CS DA response vs trials
eventRange = {'Stim','Pair'};
animalRange = 'All';
taskRange = {'Reward1','Punish1'};
conditionRange = [1,150];
signalRange = 'NAc';
conditionColors = {bluePurpleRed(1,:),[.213 .543 .324]};

stage = 2; % Plot CS only
statsType = 'stageAvg';

% Get subtrial stats
combinedStats = getGroupedTrialStats(animals,statsType,...
                            eventRange=eventRange,...
                            animalRange=animalRange,...
                            taskRange=taskRange,...
                            totalTrialRange=conditionRange,...
                            signalRange=signalRange,...
                            concatSessions=false);

for task = 1:length(taskRange)
    % Plot scatter plot and best fit line
    stats_combined = combinedStats.stats{task};
    initializeFig(.7,.7); tiledlayout('flow');
    for event = 1:length(eventRange)
        for animal = 1:length(animalList)
            nexttile;
            data = stats_combined{animal,event};
            lastTrial = 0;
            for row = 1:size(data,1)
                for col = 1:size(data,2)
                    plotData = data{row,col};
                    x = (1:size(plotData,1))+lastTrial;
                    y = plotData(:,stage);
                    scatter(x,y,100,conditionColors{col},"filled",'MarkerFaceAlpha',0.5,HandleVisibility='off'); hold on
                    lastTrial = size(plotData,1)+lastTrial;
    
                    % Calc best fit line of trial vs CS response
                    p = polyfit(x,y',1);
                    plot(x,polyval(p,x),Color=conditionColors{col},lineWidth=5);
                end
            end
            xlabel('Trials'); ylabel([signalRange,' CS response']);
            title([taskRange{task}, ': ',animalList{animal},' -> ',eventRange{event}]);
        end
    end
    % saveFigures(gcf,['Summary_CSvsTrials_',taskRange{task},'-',eventRange{event},'_',signalRange],...
    %     strcat(resultspath),...
    %     saveFIG=true,savePDF=true);
end
% autoArrangeFigures

%% 