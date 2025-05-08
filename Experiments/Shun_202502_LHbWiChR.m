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
addpath(genpath(osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Analysis/NeuroDAP/Methods')));
[~,~,~,~,~,~,bluePurpleRed] = loadColors;

% Define result directory
resultspath = osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Project valence/Results');

% Building summary struct from selected sessions
answer = questdlg('Group sessions or load combined data?','Select load sources',...
                  'Group single sessions','Load combined data','Load sample data','Load combined data');

if strcmpi(answer,'Group single sessions')
    sessionList = uipickfiles('FilterSpec',osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Project valence/Recordings'));
    groupSessions = true;
    % Update resultspath
    dirsplit = strsplit(sessionList{1},filesep); projectName = dirsplit{end-1}; 
    resultspath = strcat(resultspath,filesep,projectName);
    % Create resultspath if necessary
    if isempty(dir(resultspath)); mkdir(resultspath); end

elseif strcmpi(answer,'Load combined data')
    fileList = uipickfiles('FilterSpec',osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Project valence/Results'));
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
    resultspath = strcat(osPathSwith('/Volumes/Neurobio/MICROSCOPE/Shun/Analysis/NeuroDAP/Tutorials/Sample data/Results'),filesep,projectName);
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
        if isstring(summary(row).(stringColumnsLabels{1}))
            summary(row).(stringColumnsLabels{i}) = convertStringsToChars(summary(row).(stringColumnsLabels{i}));
        end
    end
end
disp('Finished: summary struct and trialtables loaded');

%% Optional: Make changes to summary for further analysis (first reward & punish sessions)

% Change some names if needed
for i = 1:length(summary)
    cur_task = summary(i).task;
    cur_event = summary(i).event;
    cur_date = str2double(summary(i).date);

    if strcmp('random',cur_task)
        summary(i).task = 'Random';
    elseif contains('reward pairing',cur_task,IgnoreCase=true)
        summary(i).task = 'Reward';
    elseif contains('punish pairing',cur_task,IgnoreCase=true)
        summary(i).task = 'Punish';
    end

    if contains('Stim only',cur_event,IgnoreCase=true)
        summary(i).event = 'Stim';
    elseif contains('Tone only',cur_event,IgnoreCase=true)
        summary(i).event = 'Tone';
    elseif contains('Blue stim',cur_event,IgnoreCase=true)
        summary(i).event = 'WiChR';
    end
end

%% Optional: remove not performing period

modifyRows = 623:652;
lastTrial = 80;

for i = modifyRows
    removeIdx = summary(i).trialInfo.trialTable.TrialNumber > lastTrial;
    summary(i).trialInfo.trialTable.performing(removeIdx) = false(sum(removeIdx),1);
end
disp(['Finished: ', num2str(modifyRows)]);

%% Select top 10 stim trials

modifyRows = 658:662;
topK = 10;

stimTrace = summary(modifyRows(1)).data(:,751:801);
stimTrace_smoothed = movmean(stimTrace,5,2);
if size(stimTrace,1) > 20
    stimAmp = max(stimTrace,[],2);
    [~, top20Indices] = maxk(stimAmp, topK);

    performing = false(length(stimAmp),1);
    performing(top20Indices) = true;

    for i = modifyRows
        nTrials = size(summary(i).trialInfo.trialTable.performing,1);
        performing = false(nTrials,1);
        performing(top20Indices) = true;
        summary(i).trialInfo.trialTable.performing = performing;
    end
else
    for i = modifyRows
        nTrials = size(summary(i).trialInfo.trialTable.performing,1);
        summary(i).trialInfo.trialTable.performing = true(nTrials,1);
    end
end
disp(['Finished: ', num2str(modifyRows)]);

%% Create animals struct

if isempty(dir(fullfile(resultspath,'animals*.mat'))) || groupSessions
    animals = getAnimalsStruct(summary);
end

% Add stageAmp
for i = 1:size(animals,2)
    stageMax = animals(i).stageMax.data;
    stageMin = animals(i).stageMin.data;
    
    animals(i).stageAmp = struct('data', getAmplitude(stageMax, stageMin));
end

%% Change window to 1s post stim

for i = 1:length(animals)
    if strcmpi(animals(i).event,'WiChR')
        if strcmpi(animals(i).name,'Lick')
            continue
        else
            stimTrace = animals(i).data;
            stats = analyzeStages(stimTrace,[-2,0;0:1],finalFs=50);

            % Store results
            animals(i).stageAvg.data = stats.stageAvg.data;
            animals(i).stageMax.data = stats.stageMax.data;
            animals(i).stageMin.data = stats.stageMin.data;
            animals(i).stageArea.data = stats.stageArea.data;
            animals(i).stageAmp.data = stats.stageMax.data - stats.stageMin.data;
        end
    end
end 
disp('Finished');

%% (Optional) Calculate periMax and periMin response

finalTrialsWindow = 25; % use last 50 trials to average
stageEndWindow = 120; % in ms, include these in previous stages
peakWindowWidth = 25; % in ms

for i = 1:size(animals,2)
    % Get traces for each animal and for each signal and event
    traces = animals(i).data;
    finalFs = animals(i).options.finalFs;
    eventSample = find(animals(i).timestamp==0);
    
    % Find the final session and plot the average trace for each animal
    fianlWindowStart = max(size(traces,1)-finalTrialsWindow-1,1);
    finalTrialsAvg = mean(traces(fianlWindowStart:end,:),1);
    
    % Find stage windows
    if ~isfield(animals(i).options,'stageTime')
        if strcmpi(animals(i).task,'random'); stageTime = [-2,0;0,2];
        else; stageTime = [-2,0;0,0.5;0.5,5]; end
    else
        stageTime = animals(i).options.stageTime;
    end
    stageWindow = stageTime * finalFs + eventSample;
    stageWindow(:,2) = round(stageWindow(:,2) + stageEndWindow/1000 * finalFs);
    periDelta = nan(size(traces,1),size(stageWindow,1));
    
    % Loop for each stage
    for stageIdx = 1:size(stageWindow,1)
        % Get current stage window
        curStageWindow = stageWindow(stageIdx,1):stageWindow(stageIdx,2);
        stageTraces = traces(:,curStageWindow);

        % Find the max and min window (+- 50ms)
        peakWindowWidthSamples = max(round(peakWindowWidth/1000 * finalFs),1);
        [~,maxIdx] = max(finalTrialsAvg(curStageWindow)); 
        maxWindowStart = max(1,maxIdx-peakWindowWidthSamples);
        maxWindowEnd = min(maxIdx+peakWindowWidthSamples,size(stageTraces,2));
        [~,minIdx] = min(finalTrialsAvg(curStageWindow));
        minWindowStart = max(1,minIdx-peakWindowWidthSamples);
        minWindowEnd = min(minIdx+peakWindowWidthSamples,size(stageTraces,2));
        
        % Create periMax and perMin calculation for each trace
        periMax = mean(stageTraces(:,maxWindowStart:maxWindowEnd),2);
        periMin = mean(stageTraces(:,minWindowStart:minWindowEnd),2);
        
        absDelta = abs(periMax)-abs(periMin);
        pickMaxIdx = absDelta > 0;
        pickMinIdx = absDelta < 0;
        pickAvgIdx = absDelta == 0;
        periDelta(pickMaxIdx,stageIdx) = periMax(pickMaxIdx);
        periDelta(pickMinIdx,stageIdx) = periMin(pickMinIdx);
        periDelta(pickAvgIdx,stageIdx) = periMax(pickAvgIdx) + periMin(pickAvgIdx);
        % periDelta(:,stageIdx) = periMax + periMin;

        % Test plotting
        % if stageIdx == 2 && any(strcmpi(animals(i).event,{'Stim','Pair'}))
        %     initializeFig(0.3,0.7); tiledlayout(4,1);
        %     nexttile; 
        %     plot(finalTrialsAvg); hold on;
        %     xline(curStageWindow(1),'Color','r'); hold on;
        %     xline(curStageWindow(end),'Color','r'); hold on;
        %     nexttile; 
        %     plot(finalTrialsAvg(curStageWindow)); hold on;
        %     xline(maxWindowStart,'Color','magenta'); hold on; 
        %     xline(maxWindowEnd,'Color','magenta'); hold on;
        %     xline(minWindowStart,'Color','k'); hold on;
        %     xline(minWindowEnd,'Color','k'); hold on;
        %     nexttile; 
        %     plot(periMax,'Color','magenta'); hold on
        %     plot(periMin,'Color','k'); hold on
        %     % plot(periDelta(:,2),'Color','r');
        %     nexttile;
        %     plot(periDelta(:,2),s'Color','r');
        %     disp(animals(i).name);
        %     disp(animals(i).event);
        %     disp('-------');
        %    close all;
        % end
    end

    % Save to animals
    animals(i).stageDelta.data = periDelta;
    animals(i).options.statsType{end+1} = 'stageDelta';
end

%% Save animals struct

prompt = 'Enter database notes (animals_20230326_notes.mat):';
dlgtitle = 'Save animals struct'; fieldsize = [1 45]; definput = {''};
answer = inputdlg(prompt,dlgtitle,fieldsize,definput);
today = char(datetime('today','Format','yyyyMMdd'));
filename = strcat('animals_',today,'_',answer{1});

% Save animals.mat
if ~isempty(answer)
    disp(['Ongoing: saving animals.mat (',char(datetime('now','Format','HH:mm:ss')),')']);
    save(strcat(resultspath,filesep,filename),'animals','trialTables','sessionList','-v7.3');
    disp(['Finished: saved animals.mat (',char(datetime('now','Format','HH:mm:ss')),')']);
end

%% Test: Plot traces from summary/animals struct

event = 'Stim';
animal = 'all';
task = 'Reward';
signal = 'dLight';
session = 'all';

initializeFig(0.5,0.5);
combined = combineTraces(animals,timeRange=[-0.5,3],...
                            eventRange=event,...
                            animalRange=animal,...
                            taskRange=task,...
                            totalTrialRange='All',...
                            trialRange='All',...
                            signalRange=signal,...
                            sessionRange=session);

legendList = plotGroupTraces(combined.data{1},combined.timestamp,bluePurpleRed,...
                        groupby='sessions',startIdx=combined.options.startIdx,remaining='include');

plotEvent('Stim',0,color=bluePurpleRed(500,:))
xlabel('Time (s)'); ylabel('z-score');

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

%% Plot lick trace

timeRange = [-0.5,3];
eventRange = {'Stim','Pair','Tone'};
animalRange = 'All';
totalTrialRange = 'All';
trialRange = 'All';
signalRange = 'Lick';

colorList = {bluePurpleRed(500,:),bluePurpleRed(300,:),bluePurpleRed(100,:)};
groupSizeList = [30;30;30];
nGroupsList = [5;5;5];

taskRange = {'Reward','Punish'};
ylimList = [0,7; 0,3.5];

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
                        groupby='trials',startIdx=combined.options.startIdx);
        ylim(ylimList(task,:));
        plotEvent(eventRange{event},.5,color=colorList{event});
        xlabel('Time (s)'); ylabel('Licks/s'); ylim([0 Inf]);
        legend(legendList,'Location','northeast');
    end
    saveFigures(gcf,['Summary_licking_',taskRange{task}],...
        strcat(resultspath),...
        saveFIG=true,savePDF=true);
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

%% Plot bar plot of DA final-first values for each event
% calculate the difference between the avg(final 30 trial) and avg(first 30
% trial) for each animal
% both mean and d' seems to yield same result (mean is more obvious)

trialWindow = 10; 

initializeFig(.7,.7); tiledlayout('flow');
for task = 1:length(results.stats)
    nexttile;
    task_stats = combinedStats.stats{task};
    for event = 1:length(eventRange)  
        eventData = task_stats(:,event);
        startData = cell2mat(cellfun(@(m) m(1:trialWindow,2)', eventData, UniformOutput=false));
        finalData = cell2mat(cellfun(@(m) m(end-trialWindow+1:end,2)', eventData, UniformOutput=false));

        % Option 1: mean
        diffData = mean(finalData,2) - mean(startData,2);
        % Option 2: d-prime
        % meanDiff = mean(finalData,2) - mean(startData,2);
        % varSum = var(finalData,0,2).^2 + var(startData,0,2).^2;
        % diffData = meanDiff ./ sqrt(varSum ./ 2);

        plotScatterBar(event,diffData(:),color=colorList{event},...
                       style='bar',dotSize=200,LineWidth=2);

        % Calculate significance
        if event < length(eventRange)
            for i = event+1:length(eventRange)
                otherStartData = cell2mat(cellfun(@(m) m(1:trialWindow,2)', task_stats(:,i), UniformOutput=false));
                otherFinalData = cell2mat(cellfun(@(m) m(end-trialWindow+1:end,2)', task_stats(:,i), UniformOutput=false));
                otherDiffData = mean(otherFinalData,2) - mean(otherStartData,2);
                plotStats(diffData(:),otherDiffData(:),[event i],testType='kstest');
            end
        end
    end
    xticks(1:length(eventRange)); xticklabels(eventRange);
    ylabel('Final - first (d-prime)');
end

% saveFigures(gcf,'Summary_CSvsTrialsGrouped_Pairing1_final-first_bar',...
%         strcat(resultspath),...
%         saveFIG=true,savePDF=true);


%% Plot stim amp for each animal

taskRange = {'Reward','Punish'};
statsTypes = {'stageAmp','stageAmp'}; ylabelList = {'Amp DA response during cue','Amp DA response during cue'};

animalRange = 'SL355';
conditionRange = 'All';
signalRange = 'dLight';
trialConditions = 'trials.performing == 1';

eventRange = {'Stim'};
colorList = {[.213 .543 .324]};

pairingStats_DA = getGroupedTrialStats(animals,statsTypes,...
                            eventRange=eventRange,...
                            animalRange=animalRange,...
                            taskRange=taskRange,...
                            totalTrialRange=conditionRange,...
                            signalRange=signalRange,...
                            trialConditions=trialConditions);

initializeFig(.7,.7); tiledlayout('flow');
results = plotGroupedTrialStats(pairingStats_DA,ylabelList,groupSize=10,color=colorList,plotIndividual=true);

% saveFigures(gcf,'Summary_CSvsTrialsGrouped_dLight',...
%         strcat(resultspath),...
%         saveFIG=true,savePDF=true);