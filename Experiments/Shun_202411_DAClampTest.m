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
resultspath = osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Project clamping/');

% Building summary struct from selected sessions
answer = questdlg('Group sessions or load combined data?','Select load sources',...
                  'Group single sessions','Load combined data','Load sample data','Load combined data');

if strcmpi(answer,'Group single sessions')
    sessionList = uipickfiles('FilterSpec',osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Project clamping/Recordings'));
    groupSessions = true;
    % Update resultspath
    dirsplit = strsplit(sessionList{1},filesep); projectName = dirsplit{end-1}; 
    resultspath = strcat(resultspath,filesep,projectName);
    % Create resultspath if necessary
    if isempty(dir(resultspath)); mkdir(resultspath); end

elseif strcmpi(answer,'Load combined data')
    fileList = uipickfiles('FilterSpec',osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Project clamping/Results'));
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
        if isstring(summary(row).(stringColumnsLabels{i})) 
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

    if contains('Pair',cur_event,IgnoreCase=true)
        summary(i).event = 'Tone';
    end
end

% Remove some rows if needed
eventIdx = cellfun(@(x) contains(x,'dLight',IgnoreCase=true), {summary.name});
summary(eventIdx) = [];
eventIdx = cellfun(@(x) contains(x,'Stim only',IgnoreCase=true), {summary.event});
summary(eventIdx) = [];

%% Create animals struct

if isempty(dir(fullfile(resultspath,'animals*.mat'))) || groupSessions
    animals = getAnimalsStruct(summary);
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

%% Notes: trials to remove

%{
SL325
- D1: 1-135, did not start task
- D3: satiated in later trials, don't think its necessary
%}

clampAnimals = {'SL287'};
unclampAnimals = setdiff(unique({animals.animal}),clampAnimals);

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

%% Test: Plot traces from summary/animals struct

initializeFig(0.5,0.5);
combined = combineTraces(animals,timeRange=[-0.5,3],...
                            eventRange='Water',...
                            animalRange="All",...
                            taskRange='Reward',...
                            totalTrialRange='All',...
                            trialRange='All',...
                            signalRange='dLight',...
                            sessionRange='All');
plotTraces(combined.data{1},combined.timestamp,color=bluePurpleRed(1,:));
plotEvent('Water',0,color=bluePurpleRed(1,:))
xlabel('Time (s)'); ylabel('z-score');

%% dLight: clampped vs non clampped

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

%% Plot overall to show animal learned

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

%% Plot lick trace

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
    initializeFig(1,.5); tiledlayout('flow');
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


%% Plot grouped CS DA response (grouped across animal and 10 trials)
%{
% First part
% taskRange = {'Reward1 (paAIP2)','Punish1 (paAIP2)'};
% taskRange = {'Reward2 (Ctrl)','Punish2 (paAIP2)'};
% statsTypes = {'stageMax','stageMin'}; ylabelList = {'Max DA response during cue','Max DA response during cue'};
% statsTypes = {'stageAvg','stageAvg'}; ylabelList = {'Avg DA response during cue','Avg DA response during cue'};
% ylim = [0,3];

% Second part
taskRange = {'Punish (paAIP2)','Punish (paAIP2)'};
statsTypes = {'stageAvg','stageDelta'}; ylabelList = {'Avg DA response during cue','Delta DA response during cue'};
ylim = [-0.5,1.5; -1.5,2.5];

animalRange = 'All';
conditionRange = 'All';
signalRange = 'dLight';

eventRange = {'Baseline','Pair','Tone','Stim'};
colorList = {[0.8,0.8,0.8],bluePurpleRed(300,:),bluePurpleRed(100,:),bluePurpleRed(500,:)};
% eventRange = {'Baseline','Stim','Pair'};
% colorList = {[0.8,0.8,0.8],bluePurpleRed(500,:),bluePurpleRed(300,:)};

stage = 2; % Plot CS only
groupSize = 10; % numbers of trials to calculate average
combinedStats = getGroupedTrialStats(animals,statsTypes,...
                            eventRange=eventRange,...
                            animalRange=animalRange,...
                            taskRange=taskRange,...
                            totalTrialRange=conditionRange,...
                            signalRange=signalRange);

initializeFig(.7,.7); tiledlayout('flow');
results = plotGroupedTrialStats(combinedStats,ylabelList,groupSize=10,color=colorList,xlimIdx=4,ylim=ylim);
% 
% saveFigures(gcf,'Summary_CSvsTrialsGrouped_Punish_Avg',...
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

% saveFigures(gcf,'Summary_CSvsTrialsGrouped_slope',...
%         strcat(resultspath),...
%         saveFIG=true,savePDF=true);

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

        % Calculate significance
        if event < length(eventRange)
            for i = event+1:length(eventRange)
                plotStats(animalData,mean(cur_traces{i}(:,end-3:end),2),[event i],testType='kstest');
            end
        end
    end
    xticks(1:length(eventRange)); xticklabels(eventRange);
    ylabel('Learned CS response for DA');
end

% saveFigures(gcf,'Summary_CSvsTrialsGrouped_finalVal',...
%         strcat(resultspath),...
%         saveFIG=true,savePDF=true);
%}

%% Plot grouped anticipatory lick changes

eventRange = {'Tone'};
taskRange = {'Reward'};
signalRange = 'dLight';
ylabelList = {'Anticipatory licks'};

combinedStats_ctrl = getGroupedTrialStats(animals,'nAnticipatoryLicks',...
                            eventRange=eventRange,...
                            animalRange=unclampAnimals,...
                            taskRange=taskRange,...
                            totalTrialRange='All',...
                            signalRange=signalRange);

combinedStats_pid = getGroupedTrialStats(animals,'nAnticipatoryLicks',...
                            eventRange=eventRange,...
                            animalRange=clampAnimals,...
                            taskRange=taskRange,...
                            totalTrialRange='All',...
                            signalRange=signalRange);

initializeFig(.5,.5); tiledlayout('flow');
plotGroupedTrialStats(combinedStats_pid,ylabelList,groupSize=10,color={bluePurpleRed(100,:)},...
                       eventRange=eventRange,xlimIdx=1,ylim=[0,7]);

plotGroupedTrialStats(combinedStats_ctrl,ylabelList,groupSize=10,color={[.7 .7 .7]},...
                       eventRange=eventRange,xlimIdx=1,ylim=[0,7]);

% saveFigures(gcf,'Summary_AnticipatoryLicksvsTrials',...
%         strcat(resultspath),...
%         saveFIG=true,savePDF=true);

%% Plot anticipatory lick changes

eventRange = {'Tone'};
taskRange = {'Reward'};
signalRange = 'dLight';
conditionColors = {bluePurpleRed(1,:),[.213 .543 .324]};

% Get subtrial stats
stats_ctrl = getGroupedTrialStats(animals,'nAnticipatoryLicks',...
                            eventRange=eventRange,...
                            animalRange=unclampAnimals,...
                            taskRange=taskRange,...
                            totalTrialRange='All',...
                            signalRange=signalRange,...
                            concatSessions=false);

stats_pid = getGroupedTrialStats(animals,'nAnticipatoryLicks',...
                            eventRange=eventRange,...
                            animalRange=clampAnimals,...
                            taskRange=taskRange,...
                            totalTrialRange='All',...
                            signalRange=signalRange,...
                            concatSessions=false);

% Plot scatter plot and best fit line
initializeFig(.7,.7); tiledlayout('flow');

stats = stats_pid.stats{1};
for event = 1:length(eventRange)
    for animal = 1:length(clampAnimals)
        nexttile;
        data = stats{animal,event};
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
        title([taskRange{1}, ': ',animalList{animal},' -> ',eventRange{event}]);
    end
end

stats = stats_ctrl.stats{1};
for event = 1:length(eventRange)
    for animal = 1:length(unclampAnimals)
        nexttile;
        data = stats{animal,event};
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
        title([taskRange{1}, ': ',animalList{animal},' -> ',eventRange{event}]);
    end
end
% saveFigures(gcf,['Summary_AnticipatoryLicksvsTrials_',taskRange{task},'-',eventRange{event},'_',signalRange],...
%     strcat(resultspath),...
%     saveFIG=true,savePDF=true);

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