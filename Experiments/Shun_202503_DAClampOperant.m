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
resultspath = osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Project clamping/Results');

% Building summary struct from selected sessions
answer = questdlg('Group sessions or load combined data?','Select load sources',...
                  'Group single sessions','Load combined data','Load sample data','Load combined data');

if strcmpi(answer,'Group single sessions')
    sessionList = uipickfiles('FilterSpec',osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Project clamping/Recordings'))';
    groupSessions = true;
    % Update resultspath
    dirsplit = strsplit(sessionList{1},filesep); projectName = dirsplit{end-1}; 
    resultspath = strcat(resultspath,filesep,projectName);
    % Create resultspath if necessary
    if isempty(dir(resultspath)); mkdir(resultspath); end

elseif strcmpi(answer,'Load combined data')
    fileList = uipickfiles('FilterSpec',osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Project clamping/Results'))';
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

baselineColor = [.7 .7 .7];
pidColor = [.232 .76 .58];

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
    cur_name = summary(i).name;

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
    elseif contains('Pair',cur_event,IgnoreCase=true)
        summary(i).event = 'Tone';
    end

    if strcmpi('PMT',cur_name)
        summary(i).name = 'dLight';
    end
end

% Remove some rows if needed
rowIdx = cellfun(@(x) contains(x,'Blue stim',IgnoreCase=true), {summary.event});
summary(rowIdx) = [];

rowIdx = cellfun(@(x) contains(x,'Stim',IgnoreCase=true), {summary.event});
summary(rowIdx) = [];

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

%% Define animal group

clampAnimals = {'SL287'};
unclampAnimals = setdiff(unique({animals.animal}),clampAnimals);

unclampColor = [0.7 0.7 0.7];
clampColor = [];

%% Test: Plot traces from summary/animals struct

event = 'Water';
animal = unclampAnimals;
task = 'Reward';
signal = 'dLight';
session = 'all';

initializeFig(0.5,0.5);
combined_unclamped = combineTraces(animals,timeRange=[-0.5,3],...
                            eventRange=event,...
                            animalRange=unclampAnimals,...
                            taskRange=task,...
                            totalTrialRange='All',...
                            trialRange='All',...
                            signalRange=signal,...
                            sessionRange=session);
plotTraces(combined_unclamped.data{1},combined_unclamped.timestamp,color=unclampColor);

combined_clamped = combineTraces(animals,timeRange=[-0.5,3],...
                            eventRange=event,...
                            animalRange=clampAnimals,...
                            taskRange=task,...
                            totalTrialRange='All',...
                            trialRange='All',...
                            signalRange=signal,...
                            sessionRange=session);
plotTraces(combined_clamped.data{1},combined_clamped.timestamp,color=clampColor);

plotEvent(event,0,color=bluePurpleRed(500,:))
xlabel('Time (s)'); ylabel('z-score');

%% dLight: Plot overall to show animal learned

timeRange = [-0.5,3];
eventRange = {'Stim','Pair','Tone'};
animalRange = 'All';
totalTrialRange = 'All';
trialRange = [1,150];
signalRange = 'dLight';

colorList = {bluePurpleRed(500,:),bluePurpleRed(300,:),bluePurpleRed(100,:)};
groupSizeList = [30;30;10];
nGroupsList = [15;15;15];

taskRange = {'Reward','Punish'}; % second part
ylimList = [-1,4; -1, 2.5]; % Second part

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
        ylim(ylimList(task,:));
        plotEvent(eventRange{event},.5,color=colorList{event});
        xlabel('Time (s)'); ylabel([signalRange,' z-score']);
        legend(legendList,'Location','northeast');
    end
    saveFigures(gcf,['Summary_dLight_',taskRange{task}],...
        strcat(resultspath),...
        saveFIG=true,savePDF=true);
end

%% Plot lick trace

timeRange = [-0.5,3];
eventRange = {'Tone'};
animalRange = 'All';
totalTrialRange = 'All';
trialRange = 'All';
signalRange = 'Lick';

colorList = {bluePurpleRed(500,:),bluePurpleRed(300,:),bluePurpleRed(100,:)};
groupSizeList = [30;30;30];
nGroupsList = [5;5;5];

% Second part
taskRange = {'Reward'};
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
    % saveFigures(gcf,['Summary_licking_',taskRange{task}],...
    %     strcat(resultspath),...
    %     saveFIG=true,savePDF=true);
end

%% Plot grouped CS DA response (grouped across animal and 10 trials)

taskRange = {'Reward','Punish'};
statsTypes = {'stageMax','stageMax'}; ylabelList = {'Max DA response during cue','Max DA response during cue'};
ylimit = [0,4.5; 0,4.5];

animalRange = 'All';
conditionRange = 'All';
%signalRange = 'dLight';

eventRange = {'Baseline','Pair','Tone','Stim'};
colorList = {[0.8,0.8,0.8],bluePurpleRed(300,:),bluePurpleRed(100,:),bluePurpleRed(500,:)};
% eventRange = {'Baseline','Stim','Pair'};
% colorList = {[0.8,0.8,0.8],bluePurpleRed(500,:),bluePurpleRed(300,:)};

stage = 2; % Plot CS only
groupSize = 10; % numbers of trials to calculate average
pairingStats = getGroupedTrialStats(animals,statsTypes,...
                            eventRange=eventRange,...
                            animalRange=animalRange,...
                            taskRange=taskRange,...
                            totalTrialRange=conditionRange,...
                            signalRange=signalRange);

initializeFig(.7,.7); tiledlayout('flow');
results_dLight = plotGroupedTrialStats(pairingStats,ylabelList,groupSize=10,color=colorList,xlimIdx=4,xlim=[0,170],ylim=ylimit);

% saveFigures(gcf,'Summary_CSvsTrialsGrouped_dLight',...
%         strcat(resultspath),...
%         saveFIG=true,savePDF=true);


%% Behavior: anticipatory licking vs trials

eventRange = 'Tone';
taskRange = 'Reward';
statsTypes = 'nAnticipatoryLicks';
conditionRange = 'All';
ylabels = 'Anticipatory licks';

tone_al_unclamped = getGroupedTrialStats(animals,statsTypes,...
                            eventRange=eventRange,...
                            animalRange=unclampAnimals,...
                            taskRange=taskRange,...
                            totalTrialRange=conditionRange,...
                            signalRange='Lick');
% tone_al_clamped = getGroupedTrialStats(animals,statsTypes,...
%                             eventRange=eventRange,...
%                             animalRange=clampAnimals,...
%                             taskRange=taskRange,...
%                             totalTrialRange=conditionRange,...
%                             signalRange='Lick');

results_tone_al_unclamped = plotGroupedTrialStats(tone_al_unclamped,ylabels,groupSize=10,...
                            color=baselineColor,...
                            plot=true,plotNextTile=false);
% results_tone_al_clamped = plotGroupedTrialStats(tone_al_clamped,ylabels,groupSize=10,...
%                             color=pidColor,...
%                             plot=true,plotNextTile=false);


%% Behavior: total licks vs trial

eventRange = 'Tone';
taskRange = 'Reward';
statsTypes = 'nLicks';
conditionRange = 'All';
ylabels = 'Total licks';

tone_al_unclamped = getGroupedTrialStats(animals,statsTypes,...
                            eventRange=eventRange,...
                            animalRange=unclampAnimals,...
                            taskRange=taskRange,...
                            totalTrialRange=conditionRange,...
                            signalRange='Lick');
% tone_al_clamped = getGroupedTrialStats(animals,statsTypes,...
%                             eventRange=eventRange,...
%                             animalRange=clampAnimals,...
%                             taskRange=taskRange,...
%                             totalTrialRange=conditionRange,...
%                             signalRange='Lick');

results_tone_al_unclamped = plotGroupedTrialStats(tone_al_unclamped,ylabels,groupSize=10,...
                            color=baselineColor,...
                            plot=true,plotNextTile=false);
% results_tone_al_clamped = plotGroupedTrialStats(tone_al_clamped,ylabels,groupSize=10,...
%                             color=pidColor,...
%                             plot=true,plotNextTile=false);

%% Behavior: CDF of success trials

%successData_clamped = cell(size(clampAnimals));
successData_unclamped = cell(size(unclampAnimals));
newLength = 100;

terminateAtMinTrials = false;
minTrials = inf;
maxTrials = 0;

% for a = 1:length(clampAnimals)
%     cur_animal = clampAnimals{a}; 
% 
%     % Extract anticipatory licking for the current animal
%     nAnticipatoryLicks = animals(a).trialInfo.trialTable.nAnticipatoryLicks;
% 
%     % Convert anticipatory licking to successTrials 
%     successTrials = nAnticipatoryLicks >= 3; 
% 
%     % Store success data
%     successData_clamped{a} = successTrials;
% end

for a = 1:length(unclampAnimals)
    % Find current animal row
    nAnticipatoryLicks = tone_al_unclamped.stats{1}{a};
    
    % Convert anticipatory licking to success trials
    successTrials = nAnticipatoryLicks >= 3; 
    
    % Store the processed success data
    successData_unclamped{a} = double(successTrials);
end

%cut at minimum number of trials or fill trials with 0 to max number of
%trials f
minTrials = 283;
maxTrials = 378;
unclampedData = [];
for a = 1:length(unclampAnimals)
    successTrials = successData_unclamped{a};
    % Option to truncate to minimum trials or extend to maximum trials
    if terminateAtMinTrials
        % Truncate to the minimum number of trials
        successTrials = successTrials(1:minTrials);
    else
        % Extend to the maximum number of trials by filling with 0s
        if length(successTrials) < maxTrials
            successTrials = [successTrials; false(maxTrials - length(successTrials),1)];
        end
    end
    unclampedData = [unclampedData, successTrials];
end

% Plot CDF (ask ChatGPT)
figure; hold on;

%plot
%ecdf(clampedData);
%hold on;
for a = 1:length(unclampAnimals)
    ecdf(find(unclampedData(:,a)));
    hold on;
end
