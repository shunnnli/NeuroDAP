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
loadNeuroDAP;
[~,~,~,~,~,~,bluePurpleRed] = loadColors;
clampColor = [.232 .76 .58];
unclampColor = [165, 209, 178]./255;

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

%% Optional: Create summary struct (only need to do this for initial loading)

if groupSessions    
    summary = concatAnalysis(sessionList,skipCamera=true);
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
    cur_session = summary(i).session;

    if strcmp('random',cur_session)
        summary(i).task = 'Random';
        if contains(cur_session,["unclamp","ctrl"],IgnoreCase=true)
            summary(i).task = 'Random-ctrl';
        else
            summary(i).task = 'Random-clamp';
        end
    end
end

% Change / add more details to task
keepRows = true(1, length(summary));
for i = 1:length(summary)

    cur_animal = string(summary(i).animal);
    cur_name   = string(summary(i).name);
    cur_task   = string(summary(i).task);
    cur_event  = string(summary(i).event);
    cur_date   = string(summary(i).date);

    skipGroup1 = any(strcmpi(cur_animal, "BiPOLES2")) && ...
                 any(strcmpi(cur_name, "NAc-right"));

    skipGroup2 = any(strcmpi(cur_animal, "M431")) && ...
                 any(strcmpi(cur_name, "NAc-right"));

    skipGroup3 = any(strcmpi(cur_animal, "M430")) && ...
                 any(strcmpi(cur_name, "NAc-left")) && ...
                 any(strcmpi(cur_date, "20260714"));

    if skipGroup1 || skipGroup2 || skipGroup3
        keepRows(i) = false;
    end
end
summary = summary(keepRows);

%% Optional: for ONOFF sessions only

% Change some names if needed
for i = 1:length(summary)
    cur_task = summary(i).task;
    cur_event = summary(i).event;
    cur_date = str2double(summary(i).date);
    cur_session = summary(i).session;

    % ONOFF type
    if cur_date == 20260418
        summary(i).task = 'Raw';
    elseif cur_date == 20260421 
        summary(i).task = 'withBump';
    elseif cur_date >= 20260423
        summary(i).task = 'Final';
    end
end

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

%% Optional: ONOFF only (remove artifact clamp commands)

for i = 1:size(animals,2)
    if contains('clamp',animals(i).name, IgnoreCase=true)
        data = animals(i).data;
        % data(:, [1:749, 1000:end]) = 0;
        
        % Remove trials where there's still clamp command after 6sec
    end
    
    animals(i).stageAmp = struct('data', getAmplitude(stageMax, stageMin));
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

%% Random: clamp on vs clamp off

timeRange = [-5,10];
taskRange = 'Final';
animalRange = 'All'; %{'SL431','SL432','BiPOLES2'};

trialRange = 'All'; % range of trials in each session
totalTrialRange = 'All';
signalRange = {'NAc-left','NAc-right'};
% signalRange = {'blueClamp','redClamp'};
trialConditions = 'trials.performing';

eventDuration = 5;

close all; 
initializeFig(.5,.5); tiledlayout(1,length(signalRange));
for s = 1:length(signalRange)
    nexttile;
    combined = combineTraces(animals,timeRange=timeRange,...
                                eventRange='Clamp',...
                                animalRange=animalRange,...
                                taskRange=taskRange,...
                                totalTrialRange=totalTrialRange,...
                                trialRange=trialRange,...
                                signalRange=signalRange{s},...
                                trialConditions=trialConditions);
    plotTraces(combined.data{1},combined.timestamp,color=clampColor,plotIndividual=true);
    xlabel('Time (s)'); ylabel([signalRange{s},' (\DeltaF/F)']);
    ylim([-0.3,1.4]);
    plotEvent('Clamp',eventDuration,color=clampColor);
    legend({['Clamp (n=',num2str(size(combined.data{1},1)),')']},...
            'Location','northeast');
end
% saveFigures(gcf,strcat('Summary_DA_ONOFF_',taskRange),...
%         strcat(resultspath),...
%         saveFIG=false,savePDF=true);

%% Calculate running var of clamp off vs clamp on
% use combined data above
% x axis is time, y axis is var in a window (1s)

initializeFig(.4,.5); tiledlayout(length(signalRange),4);
for s = 1:length(signalRange)
    combined = combineTraces(animals,timeRange=[0,10],...
                                eventRange='Clamp',...
                                animalRange=animalRange,...
                                taskRange=taskRange,...
                                totalTrialRange=totalTrialRange,...
                                trialRange=trialRange,...
                                signalRange=signalRange{s},...
                                trialConditions=trialConditions);
    windowSize = 3 * combined.options.finalFs;
    runningVar = movvar(combined.data{1},windowSize,0,2,'omitnan');

    nexttile((s-1)*4+1,[1 3]);
    plotTraces(runningVar,combined.timestamp,color=clampColor,plotIndividual=false);
    xlabel('Time (s)'); ylabel([signalRange{s},' var']);
    plotEvent('Clamp',eventDuration,color=clampColor);

    % Plot average variability of on vs off using plotScatterBar
    offIdx = combined.timestamp >= eventDuration;
    onIdx = combined.timestamp < eventDuration;
    varData = [mean(runningVar(:,onIdx),2,'omitnan'), ...
               mean(runningVar(:,offIdx),2,'omitnan')];

    nexttile((s-1)*4+4);
    plotScatterBar([1 2],varData,style='bar',color=[clampColor;unclampColor],...
                   plotScatter=false,connectPairs=false,dotSize=80);

    % plotStats(varData(:,1),varData(:,2),[1 2],testType='kstest');
    xticks([1 2]); xticklabels({'Off','On'});
    ylabel([signalRange{s},' var']);
end


%% Random: clamp vs unclamp
timeRange = [-0.5,3];
eventRange = {'Water','Airpuff','Tone'};
animalRange = 'All';

trialRange = 'All'; % range of trials in each session
totalTrialRange = 'All';
signalRange = {'NAc-left','NAc-right'};
trialConditions = 'trials.performing';

colorList = {bluePurpleRed(1,:),[.2,.2,.2],bluePurpleRed(500,:)};
eventDuration = [0,.2,.5];

close all; 
for i = 1:length(eventRange)
    initializeFig(.5,.5); tiledlayout(1,length(signalRange));
    for s = 1:length(signalRange)
        nexttile;
        combined = combineTraces(animals,timeRange=timeRange,...
                                    eventRange=eventRange{i},...
                                    animalRange=animalRange,...
                                    taskRange='Random-ctrl',...
                                    totalTrialRange=totalTrialRange,...
                                    trialRange=trialRange,...
                                    signalRange=signalRange{s},...
                                    trialConditions=trialConditions);
        plotTraces(combined.data{1},combined.timestamp,color=unclampColor);

        combined = combineTraces(animals,timeRange=timeRange,...
                                    eventRange=eventRange{i},...
                                    animalRange=animalRange,...
                                    taskRange='Random-clamp',...
                                    totalTrialRange=totalTrialRange,...
                                    trialRange=trialRange,...
                                    signalRange=signalRange{s},...
                                    trialConditions=trialConditions);
        plotTraces(combined.data{1},combined.timestamp,color=clampColor);
        xlabel('Time (s)'); ylabel([signalRange{s},' (\DeltaF/F)']); 
        ylim([-0.1,0.8]);
        plotEvent(eventRange{i},eventDuration(i),color=colorList{i});
        legend({[eventRange{i},' (n=',num2str(size(combined.data{1},1)),')']},...
                'Location','northeast');
    end
    saveFigures(gcf,strcat('Summary_random_',eventRange{i}),...
            strcat(resultspath),...
            saveFIG=false,savePDF=true);
end
