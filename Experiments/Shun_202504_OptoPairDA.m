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

%%
% Define result directory
resultspath = osPathSwitch('/Volumes/MICROSCOPE/Shun/Project valence/Results');

% Building summary struct from selected sessions
answer = questdlg('Group sessions or load combined data?','Select load sources',...
                  'Group single sessions','Load combined data','Load sample data','Load combined data');

if strcmpi(answer,'Group single sessions')
    sessionList = uipickfiles('FilterSpec',osPathSwitch('/Volumes/MICROSCOPE/Shun/Project valence/Recordings'))';
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

%% Optional: Make changes to summary for further analysis (first reward & punish sessions)

% Change some names if needed
for i = 1:length(summary)
    cur_task = summary(i).task;
    cur_event = summary(i).event;
    cur_date = str2double(summary(i).date);
    cur_name = summary(i).name;

    if contains('random',cur_task,IgnoreCase=true)
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
    elseif strcmpi('Blue stim',cur_event)
        summary(i).event = 'Shutter';
    end

    if contains('NAc-green',cur_name,IgnoreCase=true)
        summary(i).name = 'dLight';
    end
end

% Remove some rows if needed
% eventIdx = cellfun(@(x) contains(x,'Blue stim',IgnoreCase=true), {summary.event});
% summary(eventIdx) = [];
% 
% eventIdx = cellfun(@(x) contains(x,'iGluSnFR',IgnoreCase=true), {summary.name});
% summary(eventIdx) = [];

%% Create animals struct

if isempty(dir(fullfile(resultspath,'animals*.mat'))) || groupSessions
    animals = getAnimalsStruct(summary);
end

% Add stageAmp
for i = 1:size(animals,2)
    stageMax = animals(i).stageMax.data;
    stageMin = animals(i).stageMin.data;
    
    animals(i).stageAmp = struct('data', getAmplitude(stageMax, stageMin));

    % if strcmpi(animals(i).name,'NAc')
    %     animals(i).name = 'dLight';
    % end
end

%% (Optional) Apply trialTable of dLight to corresponding Lick

for i = 1:size(animals,2)
    cur_system = animals(i).system;
    cur_event = animals(i).event;

    if contains(cur_event,{'Stim','Pair'}) && strcmpi(cur_system,'Lick')
        disp([cur_event,' ',cur_system]);
        disp([animals(i+1).event,' ',animals(i+1).system]);
        performing = animals(i+1).trialInfo.trialTable.performing;
        animals(i).trialInfo.trialTable.performing(1:length(performing)) = performing;
    end
end
disp('Finished');

%% (Optional) Change baseline to baseline of stim trials


%% (Optional) Isolate random sessions
% Load combined_random.mat
fileList = uipickfiles('FilterSpec',osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Project valence/Results'),...
                       'Prompt','Select combined_random.mat')';
dirsplit = strsplit(fileList{1},filesep);
disp(['Ongoing: loading ',dirsplit{end}]);
load(fileList{1});
disp(['Finished: loaded ',dirsplit{end}]);

random = animals(strcmpi({animals.task},'Random'));
combined_random = [combined_random, random];

%% (Optional) remove the last stage for SL107-115

for i = 1:length(combined_random)
    cur_animal = combined_random(i).animal;
    cur_statSample = size(combined_random(i).stageAvg.data,2);

    if cur_statSample == 3
        disp(cur_animal);
        combined_random(i).stageAvg.data = combined_random(i).stageAvg.data(:,1:2);
        combined_random(i).stageAmp.data = combined_random(i).stageAmp.data(:,1:2);
        combined_random(i).stageMax.data = combined_random(i).stageMax.data(:,1:2);
        combined_random(i).stageMin.data = combined_random(i).stageMin.data(:,1:2);
    end
end

%% (Optional) Save new combined_random.mat
prompt = 'Enter database notes (combined_random_20230326_notes.mat):';
dlgtitle = 'Save combined_random struct'; fieldsize = [1 45]; definput = {''};
answer = inputdlg(prompt,dlgtitle,fieldsize,definput);
today = char(datetime('today','Format','yyyyMMdd'));
if ~isfolder(strcat(resultspath,filesep,today)); mkdir(strcat(resultspath,filesep,today)); end
if isempty(answer{1}); filename = strcat('combined_random_',today);
else; filename = strcat('combined_random_',today,'_',answer{1}); end

disp(['Ongoing: saving combined_random.mat (',char(datetime('now','Format','HH:mm:ss')),')']);
save(strcat(resultspath,filesep,today,filesep,filename),'combined_random','-v7.3');
disp(['Finished: saved combined_random.mat (',char(datetime('now','Format','HH:mm:ss')),')']);

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

initializeFig(0.5,0.5);
combined = combineTraces(animals,timeRange=[-0.5,3],...
                            eventRange='Airpuff',...
                            animalRange="All",...
                            taskRange='Random',...
                            totalTrialRange='All',...
                            trialRange='All',...
                            signalRange='dLight',...
                            sessionRange='All');
plotTraces(combined.data{1},combined.timestamp,color=bluePurpleRed(1,:));
plotEvent('Water',0,color=bluePurpleRed(1,:))
xlabel('Time (s)'); ylabel('z-score');



%% Plot overall to show animal learned

timeRange = [-0.5,3];
eventRange = {'Stim','Pair','Tone'};
animalRange = 'all';
totalTrialRange = [1 200];
trialRange = 'All';
signalRange = 'dLight';
trialConditions = 'trials.performing';

colorList = {bluePurpleRed(500,:),bluePurpleRed(300,:),bluePurpleRed(100,:)};
groupSizeList = [50;50;10];
nGroupsList = [15;15;15];

taskRange = {'Reward','Punish'};
ylimList = [-1,3; -1.5,1];


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
                                    signalRange=signalRange,...
                                    trialConditions=trialConditions);
        legendList = plotGroupTraces(combined.data{1},combined.timestamp,bluePurpleRed,...
                        groupSize=groupSizeList(event),nGroups=nGroupsList(event),...
                        groupby='trials',startIdx=combined.options.startIdx,remaining='exclude');
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
eventRange = {'Stim','Pair','Tone'};
animalRange = 'All';
totalTrialRange = 'All';
trialRange = 'All';
signalRange = 'Lick';

colorList = {bluePurpleRed(500,:),bluePurpleRed(300,:),bluePurpleRed(100,:)};
groupSizeList = [30;30;30];
nGroupsList = [5;5;5];

taskRange = {'Punish'};
ylimList = [-1.5,1; -1.5,1; -1.5,1; -1.25,1.75]; 

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

saveFigures(gcf,'Summary_CSvsTrialsGrouped_Amp',...
        strcat(resultspath),...
        saveFIG=true,savePDF=true);

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

saveFigures(gcf,'Summary_CSvsTrialsGrouped_slopeBar',...
        strcat(resultspath),...
        saveFIG=true,savePDF=true);

%% Plot bar plot of DA final-first values for each event
% calculate the difference between the avg(final 30 trial) and avg(first 30
% trial) for each animal
% both mean and d' seems to yield same result (mean is more obvious)

trialWindow = 30; 

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
    ylabel('Final - first');
end

saveFigures(gcf,'Summary_CSvsTrialsGrouped_Pairing1_final-first_bar',...
        strcat(resultspath),...
        saveFIG=true,savePDF=true);

%% Plot bar plot of DA final values
% not as obvious but significant

trialWindow = 30; 

initializeFig(.7,.7); tiledlayout('flow');
for task = 1:length(results.stats)
    nexttile;
    task_stats = combinedStats.stats{task};
    for event = 1:length(eventRange)  
        eventData = task_stats(:,event);
        finalData = cell2mat(cellfun(@(m) m(end-trialWindow+1:end,2)', eventData, UniformOutput=false));
        plotScatterBar(event,finalData(:),color=colorList{event},...
                       style='box',dotSize=200,LineWidth=2);

        % Calculate significance
        if event < length(eventRange)
            for i = event+1:length(eventRange)
                otherFinalData = cell2mat(cellfun(@(m) m(end-trialWindow+1:end,2)', task_stats(:,i), UniformOutput=false));
                plotStats(finalData(:),otherFinalData(:),[event i],testType='kstest');
            end
        end
    end
    xticks(1:length(eventRange)); xticklabels(eventRange);
    ylabel('Learned CS response for DA');
end

% saveFigures(gcf,'Summary_CSvsTrialsGrouped_finalValBox',...
%         strcat(resultspath),...
%         saveFIG=true,savePDF=true);

%%










%% Plot grouped anticipatory lick changes

taskRange = {'Reward','Punish'};
ylabelList = {'Anticipatory licks','Anticipatory licks'};
ylimList = [1.5,4.5; 0,4];

animalRange = 'all';
% animalRange = {'SL043','SL044','SL046','SL060','SL062','SL064','SL066','SL068'};%'All';
conditionRange = 'All';
signalRange = 'Lick';
trialConditions = 'trials.performing';

eventRange = {'Baseline','Stim'};
colorList = {bluePurpleRed(500,:),bluePurpleRed(300,:)};

groupSize = 10; % numbers of trials to calculate average
combinedStats = getGroupedTrialStats(animals,'nAnticipatoryLicks',...
                            eventRange=eventRange,...
                            animalRange=animalRange,...
                            taskRange=taskRange,...
                            totalTrialRange=conditionRange,...
                            signalRange=signalRange, ...
                            trialConditions=trialConditions);

initializeFig(.7,.7); tiledlayout('flow');
results = plotGroupedTrialStats(combinedStats,ylabelList,groupSize=10,color=colorList,ylim=ylimList,...
                                xlimIdx=1,xlim=[1,150],plotIndividual=false);

% saveFigures(gcf,'Summary_CSvsTrialsGrouped_Lick',...
%         strcat(resultspath),...
%         saveFIG=true,savePDF=true);

%% Plot bar plot of anticipatory lick slopes
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

% saveFigures(gcf,'Summary_anticpatoryLicks_slopeBar',...
%         strcat(resultspath),...
%         saveFIG=true,savePDF=true);

%% (Baseline) Test plotting

% animalRange = {'SL043' 'SL044' 'SL045'	'SL046'	...
%     'SL060'	'SL061'	'SL062'	'SL063'	'SL064'	'SL065'	'SL066'	'SL067'	'SL068'	'SL069'	...
%     'SL107'	'SL108'	'SL109'	'SL110'	'SL112'	'SL113'	'SL114'	'SL115'...
%     'SL316'	'SL317'	'SL318'	'SL319'	'SL320'	'SL321'	'SL322'	'SL323'};

eventRange = 'water';
animalRange = 'SL062';

timeRange = [-0.5,3];
taskRange = 'Random';
trialRange = 'All'; % range of trials in each session
totalTrialRange = 'All';
signalRange = 'dLight';
trialConditions = 'trials.performing';

close all; initializeFig(.5,.5);
combined = combineTraces(combined_random,timeRange=timeRange,...
                            eventRange=eventRange,...
                            animalRange=animalRange,...
                            taskRange=taskRange,...
                            totalTrialRange=totalTrialRange,...
                            trialRange=trialRange,...
                            signalRange=signalRange,...
                            trialConditions=trialConditions);
legendList = plotGroupTraces(combined.data{1},combined.timestamp,bluePurpleRed,...
                groupSize=20,nGroups=50,...
                groupby='sessions',startIdx=combined.options.startIdx,remaining='include');
xlabel('Time (s)'); ylabel([signalRange,' z-score']); 
legend(legendList,'Location','northeast');
% saveFigures(gcf,'Summary_random_dLight',...
%         strcat(resultspath),...
%         saveFIG=true,savePDF=true);

%% (Baseline) Remove some learning trials

targetIdx = 948:953;
for i = targetIdx
    combined_random(i).trialInfo.trialTable.performing(397:455) = 1; 
end
disp('Finished');

%% Baseline dLight: plot water, airpuff, stim, tone

timeRange = [-0.5,3];
eventRange = {'Water','Airpuff','Stim'};
animalRange = 'All';

taskRange = 'Random';
trialRange = 'All'; % range of trials in each session
totalTrialRange = 'All';
signalRange = {'dLight'};
trialConditions = 'trials.performing';

colorList = {bluePurpleRed(1,:),[.2,.2,.2],bluePurpleRed(500,:)};
eventDuration = [0,.1,.5,.5];

for s = 1:length(signalRange)
    close all; initializeFig(.5,.5);
    for i = 1:length(eventRange)
        combined = combineTraces(combined_random,timeRange=timeRange,...
                                    eventRange=eventRange{i},...
                                    animalRange=animalRange,...
                                    taskRange=taskRange,...
                                    totalTrialRange=totalTrialRange,...
                                    trialRange=trialRange,...
                                    signalRange=signalRange{s},...
                                    trialConditions=trialConditions);
        plotTraces(combined.data{1},combined.timestamp,color=colorList{i});
        xlabel('Time (s)'); ylabel([signalRange{s},' z-score']); 
        ylim([-0.8,2.1]);
        plotEvent(eventRange{i},eventDuration(i),color=colorList{i});
        legend({[eventRange{i},' (n=',num2str(size(combined.data{1},1)),')']},...
                'Location','northeast');
    end
    saveFigures(gcf,'Summary_random_dLight',...
            strcat(resultspath),...
            saveFIG=true,savePDF=true);
end

%% (Baseline) Plot bar plot of DA final values


trialConditions = 'trials.performing';

randomStats = getGroupedTrialStats(combined_random,'stageMax',...
                            eventRange='water',...
                            animalRange='All',...
                            taskRange='Random',...
                            signalRange='dLight', ...
                            trialConditions=trialConditions);
randomResults = plotGroupedTrialStats(randomStats,'',groupSize=1,plot=false);
waterData = randomResults.traces{1}{1};

randomStats = getGroupedTrialStats(combined_random,'stageAvg',...
                            eventRange='stim',...
                            animalRange='All',...
                            taskRange='Random',...
                            signalRange='dLight', ...
                            trialConditions=trialConditions);
randomResults = plotGroupedTrialStats(randomStats,'',groupSize=1,plot=false);
stimData = randomResults.traces{1}{1};

randomStats = getGroupedTrialStats(combined_random,'stageMin',...
                            eventRange='water',...
                            animalRange='All',...
                            taskRange='Random',...
                            signalRange='dLight', ...
                            trialConditions=trialConditions);
randomResults = plotGroupedTrialStats(randomStats,'',groupSize=1,plot=false);
airpuffData = randomResults.traces{1}{1};

initializeFig(0.3,0.5); tiledlayout('flow');

nexttile;
waterData_animal = mean(waterData,2,'omitnan');
stimData_animal = mean(stimData,2,'omitnan');
airpuffData_animal = mean(airpuffData,2,'omitnan');
plotScatterBar(1,waterData_animal,LineWidth=5,dotSize=200,color=bluePurpleRed(1,:),style='bar');
plotScatterBar(2,stimData_animal,LineWidth=5,dotSize=200,color=bluePurpleRed(500,:),style='bar');
plotScatterBar(3,airpuffData_animal,LineWidth=5,dotSize=200,color=[.2 .2 .2],style='bar');
% plotScatterBar([1 2],[waterData_animal,stimData_animal],connectPairs=true,...
%                 LineWidth=5,dotSize=400,color=colorList,style='bar');
% plotScatterBar([2 3],[stimData_animal,airpuffData_animal],connectPairs=true,...
%                 LineWidth=5,dotSize=400,color=colorList,style='bar');
plotStats(waterData_animal,stimData_animal,[1,2],testType='kstest');
plotStats(stimData_animal,airpuffData_animal,[2,3],testType='kstest');
xticks(1:3); 
xticklabels({['Water (n=',num2str(length(waterData_animal)),')'],...
    ['Stim (n=',num2str(length(stimData_animal)),')'],...
    ['Airpuff (',num2str(length(airpuffData_animal)),')']});
ylabel('Amp DA response during cue');

% nexttile;
% waterData_vec = rmmissing(reshape(waterData,[],1));
% stimData_vec = rmmissing(reshape(stimData,[],1));
% airpuffData_vec = rmmissing(reshape(airpuffData,[],1));
% plotScatterBar(1,waterData_vec,LineWidth=5,dotSize=200,color=bluePurpleRed(1,:),style='bar');
% plotScatterBar(2,stimData_vec,LineWidth=5,dotSize=200,color=bluePurpleRed(500,:),style='bar');
% plotScatterBar(3,airpuffData_vec,LineWidth=5,dotSize=200,color=[.2 .2 .2],style='bar');
% plotStats(waterData_vec,stimData_vec,[1,2],testType='kstest');
% plotStats(stimData_vec,airpuffData_vec,[2,3],testType='kstest');
% xticks(1:3); 
% xticklabels({['Water (n=',num2str(length(waterData_vec)),')'],...
%     ['Stim (n=',num2str(length(stimData_vec)),')'],...
%     ['Airpuff (',num2str(length(airpuffData_vec)),')']});
% ylabel('Amp DA response during cue');

saveFigures(gcf,'Summary_DA_randomVal',...
        strcat(resultspath),...
        saveFIG=true,savePDF=true);


%% (Baseline) Last session trace for each animal

timeRange = [-0.5,3];
eventRange = 'Stim';
animalRange = unique({combined_random.animal});

taskRange = 'Random';
totalTrialRange = 'All';
signalRange = 'dLight';

legendList = cell(length(animalRange),1);
colormap = bluePurpleRed;
colorList = colormap(round(linspace(1,length(colormap),length(animalRange))),:);

close all; initializeFig(.5,1);
for a = 1:length(animalRange)
    row = combineTraces(combined_random,searchRows=true,...
                                eventRange=eventRange,...
                                animalRange=animalRange{a},...
                                taskRange=taskRange,...
                                signalRange=signalRange);
    lastSessionStart = combined_random(row).options.startIdx.session{1}(end);
    lastSessionMask = zeros(size(combined_random(row).data,1),1);
    lastSessionMask(lastSessionStart:end) = 1;

    combined = combineTraces(combined_random,timeRange=timeRange,...
                                eventRange=eventRange,...
                                animalRange=animalRange{a},...
                                taskRange=taskRange,...
                                signalRange=signalRange,...
                                trialConditions='trials.performing');

    plotTraces(combined.data{1},combined.timestamp,color=colorList(a,:));
    legendList{a} = [animalRange{a},' (n=',num2str(size(combined.data{1},1)),')'];
end
xlabel('Time (s)'); ylabel([signalRange,' z-score']); ylim([-1,3]);
plotEvent(eventRange,0.5,color=bluePurpleRed(500,:));
legend(legendList,'Location','northeast');
% saveFigures(gcf,'Summary_random_dLight',...
%         strcat(resultspath),...
%         saveFIG=true,savePDF=true);