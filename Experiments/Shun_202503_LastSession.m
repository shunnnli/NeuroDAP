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
    sessionList = uipickfiles('FilterSpec',osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Project valence/Recordings'))';
    groupSessions = true;
    % Update resultspath
    projectName = 'Summary-LastSessions'; 
    resultspath = strcat(resultspath,filesep,projectName);
    % Create resultspath if necessary
    if isempty(dir(resultspath)); mkdir(resultspath); end

elseif strcmpi(answer,'Load combined data')
    fileList = uipickfiles('FilterSpec',osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Project valence/Results'))';
    groupSessions = false;
    % Update resultspath
    projectName = 'Summary-LastSessions'; 
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
    cur_name = summary(i).name;
    cur_animal = summary(i).animal;

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
        summary(i).event = 'Pair';
    end

    if strcmpi('NAc-green',cur_name)
        summary(i).name = 'dLight';
    % elseif strcmpi('PMT',cur_name)
    %     disp(i);
    %     summary(i).name = 'GCaMP8m';
    end
end

% Remove some rows if needed
% rowIdx = cellfun(@(x) contains(x,'Blue stim',IgnoreCase=true), {summary.event});
% summary(rowIdx) = [];

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


%% Test: Plot traces from summary/animals struct

event = 'Stim';
animal = 'All';
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

% plotTraces(combined.data{1},combined.timestamp,color=bluePurpleRed(500,:));
legendList = plotGroupTraces(combined.data{1},combined.timestamp,bluePurpleRed,...
                        groupby='sessions',startIdx=combined.options.startIdx,remaining='include');

plotEvent('Stim',0.5,color=bluePurpleRed(500,:))
xlabel('Time (s)'); ylabel('z-score');


%% Calculate DA trend for each animal

eventRange = 'Stim';
animalRange = unique({animals.animal});
DAtrend = struct([]);

for a = 1:length(animalRange)
    cur_animal = animalRange{a};
    cur_task = unique({animals(strcmp({animals.animal}, cur_animal)).task});
    DA_Max = getGroupedTrialStats(animals,'stageMax',...
                                eventRange=eventRange,...
                                animalRange=cur_animal,...
                                taskRange=cur_task,...
                                signalRange='dLight');
    DA_Min = getGroupedTrialStats(animals,'stageMin',...
                                eventRange=eventRange,...
                                animalRange=cur_animal,...
                                taskRange=cur_task,...
                                signalRange='dLight');
    DA_Avg = getGroupedTrialStats(animals,'stageAvg',...
                                eventRange=eventRange,...
                                animalRange=cur_animal,...
                                taskRange=cur_task,...
                                signalRange='dLight');
    % initializeFig(.7,.7); tiledlayout('flow');
    % DA_stageMax_result = plotGroupedTrialStats(DA_stageMax,'DA Max',groupSize=1,color=bluePurpleRed(500,:));

    % Save to DAtrend
    DAtrend(a).animal = cur_animal;
    DAtrend(a).task   = cur_task;
    DAtrend(a).CueMax = DA_Max.stats{1}{1}(:,2);
    DAtrend(a).CueMin = DA_Min.stats{1}{1}(:,2);
    DAtrend(a).CueAvg = DA_Avg.stats{1}{1}(:,2);
    DAtrend(a).CueAmp = getAmplitude(DA_Max.stats{1}{1}(:,2),DA_Min.stats{1}{1}(:,2));

    % Calculate DA trending stats
    fields = {'CueMax', 'CueMin', 'CueAvg', 'CueAmp'};
    nTrials = length(DAtrend(a).CueMax); % assuming all have the same length
    
    for f = 1:numel(fields)
        data = DAtrend(a).(fields{f});
        slopes = nan(nTrials, 1);
        pvals  = nan(nTrials, 1);
        
        % Loop for trending statistics: last n trials (n = 1:nTrials)
        for n = 2:nTrials
            Y = data(end-n+1:end);  % last n trials
            X = (1:n)';             % corresponding trial indices for regression
            stats = regstats(Y, X, 'linear', {'tstat'});
            slopes(n) = stats.tstat.beta(2);  % slope coefficient
            pvals(n)  = stats.tstat.pval(2);   % p-value for the slope
        end
        
        % Store the trending statistics into DAtrend
        DAtrend(a).([fields{f} '_slope']) = slopes;
        DAtrend(a).([fields{f} '_pval'])  = pvals;
        
        % Compute difference between the last 10 trials and the first 10 trials
        if nTrials >= 10
            DAtrend(a).([fields{f} '_diff']) = mean(data(end-9:end)) - mean(data(1:10));
        else
            DAtrend(a).([fields{f} '_diff']) = NaN;  % Not enough trials
        end
    end
end

%% Save DAtrend struct
prompt = 'Enter database notes (DAtrend_20230326_notes.mat):';
dlgtitle = 'Save DAtrend struct'; fieldsize = [1 45]; definput = {''};
answer = inputdlg(prompt,dlgtitle,fieldsize,definput);
today = char(datetime('today','Format','yyyyMMdd'));
filename = strcat('DAtrend_',today,'_',answer{1});

disp(['Ongoing: saving DAtrend.mat (',char(datetime('now','Format','HH:mm:ss')),')']);
save(strcat(resultspath,filesep,filename),'DAtrend','sessionList','-v7.3');
disp(['Finished: saved DAtrend.mat (',char(datetime('now','Format','HH:mm:ss')),')']);