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

    if contains('random',cur_task,IgnoreCase=true)
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

%% Change signal names based on clamp side

leftClampAnimals = {'SL431', 'SL432', 'SL433', 'BiPOLES2', 'M431'};
rightClampAnimals = {'M445', 'M446'};
NAcLSClampAnimals = {'SL478','SL479','SL480','SL481'};

for i = 1:length(summary)
    cur_animal = summary(i).animal;
    cur_name = summary(i).name;

    if any(strcmpi(cur_animal, leftClampAnimals))
        if strcmpi(cur_name, 'NAc-left')
            summary(i).name = 'NAc-clamp';
        elseif strcmpi(cur_name, 'NAc-right')
            summary(i).name = 'NAc-unclamp';
        end
    elseif any(strcmpi(cur_animal, rightClampAnimals))
        if strcmpi(cur_name, 'NAc-right')
            summary(i).name = 'NAc-clamp';
        elseif strcmpi(cur_name, 'NAc-left')
            summary(i).name = 'NAc-unclamp';
        end
    elseif any(strcmpi(cur_animal, NAcLSClampAnimals))
        if strcmpi(cur_name, 'NAc-LS')
            summary(i).name = 'NAc-clamp';
        end
    end
end


%% Change event name

% If session date is before 20260716, change event name as follows:
% If task is Random-ctrl, add "(unclamp)" to all event name except baseline
% ie "Tone" becomes "Tone (unclamp)"

% If task is Random-clamp, add "(clamp)" to all event name except baseline
% ie "Tone" becomes "Tone (clamp)"

% RandomClampMix sessions interleave clamp and unclamp trials, so their
% event names already carry the correct label and are left untouched.
mixSessionPattern = 'RandomClampMix';

for i = 1:length(summary)
    cur_task = summary(i).task;
    cur_event = summary(i).event;
    cur_date = str2double(summary(i).date);
    cur_session = summary(i).session;

    if ~(cur_date < 20260716) || strcmpi(cur_event, 'Baseline') || ...
            contains(cur_session, mixSessionPattern, IgnoreCase=true)
        continue;
    end

    if strcmpi(cur_task, 'Random-ctrl')
        eventSuffix = ' (unclamp)';
    elseif strcmpi(cur_task, 'Random-clamp')
        eventSuffix = ' (clamp)';
    else
        continue;
    end

    % Avoid adding the suffix again if this section is rerun interactively.
    if ~endsWith(cur_event, eventSuffix, IgnoreCase=true)
        summary(i).event = [cur_event, eventSuffix];
    end
end

% If animal is SL478-481, if Random-ctrl, make sure is unclamp and if task
% is Random-clamp, make sure suffix is (clamp)

% These animals are enforced regardless of session date: any existing
% (clamp)/(unclamp) suffix is stripped and replaced by the one matching the
% task, so wrongly labeled events get corrected too.
suffixAnimals = {'SL478','SL479','SL480','SL481'};

for i = 1:length(summary)
    cur_animal = summary(i).animal;
    cur_task = summary(i).task;
    cur_event = char(summary(i).event);
    cur_session = summary(i).session;

    if ~any(strcmpi(cur_animal, suffixAnimals)) || strcmpi(cur_event, 'Baseline') || ...
            contains(cur_session, mixSessionPattern, IgnoreCase=true)
        continue;
    end

    if strcmpi(cur_task, 'Random-ctrl')
        eventSuffix = ' (unclamp)';
    elseif strcmpi(cur_task, 'Random-clamp')
        eventSuffix = ' (clamp)';
    else
        continue;
    end

    % Remove any existing (clamp)/(unclamp) suffix before adding the right one
    cur_event = char(regexprep(cur_event,'\s*\((un)?clamp\)\s*$','','ignorecase'));
    summary(i).event = [cur_event, eventSuffix];
end

%% Rescale redClamp/blueClamp if neccessary
% Make sure the method is compatible with summary or animals struct

% Sessions loaded before loadSessions auto-detected the clamp full scale
% were normalized by a nominal ADC range that did not match the hardware
% (eg red topping out at 82.6% instead of 100%). rescaleClampTraces finds
% the command level each session actually reaches and maps it onto 100%.

sessionsToScale = 20260927; % scale sessions before this date

[summary,rescaleInfo] = rescaleClampTraces(summary,dateBefore=sessionsToScale);
disp(rescaleInfo);

% animals struct has no date field, so all clamp rows are rescaled:
% animals = rescaleClampTraces(animals);

%% Remove trials where redClamp or blueClamp is max for the whole trial

% In some trials the clamp command is stuck at its maximum for essentially
% the entire trial (>90% of the trace). Instead of deleting these trials,
% set trials.performing = 0 in the trialTable of every signal recorded
% during the same animal/date/session/task/event, so they are excluded by
% any analysis using trialConditions = 'trials.performing'.

clampSignals = {'redClamp','blueClamp'};
clampTasks = {'Random-clamp'};  % only screen clamp sessions, not Random-ctrl
maxFraction = 0.9;      % flag trial if clamp is at max for > this fraction of the trial
saturationPct = 100;    % clamp commands are stored as 0-100% of the clamp range
maxTolerance = 1;       % (%) how close to saturationPct still counts as "at max"

% Group rows that share the same set of trials
groupLabels = strings(length(summary),1);
for i = 1:length(summary)
    groupLabels(i) = strjoin(string({summary(i).animal, summary(i).date,...
                                     summary(i).session, summary(i).task,...
                                     summary(i).event}),'_');
end
[uniqueGroups,~,groupIdx] = unique(groupLabels);

nRemovedTotal = 0; nTrialsTotal = 0;
for g = 1:length(uniqueGroups)
    groupRows = reshape(find(groupIdx == g),1,[]);
    if ~any(strcmpi(summary(groupRows(1)).task, clampTasks)); continue; end
    signalNames = string({summary(groupRows).name});
    clampRows = groupRows(ismember(lower(signalNames),lower(clampSignals)));
    if isempty(clampRows); continue; end

    % Flag saturated trials (union across redClamp & blueClamp)
    nTrials = size(summary(clampRows(1)).data,1);
    badTrials = false(nTrials,1);
    for r = clampRows
        clampData = summary(r).data;
        if size(clampData,1) ~= nTrials
            warning(['Skipped ',char(uniqueGroups(g)),' -> ',summary(r).name,...
                     ': trial number mismatch within group']);
            continue
        end
        atMax = clampData >= (saturationPct - maxTolerance);
        fracAtMax = sum(atMax,2) ./ sum(~isnan(clampData),2);
        badTrials = badTrials | (fracAtMax > maxFraction);
    end
    nTrialsTotal = nTrialsTotal + nTrials;
    if ~any(badTrials); continue; end
    nRemovedTotal = nRemovedTotal + sum(badTrials);

    % Set performing = 0 for flagged trials in every signal of this group
    for r = groupRows
        cur_table = summary(r).trialInfo.trialTable;
        if height(cur_table) ~= nTrials
            warning(['Skipped ',char(uniqueGroups(g)),' -> ',summary(r).name,...
                     ': trialTable height does not match trial number']);
            continue
        end
        if ~ismember('performing',cur_table.Properties.VariableNames)
            cur_table.performing = ones(nTrials,1);
        end
        cur_table.performing(badTrials) = 0;
        summary(r).trialInfo.trialTable = cur_table;
    end

    disp(['Finished: flagged ',num2str(sum(badTrials)),'/',num2str(nTrials),...
          ' saturated clamp trials in ',char(uniqueGroups(g))]);
end
disp(['Finished: flagged ',num2str(nRemovedTotal),'/',num2str(nTrialsTotal),...
      ' trials (performing = 0) where redClamp or blueClamp was at max for >',...
      num2str(maxFraction*100),'% of the trial']);


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

%% Animal groups

DLS1 = {'SL479','SL480','SL481'};


%% Random: clamp vs unclamp DLS 

close all;
timeRange = [-0.5,3];
eventRange = {'Rewarded licks','Airpuff','Tone'};
animalRange = DLS1;%{'SL431','SL432','SL433','BiPOLES2'};
signalRange = {'NAc-clamp','DLS'};
trialConditions = 'trials.performing';
sessionRange = 'RandomClampMix';    % only the interleaved clamp/unclamp sessions

colorList = {bluePurpleRed(1,:),[.2,.2,.2],bluePurpleRed(100,:)};
eventDuration = [0,.2,.5];

close all; 
for i = 1:length(eventRange)
    initializeFig(.5,.5); tiledlayout(1,length(signalRange));
    for s = 1:length(signalRange)
        nexttile;
        combined = combineTraces(animals,timeRange=timeRange,...
                                    eventRange=[eventRange{i},' (unclamp)'],...
                                    animalRange=animalRange,...
                                    taskRange='Random',...
                                    sessionRange=sessionRange,...
                                    signalRange=signalRange{s},...
                                    trialConditions=trialConditions);
        plotTraces(combined.data{1},combined.timestamp,color=unclampColor);

        combined = combineTraces(animals,timeRange=timeRange,...
                                    eventRange=[eventRange{i},' (clamp)'],...
                                    animalRange=animalRange,...
                                    taskRange='Random',...
                                    sessionRange=sessionRange,...
                                    signalRange=signalRange{s},...
                                    trialConditions=trialConditions);
        plotTraces(combined.data{1},combined.timestamp,color=clampColor);
        xlabel('Time (s)'); ylabel([signalRange{s},' (\DeltaF/F)']); 
        ylim([-0.02,0.15]);
        plotEvent(eventRange{i},eventDuration(i),color=colorList{i});
        legend({[eventRange{i},' (n=',num2str(size(combined.data{1},1)),')']},...
                'Location','northeast');
    end
    % saveFigures(gcf,strcat('Summary_random_',eventRange{i}),...
    %         strcat(resultspath),...
    %         saveFIG=false,savePDF=true);
end

%% ScatterBar plot showing the stageAmp for clamp vs unclamp for each event, one panel for each region

close all;
timeRange = [-0.5,3];
eventRange = {'Water','Airpuff','Tone'};
animalRange = {'SL478','SL479','SL480','SL481'};
signalRange = {'NAc-clamp','DLS'};
trialConditions = 'trials.performing';
conditionSuffix = {' (clamp)',' (unclamp)'};
conditionColor = [clampColor; unclampColor];
stage = 2;      % stage window used for stageAmp (same default as plotGroupedTrialStats)
sessionRange = 'RandomClampMix';    % only the interleaved clamp/unclamp sessions

initializeFig(.5,.5); tiledlayout(1,length(signalRange));
for s = 1:length(signalRange)
    nexttile;
    for i = 1:length(eventRange)
        % Each dot is one trial, pooled across animals
        ampData = cell(1,length(conditionSuffix));
        for c = 1:length(conditionSuffix)
            combined = combineTraces(animals,timeRange=timeRange,...
                                        eventRange=[eventRange{i},conditionSuffix{c}],...
                                        animalRange=animalRange,...
                                        taskRange='Random',...
                                        sessionRange=sessionRange,...
                                        signalRange=signalRange{s},...
                                        statsType='stageAmp',...
                                        trialConditions=trialConditions);
            if combined.options.empty || isempty(combined.stats.stageAmp{1}); continue; end
            trialAmp = combined.stats.stageAmp{1}(:,stage);     % one value per trial
            ampData{c} = trialAmp(~isnan(trialAmp));
        end
        if any(cellfun(@isempty,ampData)); continue; end

        x = [2*i-1, 2*i];
        for c = 1:length(conditionSuffix)
            plotScatterBar(x(c),ampData{c},style='bar',color=conditionColor(c,:),...
                           dotSize=20,MarkerFaceAlpha=0.5,LineWidth=2);
        end
        plotStats(ampData{1},ampData{2},x,testType='ranksum');  % unpaired across trials
    end
    yline(0,'--',Color=[.7 .7 .7],HandleVisibility='off');
    xlim([0.5,2*length(eventRange)+0.5]); xticks(1:2*length(eventRange));
    xticklabels(reshape([strcat(eventRange,' (clamp)');strcat(eventRange,' (unclamp)')],1,[]));
    ylabel([signalRange{s},' stageAmp (\DeltaF/F)']); title(signalRange{s});
end
% saveFigures(gcf,'Summary_random_stageAmp_clampVsUnclamp',...
%         strcat(resultspath),...
%         saveFIG=false,savePDF=true);

%% Same figure but clamp trials only from RandomClampMix sessions

% combineTraces can select sessions in the animals struct as well: it matches
% options.sessionRange against the sessionList recorded in each row and keeps
% the corresponding trials. Clamp AND unclamp trials both come from the
% RandomClampMix sessions here, where the two conditions are interleaved.
close all;
timeRange = [-0.5,3];
eventRange = {'Water','Airpuff','Tone'};
animalRange = {'SL478','SL479','SL480','SL481'};
signalRange = {'NAc-clamp','DLS'};
trialConditions = 'trials.performing';
conditionSuffix = {' (clamp)',' (unclamp)'};
conditionColor = [clampColor; unclampColor];
stage = 2;
mixSessionPattern = 'RandomClampMix';
sessionRange = mixSessionPattern;

initializeFig(.5,.5); tiledlayout(1,length(signalRange));
for s = 1:length(signalRange)
    nexttile;
    for i = 1:length(eventRange)
        % Each dot is one trial, pooled across animals
        ampData = cell(1,length(conditionSuffix));
        for c = 1:length(conditionSuffix)
            combined = combineTraces(animals,timeRange=timeRange,...
                                        eventRange=[eventRange{i},conditionSuffix{c}],...
                                        animalRange=animalRange,...
                                        taskRange='Random',...
                                        sessionRange=sessionRange,...
                                        signalRange=signalRange{s},...
                                        statsType='stageAmp',...
                                        trialConditions=trialConditions);
            if combined.options.empty || isempty(combined.stats.stageAmp{1}); continue; end
            trialAmp = combined.stats.stageAmp{1}(:,stage);     % one value per trial
            ampData{c} = trialAmp(~isnan(trialAmp));
        end
        if any(cellfun(@isempty,ampData)); continue; end

        x = [2*i-1, 2*i];
        for c = 1:length(conditionSuffix)
            plotScatterBar(x(c),ampData{c},style='bar',color=conditionColor(c,:),...
                           dotSize=20,MarkerFaceAlpha=0.5,LineWidth=2);
        end
        plotStats(ampData{1},ampData{2},x,testType='ranksum');  % unpaired across trials
    end
    yline(0,'--',Color=[.7 .7 .7],HandleVisibility='off');
    xlim([0.5,2*length(eventRange)+0.5]); xticks(1:2*length(eventRange));
    xticklabels(reshape([strcat(eventRange,' (clamp)');strcat(eventRange,' (unclamp)')],1,[]));
    ylabel([signalRange{s},' stageAmp (\DeltaF/F)']);
    title([signalRange{s},' (',mixSessionPattern,' sessions)']);
end
% saveFigures(gcf,'Summary_random_stageAmp_mixSessions',...
%         strcat(resultspath),...
%         saveFIG=false,savePDF=true);

%% Same figure but clamp trials not from RandomClampMix sessions

% Same as above with the session filter inverted (ie block-design sessions
% where a whole session is either clamp or unclamp). sessionRange only
% matches, never excludes, so list every non-mix session name recorded in
% animals(i).options.sessionList.
close all;
timeRange = [-0.5,3];
eventRange = {'Water','Airpuff','Tone'};
animalRange = {'SL478','SL479','SL480','SL481'};
signalRange = {'NAc-clamp','DLS'};
trialConditions = 'trials.performing';
conditionSuffix = {' (clamp)',' (unclamp)'};
conditionColor = [clampColor; unclampColor];
stage = 2;
mixSessionPattern = 'RandomClampMix';

sessionLists = arrayfun(@(x) string(x.options.sessionList),animals,UniformOutput=false);
allSessions = unique([sessionLists{:}]);
sessionRange = allSessions(~contains(allSessions,mixSessionPattern,IgnoreCase=true));

initializeFig(.5,.5); tiledlayout(1,length(signalRange));
for s = 1:length(signalRange)
    nexttile;
    for i = 1:length(eventRange)
        % Each dot is one trial, pooled across animals
        ampData = cell(1,length(conditionSuffix));
        for c = 1:length(conditionSuffix)
            combined = combineTraces(animals,timeRange=timeRange,...
                                        eventRange=[eventRange{i},conditionSuffix{c}],...
                                        animalRange=animalRange,...
                                        taskRange='Random',...
                                        sessionRange=sessionRange,...
                                        signalRange=signalRange{s},...
                                        statsType='stageAmp',...
                                        trialConditions=trialConditions);
            if combined.options.empty || isempty(combined.stats.stageAmp{1}); continue; end
            trialAmp = combined.stats.stageAmp{1}(:,stage);     % one value per trial
            ampData{c} = trialAmp(~isnan(trialAmp));
        end
        if any(cellfun(@isempty,ampData)); continue; end

        x = [2*i-1, 2*i];
        for c = 1:length(conditionSuffix)
            plotScatterBar(x(c),ampData{c},style='bar',color=conditionColor(c,:),...
                           dotSize=20,MarkerFaceAlpha=0.5,LineWidth=2);
        end
        plotStats(ampData{1},ampData{2},x,testType='ranksum');  % unpaired across trials
    end
    yline(0,'--',Color=[.7 .7 .7],HandleVisibility='off');
    xlim([0.5,2*length(eventRange)+0.5]); xticks(1:2*length(eventRange));
    xticklabels(reshape([strcat(eventRange,' (clamp)');strcat(eventRange,' (unclamp)')],1,[]));
    ylabel([signalRange{s},' stageAmp (\DeltaF/F)']);
    title([signalRange{s},' (block sessions)']);
end
% saveFigures(gcf,'Summary_random_stageAmp_blockSessions',...
%         strcat(resultspath),...
%         saveFIG=false,savePDF=true);
