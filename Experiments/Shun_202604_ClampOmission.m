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
clampColor = [.232 .76 .58]; % #3bc294
unclampColor = [165, 209, 178]./255;
toneColor = bluePurpleRed(100,:);

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

%% Optional: for learning sessions only

% Change task to reward
for i = 1:length(summary)
    cur_session = summary(i).session;
    cur_animal  = summary(i).animal;
    cur_date    = str2double(summary(i).date);

    if contains(cur_session,"omission",IgnoreCase=true)
        if contains(cur_session,"unclamp",IgnoreCase=true)
            summary(i).task = 'Reward-Unclamp-Omission';
        elseif contains(cur_session,"clamp",IgnoreCase=true)
            summary(i).task = 'Reward-Clamp-Omission';
        else
            summary(i).task = 'Reward-Unclamp-Omission';
        end
    elseif contains(cur_session,"bringbackreward",IgnoreCase=true)
        summary(i).task = 'Reward-Unclamp-BringBackReward';
    elseif contains(cur_session,["unclamp","ctrl"],IgnoreCase=true)
        summary(i).task = 'Reward-Unclamp';
    elseif contains(cur_session,"wholeTrial",IgnoreCase=true)
        summary(i).task = 'Reward-Clamp-wholeTrial';
    elseif contains(cur_session,"delayReward",IgnoreCase=true)
        summary(i).task = 'Reward-Clamp-delayReward';
    elseif contains(cur_session,"RPE",IgnoreCase=true)
        summary(i).task = 'Reward-Clamp-withRPE';
    else
        summary(i).task = 'Reward-Unclamp';
    end

    if cur_date == 20260521 && any(strcmpi(cur_animal, ["SL438", "SL439"]))
        summary(i).task = 'Reward-Unclamp-BringBackReward';
    end
end

%% Change / add more details to task
for i = 1:length(summary)
    cur_task    = summary(i).task;
    cur_animal  = summary(i).animal;
    cur_date    = str2double(summary(i).date);
    cur_session = summary(i).session;

    % Default for all non-clamp/control animals
    if ~any(strcmpi(cur_animal, {'SL433','SL431','SL432','M431'}))
        summary(i).task = 'Reward-Ctrl';
        continue
    end

    % SL433, SL431, SL432
    if cur_date < 20260423
        summary(i).task = 'Reward-Clamp-wholeTrial';
    elseif cur_date == 20260423
        summary(i).task = 'Reward-Clamp-delayReward';
    % SL433-specific rules
    elseif strcmpi(cur_animal, 'SL433') && cur_date >= 20260424 && cur_date <= 20260426
        summary(i).task = 'Reward-Clamp-withRPE';
    elseif strcmpi(cur_animal, 'SL433') && cur_date >= 20260427 && cur_date <= 20260428
        summary(i).task = 'Reward-Unclamp';
    % SL431-specific rules
    elseif strcmpi(cur_animal, 'SL431') && cur_date >= 20260513 && cur_date <= 20260517
        summary(i).task = 'Reward-Clamp-withRPE';
    elseif strcmpi(cur_animal, 'SL431') && cur_date >= 20260518 && cur_date <= 20260520
        summary(i).task = 'Reward-Unclamp';
    else
        % Clamp animal, but date does not match any rule.
        % Keep the existing task unchanged.
        summary(i).task = cur_task;
    end
end

keepRows = true(1, length(summary));
for i = 1:length(summary)

    cur_animal = string(summary(i).animal);
    cur_name   = string(summary(i).name);
    cur_task   = string(summary(i).task);
    cur_event  = string(summary(i).event);
    cur_date   = string(summary(i).date);

    skipGroup1 = any(strcmpi(cur_animal, ["SL438", "SL439"])) && ...
                 any(strcmpi(cur_name, ["NAc-left", "NAc-rightLS"]));

    skipGroup2 = any(strcmpi(cur_animal, ["SL446","SL447","SL443","SL444","SL445","SL438","SL439","M430"])) && ...
                 any(strcmpi(cur_name, ["blueClamp", "redClamp"]));

    skipGroup3 = any(strcmpi(cur_animal, ["SL431", "SL432", "SL433"])) && ...
                 any(strcmpi(cur_task, ["Reward-Clamp-wholeTrial", ...
                                         "Reward-Clamp-delayReward", ...
                                         "Reward-Clamp-withRPE"])) && ...
                 strcmpi(cur_event, "Tone (unclamp)");

    skipGroup4 = any(strcmpi(cur_animal, ["SL443", "SL444"])) && ...
                 any(strcmpi(cur_name, "NAc-right")) && ...
                 any(strcmpi(cur_date, "20260613"));

    skipGroup5 = any(strcmpi(cur_animal, "M431")) && ...
                 any(strcmpi(cur_name, "NAc-right"));

    skipGroup6 = any(strcmpi(cur_animal, "M430")) && ...
                 any(strcmpi(cur_name, "NAc-left")) && ...
                 any(strcmpi(cur_date, "20260714"));

    if skipGroup1 || skipGroup2 || skipGroup3 || skipGroup4 || skipGroup5
        keepRows(i) = false;
    end
end
summary = summary(keepRows);

% Make sure SL431 rightLS should be NAc-left

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

%% *********************** Learning stage code ***********************
clampAnimals = {'SL431','SL432','SL433','M431'}; % 'SL431','SL432','SL433',
ctrlAnimals = {};

if isempty(ctrlAnimals)
    allAnimals = unique({animals.animal}, 'stable');
    ctrlAnimals = setdiff(allAnimals, clampAnimals, 'stable'); 
end

animalRange = {clampAnimals; ctrlAnimals};
animalTypes = {'Clamped','Ctrl'};

% Session color
clampColor_wholeTrial = clampColor;
clampColor_delayReward = [169, 178, 82]./255; % a9b252
clampColor_withRPE = [0, 118, 77]./255;

% TODO: make nSessions not hard coded
nSessions = 12; 
% TODO: make colormap change opacity within the same type of sessions
sessionColormap_clamped = [
    repmat(clampColor_wholeTrial, 3, 1)
    clampColor_delayReward
    repmat(clampColor_withRPE, 3, 1)
];
sessionColormap_unclamped = getColormap(unclampColor*255,clampColor_withRPE*255,nSessions,'midCol',clampColor*255);

%% Plot DA trace aligned to Tone (unfinished)

timeRange = [-0.5,3];
eventRange = {'Tone (clamp)','Tone (unclamp)'};
signalRange = {'blueClamp','NAc'};
taskRange = {'Reward-Clamp-wholeTrial','Reward-Ctrl'};

totalTrialRange = 'All';
trialRange = 'All';
groupSizeList = 30;
nGroupsList = 5;

initializeFig(.3,.5); tiledlayout('flow');
for type = 1:length(animalRange)
    if type == 1; sessionColormap = sessionColormap_clamped;
    else; sessionColormap = sessionColormap_unclamped; end

    nexttile;
    combined = combineTraces(animals,timeRange=timeRange,...
                                eventRange=eventRange{type},...
                                taskRange=taskRange{type},...
                                animalRange=animalRange{type},...
                                totalTrialRange=totalTrialRange,...
                                trialRange=trialRange,...
                                signalRange=signalRange{type});
    legendList = plotGroupTraces(combined.data{1},combined.timestamp,sessionColormap,...
                    groupSize=groupSizeList,nGroups=nGroupsList,...
                    groupby='session',startIdx=combined.options.startIdx);

    title(animalTypes{type});
    xlabel('Time (s)'); ylabel('dff'); %ylim([-0.1,0.5]);
    plotEvent('Tone',.5,color=toneColor);
    legend(legendList,'Location','northeast');
end
% saveFigures(gcf,['Summary_licking_',taskRange{task}],...
%     strcat(resultspath),...
%     saveFIG=true,savePDF=true);

%% Get first 10 trials of DA response for Reward-Clamp-delayReward

timeRange = [-0.5,3];
first10Trials = 1:10;
totalTrialRange = 'All';
trialRange = 'All';
ctrlSessionNum = 4;
stackSpacingScale = 0.5; % smaller = tighter stacked trial rows

initializeFig(.4,.6);
tiledlayout(1,2,Padding='compact');

combined_clamp = combineTraces(animals,timeRange=timeRange,...
                                eventRange='Tone (clamp)',...
                                taskRange='Reward-Clamp-delayReward',...
                                animalRange=animalRange{1},...
                                totalTrialRange=totalTrialRange,...
                                trialRange=trialRange,...
                                signalRange='NAc-left');

combined_ctrl = combineTraces(animals,timeRange=timeRange,...
                                eventRange='Tone (unclamp)',...
                                taskRange='Reward-Ctrl',...
                                animalRange=animalRange{2},...
                                totalTrialRange=totalTrialRange,...
                                trialRange=trialRange,...
                                signalRange='NAc');

finiteData = [combined_clamp.data{1}(isfinite(combined_clamp.data{1})); ...
              combined_ctrl.data{1}(isfinite(combined_ctrl.data{1}))];
stackSpacing = max((max(finiteData)-min(finiteData))*stackSpacingScale,0.1);

% Clamp animals: first 10 trials
nexttile;
sessionStarts = combined_clamp.options.startIdx.session{1};
nextSessionStarts = [sessionStarts(2:end); size(combined_clamp.data{1},1)+1];
for trial = first10Trials
    trialRows = sessionStarts + trial - 1;
    validRows = trialRows < nextSessionStarts & trialRows <= size(combined_clamp.data{1},1);
    plotSEM(combined_clamp.timestamp,combined_clamp.data{1}(trialRows(validRows),:) + (trial-1)*stackSpacing,...
            clampColor_delayReward,plotIndividual=true,plotPatch=false);
end
plotEvent('Tone',.5,color=toneColor);
yticks((first10Trials-1)*stackSpacing); yticklabels(string(first10Trials));
xlabel('Time (s)'); ylabel('Trial'); title('Clamp: delay reward');
xlim(timeRange); ylim([-0.35*stackSpacing,first10Trials(end)*stackSpacing]);

% Ctrl animals: first 10 trials of session 4
nexttile;
sessionStarts = combined_ctrl.options.startIdx.session{1};
animalStarts = combined_ctrl.options.startIdx.animal{1};
ctrlSessionStarts = [];
for animal = 1:numel(animalStarts)
    if animal < numel(animalStarts); animalEnd = animalStarts(animal+1)-1;
    else; animalEnd = size(combined_ctrl.data{1},1); end
    curSessionStarts = sessionStarts(sessionStarts >= animalStarts(animal) & sessionStarts <= animalEnd);
    if numel(curSessionStarts) >= ctrlSessionNum
        ctrlSessionStarts = [ctrlSessionStarts; curSessionStarts(ctrlSessionNum)];
    end
end
[~,ctrlSessionIdx] = ismember(ctrlSessionStarts,sessionStarts);
nextSessionStarts = [sessionStarts(2:end); size(combined_ctrl.data{1},1)+1];
for trial = first10Trials
    trialRows = ctrlSessionStarts + trial - 1;
    validRows = trialRows < nextSessionStarts(ctrlSessionIdx) & trialRows <= size(combined_ctrl.data{1},1);
    plotSEM(combined_ctrl.timestamp,combined_ctrl.data{1}(trialRows(validRows),:) + (trial-1)*stackSpacing,...
            unclampColor,plotIndividual=true,individualColor=addOpacity(unclampColor,0.35),...
            plotPatch=false);
end
plotEvent('Tone',.5,color=toneColor);
yticks((first10Trials-1)*stackSpacing); yticklabels(string(first10Trials));
xlabel('Time (s)'); ylabel('Trial'); title(['Ctrl: session ',num2str(ctrlSessionNum)]);
xlim(timeRange); ylim([-0.35*stackSpacing,first10Trials(end)*stackSpacing]);

% saveFigures(gcf,'Summary_first10Trials_NAc',...
%         strcat(resultspath),...
%         saveFIG=true,savePDF=true);

%% Plot lick trace

timeRange = [-0.5,3];

totalTrialRange = 'All';
trialRange = 'All';

initializeFig(.6,.5); tiledlayout('flow');

% Plot clamp animals
nexttile;
combined = combineTraces(animals,timeRange=timeRange,...
                            eventRange='Tone (clamp)',...
                            taskRange='Reward-Clamp',...
                            animalRange=animalRange{1},...
                            totalTrialRange=totalTrialRange,...
                            trialRange=trialRange,...
                            signalRange='Lick');
legendList = plotGroupTraces(combined.data{1},combined.timestamp,sessionColormap_clamped,...
                groupby='session',startIdx=combined.options.startIdx);
combined = combineTraces(animals,timeRange=timeRange,...
                            eventRange='Tone (unclamp)',...
                            taskRange='Reward-Unclamp',...
                            animalRange=animalRange{1},...
                            totalTrialRange=totalTrialRange,...
                            trialRange=trialRange,...
                            signalRange='Lick');
plotSEM(combined.timestamp,combined.data{1},unclampColor);
legendList{end+1} = strcat("Unclamp (n=",num2str(size(combined.data{1},1)),")");
plotEvent('Tone',.5,color=toneColor);

title(animalTypes{1});
xlabel('Time (s)'); ylabel('Licks/s'); ylim([0 Inf]);
legend(legendList,'Location','northeast');

% Plot ctrl animals
nexttile;
combined = combineTraces(animals,timeRange=timeRange,...
                            eventRange='Tone',...
                            taskRange='Reward-Ctrl',...
                            animalRange=animalRange{2},...
                            totalTrialRange=totalTrialRange,...
                            trialRange=trialRange,...
                            signalRange='Lick');
legendList = plotGroupTraces(combined.data{1},combined.timestamp,bluePurpleRed,...
                groupby='session',startIdx=combined.options.startIdx);
plotEvent('Tone',.5,color=toneColor);

title(animalTypes{2});
xlabel('Time (s)'); ylabel('Licks/s'); ylim([0 Inf]);
legend(legendList,'Location','northeast');

% saveFigures(gcf,['Summary_licking_',taskRange{task}],...
%     strcat(resultspath),...
%     saveFIG=true,savePDF=true);

%% Plot grouped anticipatory lick changes
% TODO: modify plotGroupedTrialsStats so that it can group based on
% sessions

taskRange = {'Reward'};
statsType = 'nAnticipatoryLicks';

conditionRange = 'All';
signalRange = 'Lick';
trialConditions = 'trials.performing';
eventRange = {'Tone'};
groupSize = 20; % numbers of trials to calculate average

if strcmpi(statsType,'nAnticipatoryLicks')
    ylabelList = 'Anticipatory licks';
else
    ylabelList = statsType;
end

initializeFig(.3,.6); tiledlayout('flow');

% Clamped animals
combinedStats = getGroupedTrialStats(animals,statsType,...
                            eventRange=eventRange,...
                            animalRange=animalRange{1},...
                            taskRange=taskRange,...
                            totalTrialRange=conditionRange,...
                            signalRange=signalRange, ...
                            trialConditions=trialConditions);
animalList_clamped = intersect(unique({animals(contains({animals.task},taskRange{1},IgnoreCase=true)).animal}),animalRange{1});
clampTrialTask = cell(numel(animalList_clamped),1);
for animal = 1:numel(animalList_clamped)
    rowIdx = find(contains({animals.animal},animalList_clamped{animal},IgnoreCase=true) & ...
                  contains({animals.task},taskRange{1},IgnoreCase=true) & ...
                  contains({animals.event},eventRange{1},IgnoreCase=true) & ...
                  contains({animals.name},signalRange,IgnoreCase=true));
    [~,sortIdx] = sort(arrayfun(@(x) min(x.options.signalRows),animals(rowIdx)));
    combinedStats.stats{1}{animal,1} = [];
    for row = rowIdx(sortIdx)
        trialIdx = 1:numel(animals(row).trialInfo.trialNumber);
        if ~(ischar(conditionRange) || isstring(conditionRange))
            trialIdx = find(animals(row).trialInfo.trialNumber >= conditionRange(1) & ...
                            animals(row).trialInfo.trialNumber <= conditionRange(2));
        end
        if ~isempty(trialConditions)
            trials = animals(row).trialInfo.trialTable;
            trialIdx = intersect(trialIdx,find(eval(trialConditions)));
        end
        combinedStats.stats{1}{animal,1} = [combinedStats.stats{1}{animal,1}; animals(row).trialInfo.trialTable.(statsType)(trialIdx)];
        clampTrialTask{animal} = [clampTrialTask{animal}; repmat(string(animals(row).task),numel(trialIdx),1)];
    end
end
results_clamped = plotGroupedTrialStats(combinedStats,ylabelList,groupSize=groupSize,color=clampColor,...
                                plotIndividual=false, plotCommonTrials=false, plot=false);

% Ctrl animals
combinedStats = getGroupedTrialStats(animals,statsType,...
                            eventRange=eventRange,...
                            animalRange=animalRange{2},...
                            taskRange=taskRange,...
                            totalTrialRange=conditionRange,...
                            signalRange=signalRange, ...
                            trialConditions=trialConditions);
results_unclamped = plotGroupedTrialStats(combinedStats,ylabelList,groupSize=groupSize,color=unclampColor,...
                                plotIndividual=false, plotCommonTrials=false, plot=false);

nexttile;
clampTaskRange = {'Reward-Clamp-wholeTrial','Reward-Clamp-delayReward','Reward-Clamp-withRPE','Reward-Unclamp','Reward'};
clampTaskLabel = {'Clamp whole trial','Clamp delay reward','Clamp with RPE','Unclamp','Reward'};
clampTaskColors = [clampColor_wholeTrial; clampColor_delayReward; clampColor_withRPE; unclampColor; clampColor];

clampGrouped = results_clamped.traces{1}{1};
x_clamp = groupSize:groupSize:groupSize*size(clampGrouped,2);
clampGroupTask = cell(numel(clampTrialTask),numel(x_clamp));
for animal = 1:numel(clampTrialTask)
    for group = 1:numel(x_clamp)
        idx = (group-1)*groupSize+1:min(group*groupSize,numel(clampTrialTask{animal}));
        if ~isempty(idx)
            clampGroupTask{animal,group} = char(mode(categorical(clampTrialTask{animal}(idx))));
        end
    end
end

clampGroupTaskMode = strings(1,numel(x_clamp));
for group = 1:numel(x_clamp)
    taskList = clampGroupTask(~cellfun(@isempty,clampGroupTask(:,group)),group);
    if ~isempty(taskList); clampGroupTaskMode(group) = string(mode(categorical(taskList))); end
end

ctrlGrouped = results_unclamped.traces{1}{1};
x_ctrl = groupSize:groupSize:groupSize*size(ctrlGrouped,2);
scatter(x_ctrl,mean(ctrlGrouped,1,'omitnan'),200,unclampColor,'filled',HandleVisibility='off'); hold on
plotSEM(x_ctrl,ctrlGrouped,unclampColor,label='Ctrl');

for task = 1:numel(clampTaskRange)
    groupIdx = find(strcmpi(clampGroupTaskMode,clampTaskRange{task}));
    if isempty(groupIdx); continue; end
    scatter(x_clamp(groupIdx),mean(clampGrouped(:,groupIdx),1,'omitnan'),200,clampTaskColors(task,:),'filled',HandleVisibility='off');
    segmentBreaks = [0,find(diff(groupIdx)>1),numel(groupIdx)];
    legendLabel = clampTaskLabel{task};
    for segment = 1:numel(segmentBreaks)-1
        idx = groupIdx(segmentBreaks(segment)+1:segmentBreaks(segment+1));
        plotSEM(x_clamp(idx),clampGrouped(:,idx),clampTaskColors(task,:),label=legendLabel);
        legendLabel = '';
    end
end
xlabel('Trials'); ylabel(ylabelList); 
xlim([0,800]);%xlim([0,max(x_clamp)]);
legend('Location','northeast');

% saveFigures(gcf,'Summary_CSvsTrialsGrouped_Lick',...
%         strcat(resultspath),...
%         saveFIG=true,savePDF=true);

%% Plot bar plot of anticipatory lick slopes (unmodified)

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


%% Plot grouped CS DA response (check correctness)

taskRange = {'Reward'};
statsTypes = {'stageAmp'}; ylabelList = {'Amp DA response during cue'};
ylimit = [-1,4; -1,5];

conditionRange = 'All';
trialRange = 'All';
trialConditions = 'trials.performing == 1';
eventRange = {'Tone'};
groupSize = 20;

initializeFig(.3,.6); tiledlayout('flow');

% Clamped animals
combinedStats = getGroupedTrialStats(animals,statsType,...
                            eventRange=eventRange,...
                            animalRange=animalRange{1},...
                            taskRange=taskRange,...
                            totalTrialRange=conditionRange,...
                            signalRange='NAc-left', ...
                            trialConditions=trialConditions);
results_clamped = plotGroupedTrialStats(combinedStats,ylabelList,groupSize=groupSize,color=clampColor,...
                                plotIndividual=false, plotCommonTrials=false);

% Unclamp animals
combinedStats = getGroupedTrialStats(animals,statsType,...
                            eventRange=eventRange,...
                            animalRange=animalRange{2},...
                            taskRange=taskRange,...
                            totalTrialRange=conditionRange,...
                            signalRange='NAc-right', ...
                            trialConditions=trialConditions);
results_unclamped = plotGroupedTrialStats(combinedStats,ylabelList,groupSize=groupSize,color=unclampColor,...
                                plotIndividual=false, plotCommonTrials=false, plotNextTile=false);

% saveFigures(gcf,'Summary_CSvsTrialsGrouped_dLight',...
%         strcat(resultspath),...
%         saveFIG=true,savePDF=true);

%% Plot bar plot of DA slopes (unmodified)

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
