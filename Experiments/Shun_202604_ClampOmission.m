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

% Change / add more details to task
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

    if skipGroup1 || skipGroup2 || skipGroup3 || skipGroup4 || skipGroup5 || skipGroup6
        keepRows(i) = false;
    end
end
summary = summary(keepRows);

%% Change NAc signal name
% Change NAc-left to NAc-clamped for M431, SL431, SL433
% Change NAc-right to NAc-ctrl for SL431, SL433
for i = 1:length(summary)
    cur_animal = string(summary(i).animal);
    cur_name   = string(summary(i).name);

    if any(strcmpi(cur_animal, ["M431", "SL431", "SL433"])) && strcmpi(cur_name, "NAc-clamped")
        summary(i).name = 'NAc-clamp';
    elseif any(strcmpi(cur_animal, ["SL431", "SL433"])) && strcmpi(cur_name, "NAc-ctrl")
        summary(i).name = 'NAc-unclamp';
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
clampAnimals = {'SL431','SL433','M431'}; % 'SL431','SL432','SL433',
ctrlAnimals = {};

if isempty(ctrlAnimals)
    allAnimals = unique({animals.animal}, 'stable');
    ctrlAnimals = setdiff(allAnimals, clampAnimals, 'stable'); 
end

animalRange = {clampAnimals; ctrlAnimals};
animalTypes = {'Clamped','Ctrl'};

% Session color
bringBackRewardColor = bluePurpleRed(50,:);
clampColor_wholeTrial = clampColor;
clampColor_delayReward = [169, 178, 82]./255; % a9b252
clampColor_withRPE = [0, 118, 77]./255;
ctrlColor = [154, 176, 166]./255; % #9ab0a6

% TODO: make nSessions not hard coded
nSessions = 12; 
% TODO: make colormap change opacity within the same type of sessions
sessionColormap_clamped = [
    repmat(clampColor_wholeTrial, 3, 1)
    clampColor_delayReward
    repmat(clampColor_withRPE, 3, 1)
];
sessionColormap_unclamped = getColormap(unclampColor*255,clampColor_withRPE*255,nSessions,'midCol',clampColor*255);

%% Plot DA trace aligned to Tone during omission

timeRange = [-0.5,3];
signalRange = {'NAc-clamp','NAc-right'};

eventRange = {'Tone (clamp omission)','Tone (unclamp omission)'};
taskRange = {'Reward-Clamp-Omission','Reward-Unclamp-Omission'};

totalTrialRange = 'All';
trialRange = 'All';
nGroupsList = 4;

initializeFig(.6,.5); tiledlayout('flow');
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
                    nGroups=nGroupsList,...
                    groupby='session',startIdx=combined.options.startIdx);

    title(animalTypes{type});
    xlabel('Time (s)'); ylabel('dff'); ylim([-0.02,0.16]);
    plotEvent('Tone',.5,color=toneColor);
    legend(legendList,'Location','northeast');
end
% saveFigures(gcf,['Summary_licking_',taskRange{task}],...
%     strcat(resultspath),...
%     saveFIG=true,savePDF=true);

%% Plot clamp trace aligned to Tone during omission

timeRange = [-0.5,3];
totalTrialRange = 'All';
trialRange = 'All';

initializeFig(.3,.5); tiledlayout('flow');
nexttile;
combined = combineTraces(animals,timeRange=timeRange,...
                            eventRange='Tone (clamp omission)',...
                            taskRange='Reward-Clamp-Omission',...
                            animalRange=animalRange{1},...
                            totalTrialRange=totalTrialRange,...
                            trialRange=trialRange,...
                            signalRange='redClamp');
legendList = plotGroupTraces(combined.data{1},combined.timestamp,sessionColormap_clamped,...
                groupby='session',startIdx=combined.options.startIdx);

title('Red clamp command');
xlabel('Time (s)'); ylabel('clamp command (%)'); %ylim([0,Inf]);
plotEvent('Tone',.5,color=toneColor);
legend(legendList,'Location','northeast');

% nexttile;
% combined = combineTraces(animals,timeRange=timeRange,...
%                             eventRange='Tone (clamp omission)',...
%                             taskRange='Reward-Clamp-Omission',...
%                             animalRange=animalRange{1},...
%                             totalTrialRange=totalTrialRange,...
%                             trialRange=trialRange,...
%                             signalRange='blueClamp');
% legendList = plotGroupTraces(combined.data{1},combined.timestamp,sessionColormap_clamped,...
%                 groupby='session',startIdx=combined.options.startIdx);
% 
% title('Blue clamp command');
% xlabel('Time (s)'); ylabel('clamp command (%)'); %ylim([0,Inf]);
% plotEvent('Tone',.5,color=toneColor);
% legend(legendList,'Location','northeast');

% saveFigures(gcf,['Summary_licking_',taskRange{task}],...
%     strcat(resultspath),...
%     saveFIG=true,savePDF=true);


%% Plot DA trace aligned to Tone back to reward

timeRange = [-0.5,3];
signalRange = {'NAc-clamp','NAc-right'};

eventRange = {'Tone (unclamp)','Tone (unclamp)'};
taskRange = {'Reward-Unclamp-BringBackReward','Reward-Unclamp-BringBackReward'};

totalTrialRange = 'All';
trialRange = 'All';
nGroupsList = 4;

initializeFig(.6,.5); tiledlayout('flow');
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
                    groupSize=20,...
                    groupby='trials',startIdx=combined.options.startIdx);

    title(animalTypes{type});
    xlabel('Time (s)'); ylabel('dff'); ylim([-0.05,0.45]);
    plotEvent('Tone',.5,color=toneColor);
    legend(legendList,'Location','northeast');
end
% saveFigures(gcf,['Summary_licking_',taskRange{task}],...
%     strcat(resultspath),...
%     saveFIG=true,savePDF=true);

%% Plot lick trace

timeRange = [-0.5,3];

totalTrialRange = 'All';
trialRange = 'All';

initializeFig(.6,.5); tiledlayout('flow');

% Plot clamp animals
nexttile;
combined = combineTraces(animals,timeRange=timeRange,...
                            eventRange='Tone (clamp omission)',...
                            taskRange='Reward-Clamp-Omission',...
                            animalRange=animalRange{1},...
                            totalTrialRange=totalTrialRange,...
                            trialRange=trialRange,...
                            signalRange='Lick');
legendList = plotGroupTraces(combined.data{1},combined.timestamp,sessionColormap_clamped,...
                groupSize=20,...
                groupby='sessions',startIdx=combined.options.startIdx);
combined = combineTraces(animals,timeRange=timeRange,...
                            eventRange='Tone (unclamp omission)',...
                            taskRange='Reward-Unclamp-Omission',...
                            animalRange=animalRange{1},...
                            totalTrialRange=totalTrialRange,...
                            trialRange=trialRange,...
                            signalRange='Lick');
plotSEM(combined.timestamp,combined.data{1},ctrlColor);
legendList{end+1} = strcat("Unclamp (n=",num2str(size(combined.data{1},1)),")");
xlabel('Time (s)'); ylabel('Licks/s'); ylim([0 3]);
plotEvent('Tone',.5,color=toneColor);
title(animalTypes{1});
legend(legendList,'Location','northeast');

% Plot ctrl animals
nexttile;
combined = combineTraces(animals,timeRange=timeRange,...
                            eventRange='Tone (unclamp omission)',...
                            taskRange='Reward-Unclamp-Omission',...
                            animalRange=animalRange{2},...
                            totalTrialRange=totalTrialRange,...
                            trialRange=trialRange,...
                            signalRange='Lick');
legendList = plotGroupTraces(combined.data{1},combined.timestamp,bluePurpleRed,...
                groupby='session',startIdx=combined.options.startIdx);
xlabel('Time (s)'); ylabel('Licks/s'); ylim([0 3]);
plotEvent('Tone',.5,color=toneColor);
title(animalTypes{2});
legend(legendList,'Location','northeast');

% saveFigures(gcf,['Summary_licking_',taskRange{task}],...
%     strcat(resultspath),...
%     saveFIG=true,savePDF=true);

%% Plot grouped anticipatory lick changes

panelTask = ["Reward-Clamp-Omission",         "Reward-Unclamp-Omission"; ...
             "Reward-Unclamp-Omission",       "Reward-Unclamp-Omission"; ...
             "Reward-Unclamp-BringBackReward","Reward-Unclamp-BringBackReward"];
% Columns are clamp and ctrl animals; session indices are within each task row.
panelSessions = {1:3, 1:3; 1, 4; 1, 1};
panelTitle = {'First 3 sessions','4th session','5th session: bring back reward'};
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

groupedResults = cell(size(panelTask));
tileSpans = ones(1,size(panelTask,1));
for panel = 1:size(panelTask,1)
    for type = 1:2
        curStats = cell(numel(animalRange{type}),numel(eventRange));

        for animal = 1:numel(animalRange{type})
            rowIdx = find(strcmpi({animals.animal},animalRange{type}{animal}) & ...
                          strcmpi({animals.task},char(panelTask(panel,type))) & ...
                          contains({animals.event},eventRange{1},IgnoreCase=true) & ...
                          contains({animals.name},signalRange,IgnoreCase=true));
            [~,sortIdx] = sort(arrayfun(@(x) min(x.options.signalRows),animals(rowIdx)));

            for row = rowIdx(sortIdx)
                sessionStart = animals(row).options.startIdx.session{1};
                sessionStop = [sessionStart(2:end)-1; numel(animals(row).trialInfo.trialNumber)];
                for session = panelSessions{panel,type}
                    if session > numel(sessionStart); continue; end
                    trialIdx = (sessionStart(session):sessionStop(session))';
                    if ~(ischar(conditionRange) || isstring(conditionRange))
                        trialNumber = animals(row).trialInfo.trialNumber(trialIdx);
                        trialIdx = trialIdx(trialNumber >= conditionRange(1) & trialNumber <= conditionRange(2));
                    end
                    if ~isempty(trialConditions)
                        trials = animals(row).trialInfo.trialTable;
                        trialIdx = intersect(trialIdx,find(eval(trialConditions)));
                    end
                    curStats{animal,1} = [curStats{animal,1}; animals(row).trialInfo.trialTable.(statsType)(trialIdx)];
                end
            end
        end

        combinedStats.stats = {curStats};
        combinedStats.options.eventRange = eventRange;
        combinedStats.options.inTrialTable = true;
        groupedResults{panel,type} = plotGroupedTrialStats(combinedStats,ylabelList,groupSize=groupSize,...
                                    plotIndividual=false, plotCommonTrials=false, plot=false);
        if ~isempty(groupedResults{panel,type}.traces{1})
            tileSpans(panel) = max(tileSpans(panel),size(groupedResults{panel,type}.traces{1}{1},2));
        end
    end
end

initializeFig(.45,.6);
t = tiledlayout(1,sum(tileSpans),TileSpacing='compact');
ax = gobjects(1,size(panelTask,1));
groupLabel = {'Clamp','Ctrl'};
groupColor = {clampColor,unclampColor};
lineStyle = ["-","--"];
tileStart = 1;
for panel = 1:size(panelTask,1)
    ax(panel) = nexttile(t,tileStart,[1 tileSpans(panel)]); hold on
    tileStart = tileStart + tileSpans(panel);

    for type = 1:2
        if isempty(groupedResults{panel,type}.traces{1}); continue; end
        grouped = groupedResults{panel,type}.traces{1}{1};
        x = (0:size(grouped,2)-1)*groupSize;
        scatter(x,mean(grouped,1,'omitnan'),200,groupColor{type},'filled',HandleVisibility='off');
        plotSEM(x,grouped,groupColor{type},label=groupLabel{type},LineStyle=lineStyle(type),plotIndividual=false);
    end

    title(panelTitle{panel}); xlabel('Trials'); xlim([0,max(1,tileSpans(panel))*groupSize]);
    if panel == 1; ylabel(ylabelList); else; yticklabels([]); end
    legend('Location','northeast');
end
linkaxes(ax,'y'); ylim([0 Inf]);

% saveFigures(gcf,'Summary_CSvsTrialsGrouped_Lick',...
%         strcat(resultspath),...
%         saveFIG=true,savePDF=true);


%% Bar plot of anticipatory licking (unfinished)

%% Plot grouped cue-evoked DA changes (use stageAmp)

panelTask = ["Reward-Clamp-Omission",         "Reward-Unclamp-Omission"; ...
             "Reward-Unclamp-Omission",       "Reward-Unclamp-Omission"; ...
             "Reward-Unclamp-BringBackReward","Reward-Unclamp-BringBackReward"];
% Columns are clamp and ctrl animals; session indices are within each task row.
panelSessions = {1:3, 1:3; 1, 4; 1, 1};
panelTitle = {'First 3 sessions','4th session','5th session: bring back reward'};
statsType = 'stageAmp';

conditionRange = 'All';
signalRange = {'NAc-clamp','NAc-right'};
trialConditions = 'trials.performing';
eventRange = {'Tone'};
groupSize = 20; % numbers of trials to calculate average

ylabelList = 'Cue-evoked DA amplitude';

groupedResults = cell(size(panelTask));
tileSpans = ones(1,size(panelTask,1));
for panel = 1:size(panelTask,1)
    for type = 1:2
        curStats = cell(numel(animalRange{type}),numel(eventRange));

        for animal = 1:numel(animalRange{type})
            rowIdx = find(strcmpi({animals.animal},animalRange{type}{animal}) & ...
                          strcmpi({animals.task},char(panelTask(panel,type))) & ...
                          contains({animals.event},eventRange{1},IgnoreCase=true) & ...
                          contains({animals.name},signalRange{type},IgnoreCase=true));
            [~,sortIdx] = sort(arrayfun(@(x) min(x.options.signalRows),animals(rowIdx)));

            for row = rowIdx(sortIdx)
                sessionStart = animals(row).options.startIdx.session{1};
                sessionStop = [sessionStart(2:end)-1; numel(animals(row).trialInfo.trialNumber)];
                for session = panelSessions{panel,type}
                    if session > numel(sessionStart); continue; end
                    trialIdx = (sessionStart(session):sessionStop(session))';
                    if ~(ischar(conditionRange) || isstring(conditionRange))
                        trialNumber = animals(row).trialInfo.trialNumber(trialIdx);
                        trialIdx = trialIdx(trialNumber >= conditionRange(1) & trialNumber <= conditionRange(2));
                    end
                    if ~isempty(trialConditions)
                        trials = animals(row).trialInfo.trialTable;
                        trialIdx = intersect(trialIdx,find(eval(trialConditions)));
                    end
                    curStats{animal,1} = [curStats{animal,1}; animals(row).(statsType).data(trialIdx,:)];
                end
            end
        end

        combinedStats.stats = {curStats};
        combinedStats.options.eventRange = eventRange;
        combinedStats.options.inTrialTable = false;
        groupedResults{panel,type} = plotGroupedTrialStats(combinedStats,ylabelList,groupSize=groupSize,...
                                    plotIndividual=false, plotCommonTrials=false, plot=false);
        if ~isempty(groupedResults{panel,type}.traces{1})
            tileSpans(panel) = max(tileSpans(panel),size(groupedResults{panel,type}.traces{1}{1},2));
        end
    end
end

initializeFig(.45,.6);
t = tiledlayout(1,sum(tileSpans),TileSpacing='compact');
ax = gobjects(1,size(panelTask,1));
groupLabel = {'Clamp','Ctrl'};
groupColor = {clampColor,unclampColor};
lineStyle = ["-","--"];
tileStart = 1;
for panel = 1:size(panelTask,1)
    ax(panel) = nexttile(t,tileStart,[1 tileSpans(panel)]); hold on
    tileStart = tileStart + tileSpans(panel);

    for type = 1:2
        if isempty(groupedResults{panel,type}.traces{1}); continue; end
        grouped = groupedResults{panel,type}.traces{1}{1};
        x = (0:size(grouped,2)-1)*groupSize;
        scatter(x,mean(grouped,1,'omitnan'),200,groupColor{type},'filled',HandleVisibility='off');
        plotSEM(x,grouped,groupColor{type},label=groupLabel{type},LineStyle=lineStyle(type),plotIndividual=false);
    end

    title(panelTitle{panel}); xlabel('Trials'); xlim([0,max(1,tileSpans(panel))*groupSize]);
    if panel == 1; ylabel(ylabelList); else; yticklabels([]); end
    legend('Location','northeast');
end
linkaxes(ax,'y'); %ylim([0 Inf]);

% saveFigures(gcf,'Summary_CSvsTrialsGrouped_stageAmp',...
%         strcat(resultspath),...
%         saveFIG=true,savePDF=true);

%% Bar plot of DA (unfinished)

%% BringBackReward analysis