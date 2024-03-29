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
answer = questdlg('Group sessions or load summary?','Select load sources',...
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

%% Optional: Make changes to summary for further analysis

for i = 1:211; summary(i).task = 'Random'; end
for i = 212:499; summary(i).task = 'Reward1'; end
for i = 500:919; summary(i).task = 'Punish1'; end
for i = 920:1279; summary(i).task = 'Reward2'; end
for i = 1280:1594; summary(i).task = 'Punish2'; end

%% Create animals struct

if isempty(dir(fullfile(resultspath,'animals*.mat'))) || groupSessions
    animals = getAnimalsStruct(summary);
end

%% Optional: update animals struct

for i = 1:length(animals)
   animals(i).options.startIdx.animal =  animals(i).options.animalStartIdx;
   animals(i).options.startIdx.session = animals(i).options.sessionStartIdx;

   animals(i).options = rmfield(animals(i).options,'animalStartIdx');
   animals(i).options = rmfield(animals(i).options,'sessionStartIdx');
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

initializeFig(0.5,0.3);
combined = combineTraces(animals,timeRange=[-0.5,3],...
                            eventRange='Water',...
                            animalRange="SL133",...
                            taskRange='Random',...
                            totalTrialRange='All',...
                            trialRange='All',...
                            signalRange='NAc');
plotSEM(combined.timestamp,shuffleTraces(combined.data{1}),[0.75,0.75,0.75]);
plotSEM(combined.timestamp,combined.data{1},bluePurpleRed(1,:));
plotEvent('Water',0,color=bluePurpleRed(1,:))
xlabel('Time (s)'); ylabel('z-score');

%% Baseline: plot water, airpuff, stim, tone

timeRange = [-0.5,3];
eventRange = {'Water','Airpuff','Stim','Tone'};
animalRange = 'All';
taskRange = 'Random';
trialRange = 'All'; % range of trials in each session
totalTrialRange = 'All';
signalRange = 'NAc';

colorList = {bluePurpleRed(1,:),[.2,.2,.2],bluePurpleRed(500,:),bluePurpleRed(350,:)};
eventDuration = [0,.1,.5,.5];

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
    plotSEM(combined.timestamp,shuffleTraces(combined.data{1}),[.75,.75,.75]);
    plotSEM(combined.timestamp,combined.data{1},colorList{i});
    xlabel('Time (s)'); ylabel([signalRange,' z-score']); ylim([-1,3]);
    plotEvent(eventRange{i},eventDuration(i),color=colorList{i});
    legend({['Shuffled (n=',num2str(size(combined.data{1},1)),')'],...
            [eventRange{i},' (n=',num2str(size(combined.data{1},1)),')']},...
            'Location','northeast');
end
saveFigures(gcf,['Summary_random_',signalRange],...
        strcat(resultspath),...
        saveFIG=true,savePDF=true);

%% Plot overall to show animal learned

timeRange = [-0.5,3];
eventRange = {'Pair'};
animalRange = {'SL133','SL135','SL136'};%'All';
taskRange = {'Reward1','Punish1','Reward2','Punish2'};
totalTrialRange = 'All';
trialRange = 'All';
signalRange = 'NAc';

groupSizeList = [50;50];
nGroupsList = [15;15];

for task = 1:length(taskRange)
    initializeFig(0.5,0.3); tiledlayout('flow');
    for event = 1:length(eventRange)
        nexttile;
        combined = combineTraces(animals,timeRange=timeRange,...
                                    eventRange=eventRange{event},...
                                    animalRange=animalRange,...
                                    taskRange=taskRange{task},...
                                    totalTrialRange=totalTrialRange,...
                                    trialRange=trialRange,...
                                    signalRange=signalRange,...
                                    trialConditions="strcmpi(trials.Outcome,'H') & trials.nAnticipatoryLicks >= 3");
        legendList = plotGroupTraces(combined.data{1},combined.timestamp,bluePurpleRed,...
                        groupSize=groupSizeList(event),nGroups=nGroupsList(event),...
                        groupby='trials',startIdx=combined.options.startIdx,remaining='separate');
        plotEvent(eventRange{event},.5,color=bluePurpleRed(500,:));
        xlabel('Time (s)'); ylabel([signalRange,' z-score']);
        legend(legendList,'Location','northeast');
    end
    % saveFigures(gcf,['Summary_pairing_',taskRange{task},'-',eventRange{event},'_',signalRange],...
    %     strcat(resultspath),...
    %     saveFIG=true,savePDF=true);
end
% autoArrangeFigures

%% Plot overall to show animal learned (w or w/o paAIP2)

timeRange = [-0.5,3];
eventRange = {'Pair'};
animalRange = {'SL133','SL135','SL136'};%'All';
taskRange = {'Reward1','Punish1','Reward2','Punish2'};
totalTrialRange = [1,60;60,150];
trialRange = 'All';
signalRange = 'NAc';

groupSizeList = [20,20;20,20];
nGroupsList = [10,10;50,50];
eventColor = {bluePurpleRed(1,:),bluePurpleRed(400,:)};

for task = 1:1%length(taskRange)
    initializeFig(0.5,0.3); tiledlayout('flow');
    for t = 1:size(totalTrialRange,1)
        for event = 1:length(eventRange)
            nexttile; 
            combined = combineTraces(animals,timeRange=timeRange,...
                                        eventRange=eventRange{event},...
                                        animalRange=animalRange,...
                                        taskRange=taskRange{task},...
                                        totalTrialRange=totalTrialRange(t,:),...
                                        trialRange=trialRange,...
                                        signalRange=signalRange);
            legendList = plotGroupTraces(combined.data{1},combined.timestamp,bluePurpleRed,...
                            groupSize=groupSizeList(t,event),nGroups=nGroupsList(t,event),...
                            groupby='trials',startIdx=combined.options.startIdx);
            plotEvent(eventRange{event},.5,color=eventColor{t});
            xlabel('Time (s)'); ylabel([signalRange,' z-score']);
            legend(legendList,'Location','northeast');
        end
    end
    % saveFigures(gcf,['Summary_paAIP2_',taskRange{task},'-',eventRange{event},'_',signalRange],...
    %     strcat(resultspath),...
    %     saveFIG=true,savePDF=true);
end
% autoArrangeFigures

%% Plot lick trace

timeRange = [-0.5,3];
eventRange = {'Stim','Pair'};
animalRange = {'SL133','SL135','SL136'};%'All';
taskRange = {'Reward1','Punish1','Reward2','Punish2'};
totalTrialRange = 'All';
trialRange = 'All';
signalRange = 'Lick';

groupSizeList = [50;50];
nGroupsList = [15;15];

for task = 1:length(taskRange)
    initializeFig(.5,.5); tiledlayout('flow');
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
        plotEvent(eventRange{event},.5,color=bluePurpleRed(500,:));
        xlabel('Time (s)'); ylabel('Licks/s'); ylim([0 inf]);
        legend(legendList,'Location','northeast');
    end
    saveFigures(gcf,['Summary_licking_',taskRange{task},'-',eventRange{event},'_',signalRange],...
        strcat(resultspath),...
        saveFIG=true,savePDF=true);
end
% autoArrangeFigures

%% Plot anticipatory lick changes

eventRange = {'Stim','Pair'};
animalRange = {'SL133','SL135','SL136'};%,'SL137'};%'All';
taskRange = {'Reward1','Punish1','Reward2','Punish2'};
conditionRange = [1,60;61,150];
signalRange = 'NAc';
conditionColors = {bluePurpleRed(1,:),[.213 .543 .324]};

% Get subtrial stats
if strcmpi(animalRange,'All'); animalRange = unique({animals.animal}); end
for task = 1:length(taskRange)
    stats_combined = cell(length(animalRange),length(eventRange));
    % animalList_task = unique({animals(strcmpi({animals.task},taskRange{task})).animal});
    % animalList = intersect(animalList_task,animalRange);
    for event = 1:length(eventRange)
        for animal = 1:length(animalRange)
            for t = 1:size(conditionRange,1)
                combined = combineTraces(animals,statsType=statsType,...
                                            eventRange=eventRange{event},...
                                            animalRange=animalRange{animal},...
                                            taskRange=taskRange{task},...
                                            totalTrialRange=conditionRange(t,:),...
                                            signalRange=signalRange);
                if combined.options.empty; continue; end
                nal = combined.trialTable{1}.nAnticipatoryLicks;
                sessionStartIdx = combined.options.startIdx.session{1};
    
                sessionStats = cell(length(sessionStartIdx),1);
    
                for session = 1:length(sessionStartIdx)
                    if session == length(sessionStartIdx); lastTrial = length(combined.trialNumber{1}); 
                    else; lastTrial = sessionStartIdx(session+1)-1; end
                    sessionWindow = sessionStartIdx(session):lastTrial;
    
                    sessionStats{session} = [sessionStats{session}; nal(sessionWindow)];
                end
                stats_combined{animal,event} = [stats_combined{animal,event}, sessionStats];
            end
        end 
    end
    
    % Plot scatter plot and best fit line
    initializeFig(.7,.7); tiledlayout('flow');
    for event = 1:length(eventRange)
        for animal = 1:length(animalRange)
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
            title([taskRange{task}, ': ',animalRange{animal},' -> ',eventRange{event}]);
        end
    end
    saveFigures(gcf,['Summary_AnticipatoryLicksvsTrials_',taskRange{task},'-',eventRange{event},'_',signalRange],...
        strcat(resultspath),...
        saveFIG=true,savePDF=true);
end
% autoArrangeFigures

%% Plot CS DA response vs trials
eventRange = {'Stim','Pair'};
animalRange = {'SL133','SL135','SL136'};%,'SL137'};%'All';
taskRange = {'Reward1','Punish1','Reward2','Punish2'};
conditionRange = [1,60;61,150];
signalRange = 'NAc';
conditionColors = {bluePurpleRed(1,:),[.213 .543 .324]};

stage = 2; % Plot CS only
statsType = 'stageAvg';

% Get subtrial stats
for task = 1:length(taskRange)
    stats_combined = cell(length(animalRange),length(eventRange));
    for event = 1:length(eventRange)
        for animal = 1:length(animalRange)
            for t = 1:size(conditionRange,1)
                combined = combineTraces(animals,statsType=statsType,...
                                            eventRange=eventRange{event},...
                                            animalRange=animalRange{animal},...
                                            taskRange=taskRange{task},...
                                            totalTrialRange=conditionRange(t,:),...
                                            signalRange=signalRange);
                if combined.options.empty; continue; end
                statsData = combined.stats.(statsType){1};
                sessionStartIdx = combined.options.startIdx.session{1};
    
                sessionStats = cell(length(sessionStartIdx),1);
    
                for session = 1:length(sessionStartIdx)
                    if session == length(sessionStartIdx); lastTrial = length(combined.trialNumber{1}); 
                    else; lastTrial = sessionStartIdx(session+1)-1; end
                    sessionWindow = sessionStartIdx(session):lastTrial;
    
                    sessionStats{session} = [sessionStats{session}; statsData(sessionWindow,:)];
                end
                stats_combined{animal,event} = [stats_combined{animal,event}, sessionStats];
            end
        end 
    end
    
    % Plot scatter plot and best fit line
    initializeFig(.7,.7); tiledlayout('flow');
    for event = 1:length(eventRange)
        for animal = 1:length(animalRange)
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
            title([taskRange{task}, ': ',animalRange{animal},' -> ',eventRange{event}]);
        end
    end
    saveFigures(gcf,['Summary_CSvsTrials_',taskRange{task},'-',eventRange{event},'_',signalRange],...
        strcat(resultspath),...
        saveFIG=true,savePDF=true);
end
% autoArrangeFigures

%% Plot slope (bar plot) for NAc during/after paAIP2

eventRange = {'Stim','Pair'};
animalRange = {'SL133','SL135','SL136'};%,'SL137'};%'All';
taskRange = {'Reward1','Punish1','Reward2','Punish2'};
conditionRange = [1,60;61,120];
signalRange = 'NAc';
conditionLabels = {'paAIP2','Control'};
conditionColors = {bluePurpleRed(1,:),[.213 .543 .324]};
statsType = 'stageAvg';
stageLabels = {'Baseline','CS','US'};
stageColors = {[.75 .75 .75],bluePurpleRed(end,:),bluePurpleRed(1,:)};
animalColors = {bluePurpleRed(50,:),bluePurpleRed(200,:),bluePurpleRed(350,:),bluePurpleRed(500,:)};

% get subtrial stats
for task = 1:length(taskRange)
    stageFit = cell(length(eventRange),size(conditionRange,1),length(animalRange),length(stageColors));
    for event = 1:length(eventRange)
        for t = 1:size(conditionRange,1)
            for animal = 1:length(animalRange)
                combined = combineTraces(animals,statsType=statsType,...
                                            eventRange=eventRange{event},...
                                            animalRange=animalRange{animal},...
                                            taskRange=taskRange{task},...
                                            totalTrialRange=conditionRange(t,:),...
                                            signalRange=signalRange);
                if combined.options.empty; continue; end
                statsData = combined.stats.(statsType){1};
                sessionStartIdx = combined.options.startIdx.session{1};
        
                % Fit stageAvg across session
                for stage = 1:size(statsData,2)
                    sessionFit = nan(length(sessionStartIdx),2);
                    for session = 1:length(sessionStartIdx)
                        if session == length(sessionStartIdx); lastTrial = length(combined.trialNumber{1}); 
                        else; lastTrial = sessionStartIdx(session+1)-1; end
                        sessionWindow = sessionStartIdx(session):lastTrial;
                        y = statsData(sessionWindow,stage)';
                        sessionFit(session,:) = polyfit(1:length(sessionWindow),y,1);
                    end
                    stageFit{event,t,animal,stage} = sessionFit;
                end
            end
        end
    end
    
    % Plot bar plot
    fitLabels = {'Slope','Intercept'};
    for i = 1:1%length(fitLabels)
        initializeFig(.5,.5); tiledlayout('flow');
        for event = 1:length(eventRange)
            for stage = 1:size(statsData,2)
                nexttile;
                allSessions = cell(size(conditionRange,1),1);
                for t = 1:size(conditionRange,1)
                    % Plot bar plot
                    allSessions{t} = cell2mat({stageFit{event,t,:,stage}}');
                    boxchart(ones(size(allSessions{t},1),1)*t,allSessions{t}(:,i),...
                            'BoxFaceColor',conditionColors{t}); hold on
                    for animal = 1:length(animalRange)
                        animalSessions = stageFit{event,t,animal,stage};
                        if isempty(animalSessions); continue; end
                        swarmchart(ones(size(animalSessions,1))*t,animalSessions(:,i),...
                                    100,animalColors{animal},'filled',...
                                    'MarkerFaceAlpha',0.8,'XJitter','density','XJitterWidth',0.01); hold on
                    end
                end
                for session = 1:size(allSessions{1},1)
                    plot([1,2],[allSessions{1}(session,i),allSessions{2}(session,i)],color=[.75,.75,.75]); hold on            
                end
                [~,p,~] = kstest2(allSessions{1}(:,i),allSessions{2}(:,i)); % plotSignificance(p,[1 2],0.5);
                ax = gca; ax.Children = [ax.Children(size(allSessions{1},1)+1:end); ax.Children(1:size(allSessions{1},1))];
                xticks([1 2]); xticklabels(conditionLabels);
                ylabel(['Subtrial ',fitLabels{i},' (',signalRange,')']); % ylim([0,20]);
                title([eventRange{event},' -> ',stageLabels{stage},' (p=',num2str(p),')']);
            end
        end
    end
    saveFigures(gcf,['Summary_paAIP2Slopes_',taskRange{task},'-',eventRange{event},'_',signalRange],...
        strcat(resultspath),...
        saveFIG=true,savePDF=true);
end
% autoArrangeFigures