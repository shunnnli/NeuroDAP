% Shun_analyzeOptoPair
% 2023/04/25

% Outputs a .mat file with aligned traces for each animal and behavior

%% Setup

clear; close all;
addpath(genpath(osPathSwitch('/Volumes/MICROSCOPE/Shun/Analysis/NeuroDAP/Methods')));
[~,~,~,~,~,~,bluePurpleRed] = loadColors;
r2p_cmap = getColormap([255, 50, 58],[0 0 0],500,'midcol',[255 255 255]);
p2r_cmap = getColormap([241 160 255],[0 0 0],500,'midcol',[255 255 255]);
today = char(datetime('today','Format','yyyyMMdd'));

% Define result directory
resultspath = osPathSwitch('/Volumes/MICROSCOPE/Shun/Project valence/Results');

% Building summary struct from selected sessions
answer = questdlg('Group sessions or load summary?','Select load sources',...
                  'Group single sessions','Load combined struct','Load combined struct');

if strcmpi(answer,'Group single sessions')
    sessionList = uipickfiles('FilterSpec',osPathSwitch('/Volumes/MICROSCOPE/Shun/Project valence/Recordings'));
    groupSessions = true;
    % Update resultspath
    dirsplit = strsplit(sessionList{1},filesep); projectName = dirsplit{end-1}; 
    resultspath = strcat(resultspath,filesep,projectName);
    % Create resultspath if necessary
    if isempty(dir(resultspath)); mkdir(resultspath); end

elseif strcmpi(answer,'Load combined struct')
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
end

%% Create/modify summary struct (only need to do this for initial loading)

regroup = true;
if groupSessions
    % Create/modify summary struct
    if isempty(dir(fullfile(resultspath,'summary*.mat'))) || regroup
        summary = struct([]);
        disp('Finished: summary.mat not found, created a new one');
    end
    
    % Group sessions to summary
    for s = 1:length(sessionList)
        % Find sessionName
        dirsplit = strsplit(sessionList{s},filesep); 
        sessionName = dirsplit{end}; clear dirsplit
    
        % Load session
        load(strcat(sessionList{s},filesep,'analysis_',sessionName,'.mat'));
    
        summary = [summary,analysis];
        disp(['Finished: session ', sessionName,' loaded (',...
              num2str(s),'/',num2str(length(sessionList)),')']);
    end
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

%% Make changes to summary

for i = 1:211; summary(i).task = 'Random'; end
for i = 212:499; summary(i).task = 'Reward1'; end
for i = 500:919; summary(i).task = 'Punish1'; end
for i = 920:1279; summary(i).task = 'Reward2'; end
for i = 1280:1594; summary(i).task = 'Punish2'; end

%% Create animals struct

regroupAnimals = true;

if isempty(dir(fullfile(resultspath,'animals*.mat'))) || regroupAnimals
    animals = getAnimalsStruct(summary);
end


%% Save animals and summary struct (time consumming!!!)

% Save animals.mat
disp(['Ongoing: saving animals.mat (',char(datetime('now','Format','HH:mm:ss')),')']);
save(strcat(resultspath,filesep,'animals_',today),'animals','sessionList','-v7.3');
disp(['Finished: saved animals.mat (',char(datetime('now','Format','HH:mm:ss')),')']);

% Save summary.mat (not recommend!! Will take forever!!)
% disp(['Ongoing: saving summary.mat (',char(datetime('now','Format','HH:mm:ss')),')']);
% save(strcat(resultspath,filesep,'summary_',today),'summary','sessionList','-v7.3');
% disp(['Finished: saved summary.mat (',char(datetime('now','Format','HH:mm:ss')),')']);

%% Test: Plot traces from summary/animals struct

timeRange = [-0.5,3];
eventRange = 'Water';
animalRange = ["SL133","SL135"];
taskRange = 'Random';
totalTrialRange = 'All';
trialRange = 'All'; % range of trials in each session
signalRange = 'NAc';

combined = combineTraces(animals,timeRange=timeRange,...
                            eventRange=eventRange,...
                            animalRange=animalRange,...
                            taskRange=taskRange,...
                            totalTrialRange=totalTrialRange,...
                            trialRange=trialRange,...
                            signalRange=signalRange);
plotSEM(combined.timestamp,combined.data{1},[.213 .543 .324]);

%% Baseline: plot water, airpuff, stim, tone

timeRange = [-0.5,3];
eventRange = {'Water','Airpuff','Stim','Tone'};
animalRange = 'All';
taskRange = 'Random';
trialRange = 'All'; % range of trials in each session
totalTrialRange = 'All';
signalRange = 'EyeArea';

colorList = {bluePurpleRed(1,:),[.2,.2,.2],bluePurpleRed(500,:),bluePurpleRed(100,:)};
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
    plotSEM(combined.timestamp,combined.data{1},colorList{i});
    xlabel('Time (s)'); ylabel([signalRange,' z-score']); ylim([-8,1]);
    plotEvent(eventRange{i},eventDuration(i),color=colorList{i});
    legend({['Shuffled (n=',num2str(size(combined.data{1},1)),')'],...
            [eventRange{i},' (n=',num2str(size(combined.data{1},1)),')']},...
            'Location','northeast');
end
% saveFigures(gcf,['Summary_random_',signalRange],...
%         strcat(resultspath),...
%         saveFIG=true,savePDF=true);

%% Plot overall to show animal learned

timeRange = [-0.5,3];
eventRange = {'Stim','Pair','Tone'};
animalRange = {'SL133','SL135','SL136'};%'All';
taskRange = {'Reward1','Punish1','Reward2','Punish2'};
totalTrialRange = 'All';
trialRange = 'All';
signalRange = 'NAc';

colorList = {bluePurpleRed(500,:),bluePurpleRed(300,:),bluePurpleRed(100,:)};
groupSizeList = [50;50;10];
nGroupsList = [15;15;15];
ylimList = [-1,5; -1,4.5;-1,4; -1,2.5];

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
        % ylim(ylimList(task,:));
        plotEvent(eventRange{event},.5,color=colorList{event});
        xlabel('Time (s)'); ylabel([signalRange,' z-score']);
        legend(legendList,'Location','northeast');
    end
    % saveFigures(gcf,['Summary_pairing_',taskRange{task},'_',signalRange],...
    %     strcat(resultspath),...
    %     saveFIG=true,savePDF=true);
end
% autoArrangeFigures

%% Plot overall to show animal learned (w or w/o paAIP2)

timeRange = [-0.5,3];
eventRange = {'Stim','Pair'};
animalRange = {'SL133','SL135','SL136'};%'All';
taskRange = {'Reward1','Punish1','Reward2','Punish2'};
totalTrialRange = [1,60;60,150];
trialRange = 'All';
signalRange = 'NAc';

groupSizeList = [20,20;20,20];
nGroupsList = [10,10;50,50];
eventColor = {bluePurpleRed(1,:),bluePurpleRed(400,:)};

for task = 1:length(taskRange)
    initializeFig(.5,.5); tiledlayout('flow');
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
                            groupby='trials',startIdx=combined.options.startIdx,...
                            remaining='include');
            plotEvent(eventRange{event},.5,color=eventColor{t});
            xlabel('Time (s)'); ylabel([signalRange,' z-score']);
            legend(legendList,'Location','northeast');
        end
    end
    saveFigures(gcf,['Summary_paAIP2_',taskRange{task},'-',eventRange{event},'_',signalRange],...
        strcat(resultspath),...
        saveFIG=true,savePDF=true);
end
% autoArrangeFigures

%% Plot lick trace

timeRange = [-0.5,3];
eventRange = {'Stim','Pair','Tone'};
animalRange = {'SL133','SL135','SL136'};%'All';
taskRange = {'Reward1','Punish1','Reward2','Punish2'};
totalTrialRange = 'All';
trialRange = 'All';
signalRange = 'Lick';

colorList = {bluePurpleRed(500,:),bluePurpleRed(300,:),bluePurpleRed(100,:)};
groupSizeList = [50;50;10];
nGroupsList = [15;15;15];
ylimList = [-1,5; -1,4.5;-1,4; -1,2.5];

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
        plotEvent(eventRange{event},.5,color=colorList{event});
        xlabel('Time (s)'); ylabel('Licks/s'); ylim([0 inf]);
        legend(legendList,'Location','northeast');
    end
    saveFigures(gcf,['Summary_licking_',taskRange{task},'_',signalRange],...
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
                nal = combined.trialTable{1}.nAnticipatoryLicks;
                sessionStartIdx = combined.options.sessionStartIdx{1};
    
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

%% Plot grouped CS DA response (grouped across animal and 10 trials)

eventRange = {'Stim','Pair','Tone'};
animalRange = 'All';
taskRange =  {'Reward2','Punish2'}; %{'Reward1','Punish1'};
signalRange = 'NAc';

colorList = {bluePurpleRed(500,:),bluePurpleRed(300,:),bluePurpleRed(100,:)};
stage = 2; % Plot CS only
statsTypes = {'stageMax','stageMin','stageMax','stageMin'};
ylabelList = {'Max DA response during cue','Min DA response during cue',...
                'Max DA response during cue','Min DA response during cue'};

groupSize = 10; % numbers of trials to calculate average

combinedStats = getGroupedTrialStats(animals,statsTypes,...
                            eventRange=eventRange,...
                            animalRange=animalRange,...
                            taskRange=taskRange,...
                            signalRange=signalRange);

initializeFig(.7,.7); tiledlayout('flow');
plotGroupedTrialStats(combinedStats,ylabelList,stage=stage,groupSize=10,color=colorList);

% saveFigures(gcf,['Summary_CSvsTrialsGrouped_Second'],...
%         strcat(resultspath),...
%         saveFIG=true,savePDF=true);

%% Plot CS DA response vs trials
eventRange = {'Stim','Pair'};
animalRange = {'SL133','SL135','SL136'};%,'SL137'};%'All';
taskRange = {'Reward1','Punish1','Reward2','Punish2'};
conditionRange = [1,60;61,150];
signalRange = 'NAc';
conditionLabels = {'paAIP2','Control'};
conditionColors = {bluePurpleRed(1,:),[.213 .543 .324]};
stageLabels = {'Baseline','CS','US'};
stageColors = {[.75 .75 .75],bluePurpleRed(end,:),bluePurpleRed(1,:)};
animalColors = {bluePurpleRed(50,:),bluePurpleRed(200,:),bluePurpleRed(350,:),bluePurpleRed(500,:)};

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
                sessionStartIdx = combined.options.sessionStartIdx{1};
    
                sessionStats = cell(length(sessionStartIdx),1);
                totalTrialNum = cell(length(sessionStartIdx),1);
    
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
                sessionStartIdx = combined.options.sessionStartIdx{1};
        
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

%% Baseline -> reward: plot best fit line for lick during/after paAIP2

%% Baseline -> reward: plot slope for lick during/after paAIP2

%% Baseline -> reward: plot hit rate (bar) during/after paAIP2

%% Baseline -> reward: plot total lick per trial during/after paAIP2 (should be same)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Not used %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
