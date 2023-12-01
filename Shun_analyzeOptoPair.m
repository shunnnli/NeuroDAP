% Shun_analyzeOptoPair
% 2023/04/25

% Outputs a .mat file with aligned traces for each animal and behavior

%% Setup

clear; close all;
addpath(genpath(osPathSwitch('/Volumes/MICROSCOPE/Shun/Analysis/NeuroDAP/Methods')));
[~,~,~,blueGreenYellow,blueWhiteRed,~,bluePurpleRed,purpleWhiteRed] = loadColors;
r2p_cmap = getColormap([255, 50, 58],[0 0 0],500,'midcol',[255 255 255]);
p2r_cmap = getColormap([241 160 255],[0 0 0],500,'midcol',[255 255 255]);
today = char(datetime('today','Format','yyyyMMdd'));

% Define result directory
resultspath = osPathSwitch('/Volumes/MICROSCOPE/Shun/Project valence/Results');

% Building summary struct from selected sessions
answer = questdlg('Group sessions or load summary?','Select load sources',...
                  'Group single sessions','Load summary struct','Load summary struct');

if strcmpi(answer,'Group single sessions')
    sessionList = uipickfiles('FilterSpec',osPathSwitch('/Volumes/MICROSCOPE/Shun/Project valence/Recordings'));
    groupSessions = true;
    % Update resultspath
    dirsplit = strsplit(sessionList{1},filesep); projectName = dirsplit{end-1}; 
    resultspath = strcat(resultspath,filesep,projectName);
    % Create resultspath if necessary
    if isempty(dir(resultspath)); mkdir(resultspath); end

elseif strcmpi(answer,'Load summary struct')
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

%% Create/modify summary struct

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

%% Make changes to summary if needed

for i = 1:1279
    if isstring(summary(i).name) 
        summary(i).name = convertStringsToChars(summary(i).name);
    end
end

for i = 1:211
    summary(i).task = 'baseline';
end

for i = 212:499
    summary(i).task = 'baseline->reward';
end

for i = 500:919
    summary(i).task = 'reward->punish';
end

for i = 920:1279
    summary(i).task = 'punish->reward';
end

%% Create animals struct

regroupAnimals = true;

if isempty(dir(fullfile(resultspath,'animals*.mat'))) || regroupAnimals
    animals = struct([]);
    animalList = unique({summary.animal});
    disp('Finished: animal.mat not found, created a new one');

    for a = 1:length(animalList)
        animalIdx = find(cellfun(@(x) contains(x,animalList{a}), {summary.animal}));
        taskList = unique({summary(animalIdx).task});
        
        for task = 1:length(taskList)
            taskIdx = cellfun(@(x) strcmpi(x,taskList{task}), {summary(animalIdx).task});
            taskRows = animalIdx(taskIdx);
            eventList = unique({summary(taskRows).event});
    
            for event = 1:length(eventList)
                eventIdx = cellfun(@(x) contains(x,eventList{event},IgnoreCase=true), {summary(taskRows).event});
                eventRows = taskRows(eventIdx);
                signalList = unique({summary(eventRows).name});
    
                for signal = 1:length(signalList)
                    disp(['Ongoing: ',animalList{a},' -> ',taskList{task},' -> ',eventList{event},...
                        ' -> ',signalList{signal}]);
    
                    combined = combineTraces(summary,animalRange=animalList{a},...
                                    eventRange=eventList{event},...
                                    taskRange=taskList{task},...
                                    signalRange=signalList{signal},...
                                    statsType='All');
    
                    row = size(animals,2) + 1;
                    animals(row).animal = animalList{a};
                    animals(row).task = taskList{task};
                    animals(row).event = eventList{event};
                    animals(row).name = signalList{signal};
                    animals(row).system = combined.options.system;
                    animals(row).data = combined.data{1};
                    animals(row).stageAvg.data = combined.stats.stageAvg{1};
                    animals(row).stageMax.data = combined.stats.stageMax{1};
                    animals(row).stageMin.data = combined.stats.stageMin{1};
                    animals(row).timestamp = combined.timestamp;
                    animals(row).timeRange = combined.options.timeRange;
                    animals(row).finalFs = combined.options.finalFs;
                    animals(row).trialInfo.trialNumber = combined.trialNumber{1};
                    animals(row).trialInfo.trialTable = combined.trialTable{1};
                    animals(row).options = combined.options;
                end
            end
        end
    end
end


%% Save animals and summary struct (time consumming!!!)

% Save animals.mat
disp(['Ongoing: saving animals.mat (',char(datetime('now','Format','HH:mm:ss')),')']);
save(strcat(resultspath,filesep,'animals_',today),'animals','sessionList','-v7.3');
disp(['Finished: saved animals.mat (',char(datetime('now','Format','HH:mm:ss')),')']);

% Save summary.mat
disp(['Ongoing: saving summary.mat (',char(datetime('now','Format','HH:mm:ss')),')']);
save(strcat(resultspath,filesep,'summary_',today),'summary','sessionList','-v7.3');
disp(['Finished: saved summary.mat (',char(datetime('now','Format','HH:mm:ss')),')']);

%% Test: Plot traces from summary struct

timeRange = [-0.5,3];
eventRange = 'Water';
animalRange = ["SL133","SL135"];
taskRange = 'baseline';
totalTrialRange = 'All';
trialRange = 'All'; % range of trials in each session
signalRange = 'NAc';

combined = combineTraces(summary,timeRange=timeRange,...
                            eventRange=eventRange,...
                            animalRange=animalRange,...
                            taskRange=taskRange,...
                            totalTrialRange=totalTrialRange,...
                            trialRange=trialRange,...
                            signalRange=signalRange);

plotSEM(combined.timestamp,combined.data{1},[.213 .543 .324]);

%% Test: Plot traces from animals struct

timeRange = [-0.5,3];
eventRange = 'Water';
animalRange = ["SL133","SL135"];
taskRange = 'baseline';
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
taskRange = 'baseline';
trialRange = 'All'; % range of trials in each session
totalTrialRange = 'All';
signalRange = 'NAc';

colorList = {bluePurpleRed(1,:),[.2,.2,.2],bluePurpleRed(500,:),bluePurpleRed(350,:)};
eventDuration = [0,.1,.5,.5];

initializeFig(.5,.5); tiledlayout(2,2);
for i = 1:length(eventRange)
    nexttile;
    combined = combineTraces(summary,timeRange=timeRange,...
                                eventRange=eventRange{i},...
                                animalRange=animalRange,...
                                taskRange=taskRange,...
                                totalTrialRange=trialRange,...
                                trialRange=trialRange,...
                                signalRange=signalRange);
    plotSEM(combined.timestamp,shuffleTraces(combined.data{1}),[.75,.75,.75]);
    plotSEM(combined.timestamp,combined.data{1},colorList{i});
    plotEvent(eventRange{i},eventDuration(i),color=colorList{i});
    xlabel('Time (s)'); ylabel([signalRange,' z-score']);
    legend({['Shuffled (n=',num2str(size(combined.data{1},1)),')'],...
            [eventRange{i},' (n=',num2str(size(combined.data{1},1)),')']},...
            'Location','northeast');
end

%% Baseline -> reward: plot overall to show animal learned

timeRange = [-0.5,3];
eventRange = {'Stim','Pair'};
animalRange = {'SL133'};%,'SL135','SL136'};%'All';
taskRange = 'punish->reward';
totalTrialRange = 'All';
trialRange = {[1,10;11,20;21,30;31,40;41,50],...
              [1,20;21,40;41,60;61,80;81,100],...
              [1,10;11,20;21,30;31,40;41,50]};
signalRange = 'NAc';

colorList = [1,100,200,300,400,500];

initializeFig(.5,.5); tiledlayout('flow');
for event = 1:length(eventRange)
    nexttile; legendList = cell(size(trialRange,1),1);
    for t = 1:size(trialRange{event},1)
        combined = combineTraces(animals,timeRange=timeRange,...
                                    eventRange=eventRange{event},...
                                    animalRange=animalRange,...
                                    taskRange=taskRange,...
                                    totalTrialRange=totalTrialRange,...
                                    trialRange=trialRange{event}(t,:),...
                                    signalRange=signalRange);
        plotSEM(combined.timestamp,combined.data{1},bluePurpleRed(colorList(t),:));
        legendList{t} = ['Trial ', num2str(trialRange{event}(t,1)),'-',num2str(trialRange{event}(t,2)),...
                        ' (n=',num2str(size(combined.data{1},1)),')'];
    end
    plotEvent(eventRange{event},.5,color=bluePurpleRed(500,:));
    xlabel('Time (s)'); ylabel([signalRange,' z-score']);
    legend(legendList,'Location','northeast');
end

%% Baseline -> reward: plot overall to show animal learned (w or w/o paAIP2)

timeRange = [-0.5,3];
eventRange = {'Stim','Pair'};
animalRange = {'SL133','SL135','SL136'};%'All';
taskRange = 'baseline->reward';
totalTrialRange = [1,60;60,150];
trialRange = 'All';
signalRange = 'NAc';

groupSizeList = [10,20;10,20];
nGroupsList = [10,10;50,50];
eventColor = {blueWhiteRed(1,:),bluePurpleRed(350,:)};

initializeFig(.5,.5); tiledlayout('flow');
for t = 1:size(totalTrialRange,1)
    for event = 1:length(eventRange)
        nexttile; 
        combined = combineTraces(animals,timeRange=timeRange,...
                                    eventRange=eventRange{event},...
                                    animalRange=animalRange,...
                                    taskRange=taskRange,...
                                    totalTrialRange=totalTrialRange(t,:),...
                                    trialRange=trialRange,...
                                    signalRange=signalRange);
        legendList = plotGroupTraces(combined.data{1},combined.timestamp,bluePurpleRed,...
                        groupSize=groupSizeList(t,event),nGroups=nGroupsList(t,event),...
                        animalStartIdx=combined.options.animalStartIdx);
        plotEvent(eventRange{event},.5,color=eventColor{t});
        xlabel('Time (s)'); ylabel([signalRange,' z-score']);
        legend(legendList,'Location','northeast');
    end
end

%% Baseline -> reward: plot lick trace

timeRange = [-0.5,3];
eventRange = {'Stim','Pair'};
animalRange = {'SL133','SL135','SL136'};%'All';
taskRange = 'punish->reward';
trialRange = {[1,10;11,20;21,30;31,40;41,50],...
              [1,20;21,40;41,60;61,80;81,100],...
              [1,10;11,20;21,30;31,40;41,50]};
% trialRange = {[1,20;21,50;51,70;71,100;101,110],...
%               [1,40;41,100;101,140;141,200;201,240],...
%               [1,20;21,50;51,70;71,100;101,110]};
signalRange = 'lick';

colorList = [1,100,200,300,400,500];
eventDuration = [0,.1,.5,.5];

initializeFig(.5,.5); tiledlayout('flow');
for event = 1:length(eventRange)
    nexttile; legendList = cell(size(trialRange,1),1);
    for t = 1:size(trialRange{event},1)
        combined = combineTraces(animals,timeRange=timeRange,...
                                    eventRange=eventRange{event},...
                                    animalRange=animalRange,...
                                    taskRange=taskRange,...
                                    trialRange=trialRange{event}(t,:),...
                                    signalRange=signalRange);
        plotSEM(combined.timestamp,combined.data{1},bluePurpleRed(colorList(t),:));
        legendList{t} = ['Trial ', num2str(trialRange{event}(t,1)),'-',num2str(trialRange{event}(t,2)),...
                        ' (n=',num2str(size(combined.data{1},1)),')'];
    end
    plotEvent(eventRange{event},.5,color=bluePurpleRed(500,:));
    xlabel('Time (s)'); ylabel([signalRange,' z-score']);
    legend(legendList,'Location','northeast');
end

%% Baseline -> reward: plot anticipatory lick changes

timeRange = [-0.5,3];
eventRange = {'Stim','Pair'};
animalRange = {'SL133','SL135','SL136'};%'All';
taskRange = 'punish->reward';
trialRange = {[1,10;11,20;21,30;31,40;41,50],...
              [1,20;21,40;41,60;61,80;81,100],...
              [1,10;11,20;21,30;31,40;41,50]};
% trialRange = {[1,20;21,50;51,70;71,100;101,110],...
%               [1,40;41,100;101,140;141,200;201,240],...
%               [1,20;21,50;51,70;71,100;101,110]};
signalRange = 'lick';

colorList = [1,100,200,300,400,500];
eventDuration = [0,.1,.5,.5];

initializeFig(.5,.5); tiledlayout('flow');
for event = 1:length(eventRange)
    nexttile;
    for t = 1:size(trialRange{event},1)
        combined = combineTraces(animals,timeRange=timeRange,...
                                    eventRange=eventRange{event},...
                                    animalRange=animalRange,...
                                    taskRange=taskRange,...
                                    trialRange=trialRange{event}(t,:),...
                                    signalRange=signalRange);
        % al = combined.
        
    end
    plotEvent(eventRange{event},.5,color=bluePurpleRed(500,:));
    xlabel('Time (s)'); ylabel([signalRange,' z-score']);
    legend(legendList,'Location','northeast');
end

%% Baseline -> reward: plot stage scatter for NAc during/after paAIP2

%% Baseline -> reward: plot best fit line for NAc during/after paAIP2

%% Baseline -> reward: plot slope (bar plot) for NAc during/after paAIP2

eventRange = {'Stim','Pair'};
animalRange = {'SL133','SL135','SL136'};%,'SL137'};%'All';
taskRange = 'baseline->reward'; 
% taskRange = 'punish->reward';
conditionRange = [1,60;61,150];
signalRange = 'NAc';
conditionLabels = {'paAIP2','Control'};
conditionColors = {blueWhiteRed(50,:),[.75,.75,.75]};
statsType = 'stageAvg';
stageLabels = {'Baseline','CS','US'};
stageColors = {[.75 .75 .75],bluePurpleRed(end,:),bluePurpleRed(1,:)};
animalColors = {bluePurpleRed(50,:),bluePurpleRed(200,:),bluePurpleRed(350,:),bluePurpleRed(500,:)};

% get subtrial stats
stageFit = cell(length(eventRange),size(conditionRange,1),length(animalRange),length(stageColors));
for event = 1:length(eventRange)
    for t = 1:size(conditionRange,1)
        for animal = 1:length(animalRange)
            combined = combineTraces(summary,statsType=statsType,...
                                        eventRange=eventRange{event},...
                                        animalRange=animalRange{animal},...
                                        taskRange=taskRange,...
                                        totalTrialRange=conditionRange(t,:),...
                                        signalRange=signalRange);
            statsData = combined.stats.(statsType){1};
            trialWindowLength = conditionRange(t,2)-conditionRange(t,1) + 1;
            nSessions = sum(combined.options.signalRows);
    
            % Fit stageAvg across session
            for stage = 1:size(statsData,2)
                sessionFit = nan(nSessions,2);
                for session = 1:nSessions
                    sessionWindow = trialWindowLength*(session-1)+1 : min(trialWindowLength*session,size(statsData,1));
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
                    swarmchart(ones(size(animalSessions,1))*t,animalSessions(:,i),...
                                100,animalColors{animal},'filled',...
                                'MarkerFaceAlpha',0.8,'XJitter','density','XJitterWidth',0.1); hold on
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

%% Baseline -> reward: plot best fit line for lick during/after paAIP2

%% Baseline -> reward: plot slope for lick during/after paAIP2

%% Baseline -> reward: plot hit rate (bar) during/after paAIP2

%% Baseline -> reward: plot total lick per trial during/after paAIP2 (should be same)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Not used %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
