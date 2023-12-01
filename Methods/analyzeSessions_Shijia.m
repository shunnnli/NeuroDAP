function analyzeSessions_Shijia(sessionpath,options)

arguments
    sessionpath string
    options.analyzeTraces logical = true
    options.redo logical = true % Recalculate trial table and all preprocessing
    
    options.plotPhotometry logical = true % Plot photometry summary plot
    options.plotLicks logical = true % Plot lick raster summary plot

    options.behaviorFs double = 2000

    options.lick_binSize double = 0.1
end

%% Notes
% Modified from Shun_analyzeBehavior_optoPair
% Shun Li, 11/20/2023
% 02/14/2023: tidied up code, renamed to analyzeBehavior_optoPair
% 2023/07/28: packaged trial table into a function
% 2023/09/02: added camera plotting
% 2023/09/05: changed baselineIdx to selecting baseline licks 
% 2023/10/23: changed how to plot photometry signal, assume everything
% recorded in labjack

%% Load data

[~,~,~,~,~,~,bluePurpleRed] = loadColors;
             
% 1. Select session via uigetdir
dirsplit = strsplit(sessionpath,filesep); 
sessionName = dirsplit{end}; 
if ispc; projectPath = strcat('\\',fullfile(dirsplit{2:end-1}));
elseif isunix; projectPath = strcat('/',fullfile(dirsplit{2:end-1}));
end
% Get animal name and session date
dirsplit = strsplit(sessionName,'_'); % 20231122 shijia changed '-' to '_'
date = dirsplit{1}; animal = dirsplit{2};
clear dirsplit

disp(strcat('**********',sessionName,'**********'));
load(strcat(sessionpath,filesep,'timeseries_',sessionName,'.mat'));
load(strcat(sessionpath,filesep,'data_',sessionName,'.mat'));
load(strcat(sessionpath,filesep,'behavior_',sessionName,'.mat'));
load(strcat(sessionpath,filesep,'sync_',sessionName,'.mat'));

if ~isfield(params.session,'name'); params.session.name = sessionName; end
if ~isfield(params.session,'date'); params.session.date = date; end
if ~isfield(params.session,'animal'); params.session.animal = animal; end
if ~isfield(params.session,'projectPath'); params.session.projectPath = projectPath; end

% Create analysis.mat
if ~isempty(dir(fullfile(sessionpath,"analysis_*.mat")))
    load(strcat(sessionpath,filesep,'analysis_',sessionName,'.mat'));
else
    save(strcat(sessionpath,filesep,'analysis_',sessionName),'sessionName','-v7.3');
    disp('Finished: analysis_.mat not found, created a new one');
end
disp(['Finished: Session ',sessionName,' loaded']);

%% Generate trial table

cueIdx = find(labjack.cue); 
waterIdx = find(labjack.solenoid);
lickIdx = find(labjack.lick);
% laserIdx = find(labjack.sync);

% Get allTrials (delete manual water)
manualWaterIdx = [];
for i = 1:length(cueIdx)
    cloestWater = min(abs(cueIdx(i) - waterIdx))/options.behaviorFs;
    if cloestWater <= 0.05; manualWaterIdx = [manualWaterIdx; i]; end
end
allTrials = cueIdx; allTrials(manualWaterIdx) = [];

events{1} = allTrials; events{2} = waterIdx;
events{3} = lickIdx;

trials = getTrialTable_shijia(events);

% For converting to datajoint
disp('Ongoing: create tables for datajoint pipeline');
% trialTable
trialTable = trials(:,1:end-5);
trialTable.block = ones(size(trials,1),1);
trialTable.session_position = (1:size(trials,1))';
trialTable = replaceNaN(trialTable,-1);
parquetwrite(strcat(sessionpath,filesep,'trialTable.parquet'),trialTable);
% eventTable
% eventTable = getEventTable(events,params);
% parquetwrite(strcat(sessionpath,filesep,'eventTable.parquet'),eventTable);
% blockTable
firstTrial = 1; lastTrial = size(trials,1);
blockTable = table(firstTrial,lastTrial);
parquetwrite(strcat(sessionpath,filesep,'blockTable.parquet'),blockTable);


% Save to behavior_.mat
% save(strcat(sessionpath,filesep,'behavior_',params.session.name),'trials',...
%     'eventTable','trialTable','blockTable','-append');
disp('Finished: trial table saved');

%% Select events of interest

rewardedTrials = trials{strcmpi(trials.TrialType,'Rewarded'),["TrialNumber","CueTime"]};
rewardedIdx = rewardedTrials(:,2);
rewarded_choiceLicks = trials{strcmpi(trials.TrialType,'Rewarded'),"ChoiceLicks"};
rewarded_outcomeLicks = trials{strcmpi(trials.TrialType,'Rewarded'),"OutcomeLicks"};
rewarded_baselineLicks = trials{strcmpi(trials.TrialType,'Rewarded'),"BaselineLicks"};
rewarded_choiceLicks = rewarded_choiceLicks(~cellfun('isempty',rewarded_choiceLicks));
rewarded_outcomeLicks = rewarded_outcomeLicks(~cellfun('isempty',rewarded_outcomeLicks));
rewarded_baselineLicks = rewarded_baselineLicks(~cellfun('isempty',rewarded_baselineLicks));
rewarded_firstLick = cellfun(@(x) x(1),rewarded_choiceLicks);
rewarded_fourthLick = cellfun(@(x) x(1),rewarded_outcomeLicks);
rewarded_firstBoutLastLick = cellfun(@(x) x(end),rewarded_outcomeLicks);
rewarded_trialLastLick = cellfun(@(x) x(end),rewarded_baselineLicks);

omissionTrials = trials{strcmpi(trials.TrialType,'Omission'),["TrialNumber","CueTime"]};
omissionIdx = omissionTrials(:,2);
omission_choiceLicks = trials{strcmpi(trials.TrialType,'Omission'),"ChoiceLicks"};
omission_outcomeLicks = trials{strcmpi(trials.TrialType,'Omission'),"OutcomeLicks"};
omission_baselineLicks = trials{strcmpi(trials.TrialType,'Omission'),"BaselineLicks"};
omission_choiceLicks = omission_choiceLicks(~cellfun('isempty',omission_choiceLicks));
omission_outcomeLicks = omission_outcomeLicks(~cellfun('isempty',omission_outcomeLicks));
omission_baselineLicks = omission_baselineLicks(~cellfun('isempty',omission_baselineLicks));
omission_firstLick = cellfun(@(x) x(1),omission_choiceLicks);
omission_fourthLick = cellfun(@(x) x(1),omission_outcomeLicks);
omission_firstBoutLastLick = cellfun(@(x) x(end),omission_outcomeLicks);
omission_trialLastLick = cellfun(@(x) x(end),omission_baselineLicks);

missTrials = trials{strcmpi(trials.TrialType,'Miss'),["TrialNumber","CueTime"]};
missIdx = missTrials(:,2);
miss_choiceLicks = trials{strcmpi(trials.TrialType,'Miss'),"ChoiceLicks"};
miss_outcomeLicks = trials{strcmpi(trials.TrialType,'Miss'),"OutcomeLicks"};
miss_baselineLicks = trials{strcmpi(trials.TrialType,'Miss'),"BaselineLicks"};
miss_choiceLicks = miss_choiceLicks(~cellfun('isempty',miss_choiceLicks));
miss_outcomeLicks = miss_outcomeLicks(~cellfun('isempty',miss_outcomeLicks));
miss_baselineLicks = miss_baselineLicks(~cellfun('isempty',miss_baselineLicks));
miss_firstLick = cellfun(@(x) x(1),miss_choiceLicks);
miss_fourthLick = cellfun(@(x) x(1),miss_outcomeLicks);
miss_firstBoutLastLick = cellfun(@(x) x(end),miss_outcomeLicks);
miss_trialLastLick = cellfun(@(x) x(end),miss_baselineLicks);


% Get baseline (random sample from the session)
randomMinSample = 15*params.sync.behaviorFs;
randomMaxSample = length(params.sync.timePhotometry) - (15*params.sync.behaviorFs);
baselineIdx = randi([randomMinSample,randomMaxSample],100,1);

%% Analyze traces
% Figure set 1: for cue and water events
analysisEvents_1 = {waterIdx,rewardedIdx,omissionIdx,missIdx,baselineIdx};
eventTrialNum_1 = {findTrials(waterIdx,trials),...
                 rewardedTrials(:,1),omissionTrials(:,1),missTrials(:,1),...
                 findTrials(baselineIdx,trials)};
analysisLabels_1 = {'Water','Rewarded','Omission','Miss','Baseline'};
taskLegend_1 = getLegend(analysisEvents_1,analysisLabels_1);

% Figure set 2: for first lick
analysisEvents_2 = {rewarded_firstLick,omission_firstLick};
eventTrialNum_2 = {findTrials(rewarded_firstLick,trials),...
                 findTrials(omission_firstLick,trials)};
analysisLabels_2 = {'Rewarded first lick','Omission first lick'};
taskLegend_2 = getLegend(analysisEvents_2,analysisLabels_2);

% Figure set 3: for fourth lick
analysisEvents_3 = {rewarded_fourthLick,omission_fourthLick};
eventTrialNum_3 = {findTrials(rewarded_fourthLick,trials),...
                 findTrials(omission_fourthLick,trials)};
analysisLabels_3 = {'Rewarded fourth lick','Omission fourth lick'};
taskLegend_3 = getLegend(analysisEvents_3,analysisLabels_3);


% Figure set 4: for first bout last lick
analysisEvents_4 = {rewarded_firstBoutLastLick,omission_firstBoutLastLick};
eventTrialNum_4 = {findTrials(rewarded_firstBoutLastLick,trials),...
                 findTrials(omission_firstBoutLastLick,trials)};
analysisLabels_4 = {'Rewarded first bout last lick','Omission first bout last lick'};
taskLegend_4 = getLegend(analysisEvents_4,analysisLabels_4);


if options.analyzeTraces
    analysis = analyzeTraces(timeSeries,labjack.lick,analysisEvents,analysisLabels,params,...
                            trialNumber=eventTrialNum,trialTable=trials);
end

%% Test plotting traces

timeRange = [-0.5,3];

% Find the number of photometry channels
photometryIdx = find(cellfun(@(x) contains(x,["NI","LJ"],"IgnoreCase",true), {timeSeries.system}));
% photometryName = cellfun(@(x) unique(x,'rows'), {timeSeries(photometryIdx).name},'UniformOutput',false);
nSignals = length(photometryIdx);
disp(['Finished: found ', num2str(nSignals),' photometry signals']);

if options.plotPhotometry

    for photometry = 1:nSignals
        % Load signal of interest
        path = photometryIdx(photometry);
        signal = timeSeries(path).data;
        finalFs = timeSeries(path).finalFs;
        system = timeSeries(path).system;

%         initializeFig(.5,.5); tiledlayout('flow');
        initializeFig(.5,1); tiledlayout(4,2);

        % Figure set 1: for cue and water events
        nexttile
        [~,~] = plotTraces(waterIdx,timeRange,signal,bluePurpleRed(1,:),params,...
                                signalFs=finalFs,...
                                signalSystem=system,eventSystem=params.session.baselineSystem);
        [~,~] = plotTraces(rewardedIdx,timeRange,signal,bluePurpleRed(100,:),params,...
                        signalFs=finalFs,...
                        signalSystem=system,eventSystem=params.session.baselineSystem);
        [~,~] = plotTraces(omissionIdx,timeRange,signal,bluePurpleRed(350,:),params,...
                        signalFs=finalFs,...
                        signalSystem=system,eventSystem=params.session.baselineSystem);
        [~,~] = plotTraces(missIdx,timeRange,signal,bluePurpleRed(500,:),params,...
                        signalFs=finalFs,...
                        signalSystem=system,eventSystem=params.session.baselineSystem);
        [~,~] = plotTraces(baselineIdx,timeRange,signal,[.5,.5,.5],params,...
                signalFs=finalFs,...
                signalSystem=system,eventSystem=params.session.baselineSystem);
        plotEvent('',0);
        xlabel('Time to cue/water (s)'); ylabel('z-score');
        legend(taskLegend_1,'Location','northeast');
        title('Cue/water, neural');
        
        nexttile
        plotLicks(waterIdx,timeRange,options.lick_binSize,bluePurpleRed(1,:),[],labjack.lick,params);
        plotLicks(rewardedIdx,timeRange,options.lick_binSize,bluePurpleRed(100,:),[],labjack.lick,params);
        plotLicks(omissionIdx,timeRange,options.lick_binSize,bluePurpleRed(350,:),[],labjack.lick,params);
        plotLicks(missIdx,timeRange,options.lick_binSize,bluePurpleRed(500,:),[],labjack.lick,params);
        plotLicks(baselineIdx,timeRange,options.lick_binSize,[.5,.5,.5],[],labjack.lick,params);
        plotEvent('',0);
        xlabel('Time to cue/water (s)'); ylabel('Licks/s');
        legend(taskLegend_1,'Location','northeast');
        title('Cue/water, lick');

        % Figure set 2: for first lick
        nexttile
        [~,~] = plotTraces(rewarded_firstLick,timeRange,signal,bluePurpleRed(1,:),params,...
                                signalFs=finalFs,...
                                signalSystem=system,eventSystem=params.session.baselineSystem);
        [~,~] = plotTraces(omission_firstLick,timeRange,signal,bluePurpleRed(500,:),params,...
                        signalFs=finalFs,...
                        signalSystem=system,eventSystem=params.session.baselineSystem);
        
        plotEvent('',0);
        xlabel('Time to first lick (s)'); ylabel('z-score');
        legend(taskLegend_2,'Location','northeast');
        title('First lick, neural');

        nexttile
        plotLicks(rewarded_firstLick,timeRange,options.lick_binSize,bluePurpleRed(1,:),[],labjack.lick,params);
        plotLicks(omission_firstLick,timeRange,options.lick_binSize,bluePurpleRed(500,:),[],labjack.lick,params);
       
        plotEvent('',0);
        xlabel('Time to first lick (s)'); ylabel('Licks/s');
        legend(taskLegend_2,'Location','northeast');
        title('First lick, lick');


        % Figure set 3: for fourth lick
        nexttile
        [~,~] = plotTraces(rewarded_fourthLick,timeRange,signal,bluePurpleRed(1,:),params,...
                                signalFs=finalFs,...
                                signalSystem=system,eventSystem=params.session.baselineSystem);
        [~,~] = plotTraces(omission_fourthLick,timeRange,signal,bluePurpleRed(500,:),params,...
                        signalFs=finalFs,...
                        signalSystem=system,eventSystem=params.session.baselineSystem);
        
        plotEvent('',0);
        xlabel('Time to 4th lick (s)'); ylabel('z-score');
        legend(taskLegend_3,'Location','northeast');
        title('Fourth lick, neural');

        
        nexttile
        plotLicks(rewarded_fourthLick,timeRange,options.lick_binSize,bluePurpleRed(1,:),[],labjack.lick,params);
        plotLicks(omission_fourthLick,timeRange,options.lick_binSize,bluePurpleRed(500,:),[],labjack.lick,params);
       
        plotEvent('',0);
        xlabel('Time to 4th lick (s)'); ylabel('Licks/s');
        legend(taskLegend_3,'Location','northeast');
        title('Fourth lick, lick');

        % Figure set 4: for first bout last lick
        nexttile
        [~,~] = plotTraces(rewarded_firstBoutLastLick,timeRange,signal,bluePurpleRed(1,:),params,...
                                signalFs=finalFs,...
                                signalSystem=system,eventSystem=params.session.baselineSystem);
        [~,~] = plotTraces(omission_firstBoutLastLick,timeRange,signal,bluePurpleRed(500,:),params,...
                        signalFs=finalFs,...
                        signalSystem=system,eventSystem=params.session.baselineSystem);
        
        plotEvent('',0);
        xlabel('Time to last lick (s)'); ylabel('z-score');
        legend(taskLegend_4,'Location','northeast');
        title('First bout last lick, neural');

        nexttile
        plotLicks(rewarded_firstBoutLastLick,timeRange,options.lick_binSize,bluePurpleRed(1,:),[],labjack.lick,params);
        plotLicks(omission_firstBoutLastLick,timeRange,options.lick_binSize,bluePurpleRed(500,:),[],labjack.lick,params);
       
        plotEvent('',0);
        xlabel('Time to last lick (s)'); ylabel('Licks/s');
        legend(taskLegend_4,'Location','northeast');
        title('First bout last lick, lick');

        saveas(gcf,strcat(sessionpath,filesep,'Summary_events_',timeSeries(path).name,'.png'));
    end
end

%% Plot lick scatter plot

