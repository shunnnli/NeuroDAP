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
% manualWaterIdx = [];
% for i = 1:length(cueIdx)
%     cloestWater = min(abs(cueIdx(i) - waterIdx))/options.behaviorFs;
%     if cloestWater <= 0.05; manualWaterIdx = [manualWaterIdx; i]; end
% end

manualWaterIdx_cue = [];
manualWaterIdx_water = [];

for i = 1:length(cueIdx)
    [minDiff,minIdx] =  min(abs(cueIdx(i) - waterIdx));
    closestWater = minDiff/options.behaviorFs;
    if closestWater <= 0.05 
        manualWaterIdx_cue = [manualWaterIdx_cue; i]; 
        manualWaterIdx_water = [manualWaterIdx_water; minIdx]; 
    end
end

allTrials = cueIdx; allTrials(manualWaterIdx_cue) = [];
waterIdx(manualWaterIdx_water) = [];

%% Debugging figure 1: plot raw rising edges for cue and water
% % tmp = rightSolenoidON(rightSolenoidON>=cur_cue & rightSolenoidON<next_cue);
% % plot the licks to see if they are good
% figure;
% data = zeros(size(allTrials));
% data(allTrials)=1;
% plot(1:length(data),data,'b');
% legend('cue');
% hold on;
% data2 = zeros(size(waterIdx));
% data2(waterIdx)=1;
% plot(1:length(data2),data2,'r');
% legend('water');
% hold on;
% data3 = zeros(size(tmp));
% data3(tmp)=1;
% plot(1:length(data3),data3,'g');
% legend('rewardTimeRaw');

%%
events{1} = allTrials; events{2} = waterIdx;
events{3} = lickIdx;

trials = getTrialTable_shijiaCatch(events);

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
rewardedLickTime = trials{strcmpi(trials.TrialType,'Rewarded'),"RewardLickTime"};

rewarded_choiceLicks = trials{strcmpi(trials.TrialType,'Rewarded'),"ChoiceLicks"};
rewarded_outcomeLicks = trials{strcmpi(trials.TrialType,'Rewarded'),"OutcomeLicks"};
rewarded_spontaneousLicks = trials{strcmpi(trials.TrialType,'Rewarded'),"SpontaneousLicks"};
rewarded_choiceLicks = rewarded_choiceLicks(~cellfun('isempty',rewarded_choiceLicks));
rewarded_outcomeLicks = rewarded_outcomeLicks(~cellfun('isempty',rewarded_outcomeLicks));
rewarded_spontaneousLicks = rewarded_spontaneousLicks(~cellfun('isempty',rewarded_spontaneousLicks));
rewarded_firstLick = cellfun(@(x) x(1),rewarded_choiceLicks);
rewarded_fourthLick = cellfun(@(x) x(1),rewarded_outcomeLicks);
rewarded_firstBoutLastLick = cellfun(@(x) x(end),rewarded_outcomeLicks);
rewarded_trialLastLick = cellfun(@(x) x(end),rewarded_spontaneousLicks);

omissionTrials = trials{strcmpi(trials.TrialType,'Omission'),["TrialNumber","CueTime"]};
omissionIdx = omissionTrials(:,2);
omission_choiceLicks = trials{strcmpi(trials.TrialType,'Omission'),"ChoiceLicks"};
omission_outcomeLicks = trials{strcmpi(trials.TrialType,'Omission'),"OutcomeLicks"};
omission_spontaneousLicks = trials{strcmpi(trials.TrialType,'Omission'),"SpontaneousLicks"};
omission_choiceLicks = omission_choiceLicks(~cellfun('isempty',omission_choiceLicks));
omission_outcomeLicks = omission_outcomeLicks(~cellfun('isempty',omission_outcomeLicks));
omission_spontaneousLicks = omission_spontaneousLicks(~cellfun('isempty',omission_spontaneousLicks));
omission_firstLick = cellfun(@(x) x(1),omission_choiceLicks);
omission_fourthLick = cellfun(@(x) x(1),omission_outcomeLicks);
omission_firstBoutLastLick = cellfun(@(x) x(end),omission_outcomeLicks);
omission_trialLastLick = cellfun(@(x) x(end),omission_spontaneousLicks);

missTrials = trials{strcmpi(trials.TrialType,'Miss'),["TrialNumber","CueTime"]};
missIdx = missTrials(:,2);
miss_choiceLicks = trials{strcmpi(trials.TrialType,'Miss'),"ChoiceLicks"};
miss_outcomeLicks = trials{strcmpi(trials.TrialType,'Miss'),"OutcomeLicks"};
miss_spontaneousLicks = trials{strcmpi(trials.TrialType,'Miss'),"SpontaneousLicks"};
miss_choiceLicks = miss_choiceLicks(~cellfun('isempty',miss_choiceLicks));
miss_outcomeLicks = miss_outcomeLicks(~cellfun('isempty',miss_outcomeLicks));
miss_spontaneousLicks = miss_spontaneousLicks(~cellfun('isempty',miss_spontaneousLicks));
miss_firstLick = cellfun(@(x) x(1),miss_choiceLicks);
miss_fourthLick = cellfun(@(x) x(1),miss_outcomeLicks);
miss_firstBoutLastLick = cellfun(@(x) x(end),miss_outcomeLicks);
miss_trialLastLick = cellfun(@(x) x(end),miss_spontaneousLicks);


% Get baseline (random sample from the session)
randomMinSample = 15*params.sync.behaviorFs;
randomMaxSample = length(params.sync.timePhotometry) - (15*params.sync.behaviorFs);
baselineIdx = randi([randomMinSample,randomMaxSample],100,1);

%% Catch Trial Experiment
%% For catch trials - normal lick 4

rewardedTrials_normal4 = trials{strcmpi(trials.RewardType,'lick4'),["TrialNumber","CueTime"]};
rewarded_normal4_rewardTime = trials{strcmpi(trials.RewardType,'lick4'),"RewardLickTime"}; %% THIS 'RewardLickTime' (LICK REWARD+1 if exist) IS USED FOR PLOTTING ACROSS EXPERIMENTS!!

rewarded_normal4Idx = rewardedTrials_normal4(:,2);
rewarded_normal4_choiceLicks = trials{strcmpi(trials.RewardType,'lick4'),"ChoiceLicks"};
rewarded_normal4_outcomeLicks = trials{strcmpi(trials.RewardType,'lick4'),"OutcomeLicks"};
rewarded_normal4_spontaneousLicks = trials{strcmpi(trials.RewardType,'lick4'),"SpontaneousLicks"};
rewarded_normal4_choiceLicks = rewarded_normal4_choiceLicks(~cellfun('isempty',rewarded_normal4_choiceLicks));
rewarded_normal4_outcomeLicks = rewarded_normal4_outcomeLicks(~cellfun('isempty',rewarded_normal4_outcomeLicks));
rewarded_normal4_spontaneousLicks = rewarded_normal4_spontaneousLicks(~cellfun('isempty',rewarded_normal4_spontaneousLicks));
rewarded_normal4_firstLick = cellfun(@(x) x(1),rewarded_normal4_choiceLicks);
% rewarded_normal4_secondLick = cellfun(@(x) x(2),rewarded_normal4_choiceLicks);
rewarded_normal4_fourthLick = cellfun(@(x) x(1),rewarded_normal4_outcomeLicks);

% rewarded_normal4_fifthLick = cellfun(@(x) x(2),rewarded_normal4_outcomeLicks);

rewarded_normal4_firstBoutLastLick = cellfun(@(x) x(end),rewarded_normal4_outcomeLicks);
rewarded_normal4_trialLastLick = cellfun(@(x) x(end),rewarded_normal4_spontaneousLicks);

%% For catch trials - catch 2

rewardedTrials_catch2 = trials{strcmpi(trials.RewardType,'lick2'),["TrialNumber","CueTime"]};
rewarded_catch2_rewardTime = trials{strcmpi(trials.RewardType,'lick2'),"RewardLickTime"};

rewarded_catch2Idx = rewardedTrials_catch2(:,2);
rewarded_catch2_choiceLicks = trials{strcmpi(trials.RewardType,'lick2'),"ChoiceLicks"};
rewarded_catch2_outcomeLicks = trials{strcmpi(trials.RewardType,'lick2'),"OutcomeLicks"};
rewarded_catch2_spontaneousLicks = trials{strcmpi(trials.RewardType,'lick2'),"SpontaneousLicks"};
rewarded_catch2_choiceLicks = rewarded_catch2_choiceLicks(~cellfun('isempty',rewarded_catch2_choiceLicks));
rewarded_catch2_outcomeLicks = rewarded_catch2_outcomeLicks(~cellfun('isempty',rewarded_catch2_outcomeLicks));
rewarded_catch2_spontaneousLicks = rewarded_catch2_spontaneousLicks(~cellfun('isempty',rewarded_catch2_spontaneousLicks));
rewarded_catch2_firstLick = cellfun(@(x) x(1),rewarded_catch2_choiceLicks);
rewarded_catch2_secondLick = cellfun(@(x) x(2),rewarded_catch2_choiceLicks);

% rewarded_catch2_thirdLick = cellfun(@(x) x(3),rewarded_catch2_choiceLicks);

rewarded_catch2_fourthLick = cellfun(@(x) x(1),rewarded_catch2_outcomeLicks);
rewarded_catch2_firstBoutLastLick = cellfun(@(x) x(end),rewarded_catch2_outcomeLicks);
rewarded_catch2_trialLastLick = cellfun(@(x) x(end),rewarded_catch2_spontaneousLicks);

%% For catch trials - catch 6

rewardedTrials_catch6 = trials{strcmpi(trials.RewardType,'lick6'),["TrialNumber","CueTime"]};
rewarded_catch6_rewardTime = trials{strcmpi(trials.RewardType,'lick6'),"RewardLickTime"};

rewarded_catch6Idx = rewardedTrials_catch6(:,2);
rewarded_catch6_choiceLicks = trials{strcmpi(trials.RewardType,'lick6'),"ChoiceLicks"};
rewarded_catch6_outcomeLicks = trials{strcmpi(trials.RewardType,'lick6'),"OutcomeLicks"};
rewarded_catch6_spontaneousLicks = trials{strcmpi(trials.RewardType,'lick6'),"SpontaneousLicks"};
rewarded_catch6_choiceLicks = rewarded_catch6_choiceLicks(~cellfun('isempty',rewarded_catch6_choiceLicks));
rewarded_catch6_outcomeLicks = rewarded_catch6_outcomeLicks(~cellfun('isempty',rewarded_catch6_outcomeLicks));
rewarded_catch6_spontaneousLicks = rewarded_catch6_spontaneousLicks(~cellfun('isempty',rewarded_catch6_spontaneousLicks));
rewarded_catch6_firstLick = cellfun(@(x) x(1),rewarded_catch6_choiceLicks);
rewarded_catch6_secondLick = cellfun(@(x) x(2),rewarded_catch6_choiceLicks);
rewarded_catch6_fourthLick = cellfun(@(x) x(1),rewarded_catch6_outcomeLicks);
rewarded_catch6_sixthLick = cellfun(@(x) x(3),rewarded_catch6_outcomeLicks);

% rewarded_catch6_seventhLick = cellfun(@(x) x(4),rewarded_catch6_outcomeLicks);

rewarded_catch6_firstBoutLastLick = cellfun(@(x) x(end),rewarded_catch6_outcomeLicks);
rewarded_catch6_trialLastLick = cellfun(@(x) x(end),rewarded_catch6_spontaneousLicks);



%% Debugging figure 2: plot water reward timing
% % RewardTimeDebugPlotIdx = rewardedTrials_catch6(:,2);
% RewardTimeDebugPlotIdx = trials{strcmpi(trials.RewardType,'lick6'),"RewardTimeDebugPlot"};
% analysisEvents_0 = {RewardTimeDebugPlotIdx};
% analysisLabels_0 = {'DebugWaterReward'};
% taskLegend_0 = getLegend(analysisEvents_0,analysisLabels_0);
% 
% timeRange = [-0.5,3];
% 
% % Find the number of photometry channels
% % photometryIdx = find(cellfun(@(x) contains(x,["NI","LJ"],"IgnoreCase",true), {timeSeries.system}));
% photometryIdx = find(cellfun(@(x) contains(x,["Green","Iso"],"IgnoreCase",true), {timeSeries.name}));
% 
% % photometryName = cellfun(@(x) unique(x,'rows'), {timeSeries(photometryIdx).name},'UniformOutput',false);
% nSignals = length(photometryIdx);
% disp(['Finished: found ', num2str(nSignals),' photometry signals']);
% 
% if options.plotPhotometry
% 
%     for photometry = 1:nSignals
%         % Load signal of interest
%         path = photometryIdx(photometry);
%         signal = timeSeries(path).data;
%         finalFs = timeSeries(path).finalFs;
%         system = timeSeries(path).system;
% 
% %         initializeFig(.5,.5); tiledlayout('flow');
%         initializeFig(.5,1); 
%         % tiledlayout(5,2);
% 
%         % Figure set 0: for cue and water events
%         nexttile
%         [~,~] = plotTraces(RewardTimeDebugPlotIdx,timeRange,signal,bluePurpleRed(1,:),params,...
%                                 signalFs=finalFs,...
%                                 signalSystem=system,eventSystem=params.session.baselineSystem);
%         plotEvent('',0);
%         xlabel('Time to rewardTimeDebug (s)'); ylabel('z-score');
%         legend(taskLegend_0,'Location','northeast');
%         title('rewardTimeDebug, neural');
% 
%         nexttile
%         plotLicks(RewardTimeDebugPlotIdx,timeRange,options.lick_binSize,bluePurpleRed(1,:),[],labjack.lick,params);
%         plotEvent('',0);
%         xlabel('Time to rewardTimeDebug (s)'); ylabel('Licks/s');
%         legend(taskLegend_0,'Location','northeast');
%         title('rewardTimeDebug, lick');
% 
%         % Figure set 1: for cue and water events
%         nexttile
%         [~,~] = plotTraces(waterIdx,timeRange,signal,bluePurpleRed(1,:),params,...
%                                 signalFs=finalFs,...
%                                 signalSystem=system,eventSystem=params.session.baselineSystem);
%         [~,~] = plotTraces(rewardedIdx,timeRange,signal,bluePurpleRed(100,:),params,...
%                         signalFs=finalFs,...
%                         signalSystem=system,eventSystem=params.session.baselineSystem);
%         [~,~] = plotTraces(omissionIdx,timeRange,signal,bluePurpleRed(350,:),params,...
%                         signalFs=finalFs,...
%                         signalSystem=system,eventSystem=params.session.baselineSystem);
% 
%         plotEvent('',0);
%         xlabel('Time to cue/water (s)'); ylabel('z-score');
%         legend(taskLegend_1,'Location','northeast');
%         title('Cue/water, neural');
% 
%         nexttile
%         plotLicks(waterIdx,timeRange,options.lick_binSize,bluePurpleRed(1,:),[],labjack.lick,params);
%         plotLicks(rewardedIdx,timeRange,options.lick_binSize,bluePurpleRed(100,:),[],labjack.lick,params);
%         plotLicks(omissionIdx,timeRange,options.lick_binSize,bluePurpleRed(350,:),[],labjack.lick,params);
%         plotEvent('',0);
%         xlabel('Time to cue/water (s)'); ylabel('Licks/s');
%         legend(taskLegend_1,'Location','northeast');
%         title('Cue/water, lick');
%     end
% end


%% Analyze traces

% Figure set 1: for cue and water events
analysisEvents_1 = {rewardedLickTime,rewardedIdx,omissionIdx,missIdx,baselineIdx};
eventTrialNum_1 = {findTrials(rewardedLickTime,trials),...
                 rewardedTrials(:,1),omissionTrials(:,1),missTrials(:,1),...
                 findTrials(baselineIdx,trials)};
analysisLabels_1 = {'FirstWaterCollectionLick','Rewarded','Omission','Miss','Baseline'};
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

% % Figure set 5: for catch trials
% analysisEvents_5 = {rewarded_normal4_fourthLick,rewarded_catch2_secondLick};
% eventTrialNum_5 = {findTrials(rewarded_normal4_fourthLick,trials),...
%                  findTrials(rewarded_catch2_secondLick,trials)};
% analysisLabels_5 = {'Normal 4','Catch 2'};
% taskLegend_5 = getLegend(analysisEvents_5,analysisLabels_5);

% Figure set 5: for catch trials
analysisEvents_5 = {rewarded_normal4_rewardTime,rewarded_catch2_rewardTime,rewarded_catch6_rewardTime};
eventTrialNum_5 = {findTrials(rewarded_normal4_rewardTime,trials),...
                 findTrials(rewarded_catch2_rewardTime,trials),...
                 findTrials(rewarded_catch6_rewardTime,trials)};
analysisLabels_5 = {'Normal 4','Catch 2','Catch 6'};
taskLegend_5 = getLegend(analysisEvents_5,analysisLabels_5);

% % Figure set 6: for catch trials, first water collection lick
% analysisEvents_5 = {rewarded_normal4_rewardTime,rewarded_catch2_rewardTime,rewarded_catch6_rewardTime};
% eventTrialNum_5 = {findTrials(rewarded_normal4_rewardTime,trials),...
%                  findTrials(rewarded_catch2_rewardTime,trials),...
%                  findTrials(rewarded_catch6_rewardTime,trials)};
% analysisLabels_5 = {'Normal 4','Catch 2','Catch 6'};
% taskLegend_5 = getLegend(analysisEvents_5,analysisLabels_5);
% 

% to-do: concatenate analysis events
analysisEvents = [analysisEvents_1, analysisEvents_2, analysisEvents_3, analysisEvents_4, analysisEvents_5];
analysisLabels = [analysisLabels_1, analysisLabels_2, analysisLabels_3, analysisLabels_4, analysisLabels_5]; 
eventTrialNum = [eventTrialNum_1, eventTrialNum_2, eventTrialNum_3, eventTrialNum_4, eventTrialNum_5];

% clean licks with ITI 90ms, somehow it's not really working for plotting 
lickFallingEdges = find(labjack.lick);
previousLickOnset = lickFallingEdges(1);
lickITICutoffArduino = 90;
lickFallingEdges_cleaned = [previousLickOnset];
for p = 1:length(lickFallingEdges)
    if lickFallingEdges(p)-previousLickOnset>(lickITICutoffArduino/1000*labjack.samplerate)
        previousLickOnset = lickFallingEdges(p);
        lickFallingEdges_cleaned = [lickFallingEdges_cleaned previousLickOnset];
    end
end
lick_labjack_zeroed = zeros(size(labjack.lick));
lick_labjack_zeroed(lickFallingEdges_cleaned)=1;
labjack.lick = lick_labjack_zeroed;
% figure;plot(lick_labjack_tmp(25000:30000),'b');hold on;plot(lick_labjack_zeroed(25000:30000),'r');

if options.analyzeTraces
    analysis = analyzeTraces(timeSeries,labjack.lick,analysisEvents,analysisLabels,params,...
                            trialNumber=eventTrialNum,trialTable=trials);
end

%% Test plotting traces

timeRange = [-0.5,3];

% Find the number of photometry channels
% photometryIdx = find(cellfun(@(x) contains(x,["NI","LJ"],"IgnoreCase",true), {timeSeries.system}));
photometryIdx = find(cellfun(@(x) contains(x,["Green","Iso"],"IgnoreCase",true), {timeSeries.name}));

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
        initializeFig(.5,1); tiledlayout(5,2);

        % Figure set 1: for cue and water events
        nexttile
        [~,~] = plotTraces(rewardedLickTime,timeRange,signal,bluePurpleRed(1,:),params,...
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
        xlabel('Time to cue/waterlick (s)'); ylabel('z-score');
        legend(taskLegend_1,'Location','northeast');
        title('Cue/water, neural');
        
        nexttile
        plotLicks(rewardedLickTime,timeRange,options.lick_binSize,bluePurpleRed(1,:),[],labjack.lick,params);
        plotLicks(rewardedIdx,timeRange,options.lick_binSize,bluePurpleRed(100,:),[],labjack.lick,params);
        plotLicks(omissionIdx,timeRange,options.lick_binSize,bluePurpleRed(350,:),[],labjack.lick,params);
        plotLicks(missIdx,timeRange,options.lick_binSize,bluePurpleRed(500,:),[],labjack.lick,params);
        plotLicks(baselineIdx,timeRange,options.lick_binSize,[.5,.5,.5],[],labjack.lick,params);
        plotEvent('',0);
        xlabel('Time to cue/waterlick (s)'); ylabel('Licks/s');
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

%         % Figure set 5: for catch trial reward lick 
%         nexttile
%         [~,~] = plotTraces(rewarded_normal4_fourthLick,timeRange,signal,bluePurpleRed(1,:),params,...
%                                 signalFs=finalFs,...
%                                 signalSystem=system,eventSystem=params.session.baselineSystem);
%         [~,~] = plotTraces(rewarded_catch2_secondLick,timeRange,signal,bluePurpleRed(500,:),params,...
%                         signalFs=finalFs,...
%                         signalSystem=system,eventSystem=params.session.baselineSystem);
%         
%         plotEvent('',0);
%         xlabel('Time to reward lick (s)'); ylabel('z-score');
%         legend(taskLegend_5,'Location','northeast');
%         title('Catch Experiment: reward lick, neural');
% 
%         nexttile
%         plotLicks(rewarded_normal4_fourthLick,timeRange,options.lick_binSize,bluePurpleRed(1,:),[],labjack.lick,params);
%         plotLicks(rewarded_catch2_secondLick,timeRange,options.lick_binSize,bluePurpleRed(500,:),[],labjack.lick,params);
%        
%         plotEvent('',0);
%         xlabel('Time to reward lick (s)'); ylabel('Licks/s');
%         legend(taskLegend_5,'Location','northeast');
%         title('Catch Experiment: reward lick, lick');

        % Figure set 5: for catch trial REWARD LICK onset 
        nexttile
        [~,~] = plotTraces(rewarded_normal4_rewardTime,timeRange,signal,bluePurpleRed(1,:),params,...
                                signalFs=finalFs,...
                                signalSystem=system,eventSystem=params.session.baselineSystem);
        [~,~] = plotTraces(rewarded_catch2_rewardTime,timeRange,signal,bluePurpleRed(250,:),params,...
                        signalFs=finalFs,...
                        signalSystem=system,eventSystem=params.session.baselineSystem);
        [~,~] = plotTraces(rewarded_catch6_rewardTime,timeRange,signal,bluePurpleRed(500,:),params,...
                signalFs=finalFs,...
                signalSystem=system,eventSystem=params.session.baselineSystem);
        
        plotEvent('',0);
        xlabel('Time to first water collection lick onset (s)'); ylabel('z-score');
        legend(taskLegend_5,'Location','northeast');
        title('Catch Experiment: first water collection lick onset, neural');

        nexttile
        plotLicks(rewarded_normal4_rewardTime,timeRange,options.lick_binSize,bluePurpleRed(1,:),[],labjack.lick,params);
        plotLicks(rewarded_catch2_rewardTime,timeRange,options.lick_binSize,bluePurpleRed(250,:),[],labjack.lick,params);
        plotLicks(rewarded_catch6_rewardTime,timeRange,options.lick_binSize,bluePurpleRed(500,:),[],labjack.lick,params);

        plotEvent('',0);
        xlabel('Time to first water collection lick onset (s)'); ylabel('Licks/s');
        legend(taskLegend_5,'Location','northeast');
        title('Catch Experiment: first water collection lick onset, lick');

        saveas(gcf,strcat(sessionpath,filesep,'Summary_events_',timeSeries(path).name,'.png'));

    end
end
1;

%% Plot lick scatter plot

