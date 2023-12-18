function analyzeSessions_optoPair(sessionpath,options)

arguments
    sessionpath string
    
    options.task
    options.stimPattern cell

    options.pavlovian logical = false
    options.reactionTime double = 1
    options.minLicks double = 2 % min licks to get reward

    options.analyzeTraces logical = true

    options.redo logical = true % Recalculate trial table and all preprocessing
    options.round logical = false % Round reward/airpuff/tone to get duration data
    options.performing logical = false % Only plot traces where the animal performs
    
    options.plotPhotometry logical = true % Plot photometry summary plot
    options.plotBehavior logical = true % Plot lick raster summary plot

    options.lick_binSize double = 0.1
end

%% Notes
% Shun_analyzeBehavior_optoPair
% Shun Li, 11/8/2022
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
dirsplit = strsplit(sessionName,'-');
date = dirsplit{1}; animal = dirsplit{2}; sessionTask = dirsplit{3};
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

% Load behaivor params
if isfield(options,'stimPattern') && (options.redo || ~isfield(params,'stim'))
    params.stim.pulseFreq = str2double(options.stimPattern{1}); 
    params.stim.pulseDuration = str2double(options.stimPattern{2}); 
    params.stim.stimDuration = str2double(options.stimPattern{3});
    params.stim.nPulsesPerStim = (params.stim.stimDuration/1000) * params.stim.pulseFreq;
    params.stim.pulseInterval = (1/params.stim.pulseFreq) - params.stim.pulseDuration;
elseif ~isfield(options,'stimPattern') && ~isfield(params,'stim') 
    params.stim.pulseFreq = 20; 
    params.stim.pulseDuration = 5; 
    params.stim.stimDuration = 500;
    params.stim.nPulsesPerStim = (params.stim.stimDuration/1000) * params.stim.pulseFreq;
    params.stim.pulseInterval = (1/params.stim.pulseFreq) - params.stim.pulseDuration;
end

% Load task
if ~isfield(params.session,'task') && isfield(options,'task')
    params.session.task = options.task; 
elseif isfield(params.session,'task') && isfield(options,'task')
    if ~strcmpi(options.task,params.session.task)
        % options.task = params.session.task;
        params.session.task = options.task;
        disp('Finished: options.task not provided or differed, use the original one');
    end
elseif ~isfield(params.session,'task') && ~isfield(options,'task')
    if contains(sessionTask,"P",IgnoreCase=false)
        options.task = 'punish pairing';
        params.session.task = options.task;
    elseif contains(sessionTask,"Random",IgnoreCase=false)
        options.task = 'random';
        params.session.task = options.task;
    elseif contains(sessionTask,"R",IgnoreCase=false)
        options.task = 'reward pairing';
        params.session.task = options.task;
    else % contains(sessionTask,"Random",IgnoreCase=false)
        options.task = 'random';
        params.session.task = options.task;
    end
    disp(['Finished: parsing session task as ',options.task]);
elseif isfield(params.session,'task') && ~isfield(options,'task')
    options.task = params.session.task;
end

% Define baselineSystem
if ~isfield(params.session,'baselineSystem')
    params.session.baselineSystem = 'NI';
end

% Load analyze options
if options.redo || ~isfield(params,'analyze')
    params.analyze = options;
    options.pavlovian = params.analyze.pavlovian;
    options.reactionTime = params.analyze.reactionTime;
    options.minLicks = params.analyze.minLicks;
    save(strcat(sessionpath,filesep,'sync_',params.session.name),'params','-append');
elseif isfield(params,'analyze')
    options.pavlovian = params.analyze.pavlovian;
    options.reactionTime = params.analyze.reactionTime;
    options.minLicks = params.analyze.minLicks;
end

disp(['Behavior params: task = ',params.session.task]);
disp(['Behavior params: pavlovian = ',num2str(params.analyze.pavlovian)]);
disp(['Behavior params: reactionTime = ',num2str(params.analyze.reactionTime)]);
disp(['Behavior params: minLicks = ',num2str(params.analyze.minLicks)]);

% Load camera signal
eyeAreaIdx = find(cellfun(@(x) strcmpi(x,'eyeArea'), {timeSeries.name}));
if ~isempty(eyeAreaIdx); eyeArea_detrend = timeSeries(eyeAreaIdx).data; end
pupilAreaIdx = find(cellfun(@(x) strcmpi(x,'pupilArea'), {timeSeries.name}));
if ~isempty(pupilAreaIdx); pupilArea_detrend = timeSeries(pupilAreaIdx).data; end

save(strcat(sessionpath,filesep,'sync_',params.session.name),'params','-append');

%% Preprocess outcome and opto data

disp('Ongoing: preprocess outcome and opto data');

% Reward/punishment params
rewardUnit = 0.012; % 8ms opening to dispense 1ul water
rewardList = [0 3 10]; % in ul
% punishList = [0 0.1];
% toneList = [0 0.5 1]; % in sec
if options.round
    if ~exist('rightSolenoid_rounded','var') || options.redo
        % Round reward and tone
        rightSolenoid = rightSolenoid ./ rewardUnit;
        rightSolenoid_rounded = roundToTarget(rightSolenoid, rewardList); disp('Finished rounding: rightSolenoid');
        airpuff_rounded = roundToTarget(airpuff,punishList); disp('Finished rounding: airpuff');
        
        % Save rounded data
        save(strcat(sessionpath,filesep,'timeseries_',params.session.name),'rightSolenoid_rounded','airpuff_rounded','-append');
        disp('Finished: rounding cue/outcome data');
    end
else
    rightSolenoid_rounded = rightSolenoid;
    airpuff_rounded = airpuff;
end


% Find start of opto cue
if ~exist('optoCue','var') || options.redo
    if ~isempty(find(redLaser, 1))
        % Find the first pulse of each stim pattern if nPulsePerPattern>1
        if params.stim.nPulsesPerStim > 1
            allPulses = find(redLaser);
            intervalThreshold = 10000;
            temp_interval = [100000,diff(allPulses)];
            optoCue = allPulses(temp_interval > intervalThreshold);
    
            % Save first pulse data
            save(strcat(sessionpath,filesep,'timeseries_',params.session.name),"optoCue",'-append');
        else
            optoCue = find(redLaser);
            save(strcat(sessionpath,filesep,'timeseries_',params.session.name),"optoCue",'-append');
        end
        disp('Finished: saved optoCue data');
    else 
        optoCue = [];
    end
    save(strcat(sessionpath,filesep,'timeseries_',params.session.name),"optoCue",'-append');
end


% Find start of lick bout
rightLickON = find(rightLick);
if ~exist('lickBout','var') || options.redo
    % Get lick bout start time (ILI < 0.5s)
    lickBout = getLickBout(rightLickON);
    save(strcat(sessionpath,filesep,'timeseries_',params.session.name),"lickBout",'-append');
end


% Combine stim&tone to form trial start
waterIdx = find(rightSolenoid_rounded);  
airpuffIdx = find(airpuff_rounded);

if ~exist('trials','var') || options.redo
    if strcmp(options.task,'random')
        [allTrials,~] = getTrials(find(leftTone),optoCue,...
                             waterIdx,airpuffIdx);
    elseif contains(options.task,'punish')
        [allTrials,~] = getTrials(find(leftTone),optoCue,waterIdx);
    else
        [allTrials,~] = getTrials(find(leftTone),optoCue);
    end
    save(strcat(sessionpath,filesep,'timeseries_',params.session.name),"allTrials",'-append');
end

% Find water lick (first lick in response to water)
waterLickIdx = nan(size(waterIdx));
for i = 1:length(waterIdx)
    nextLick = rightLickON(find(rightLickON>=waterIdx(i),1));
    if ~isempty(nextLick); waterLickIdx(i) = nextLick; end
end
waterLickIdx = rmmissing(waterLickIdx);

disp('Finished: preprocess outcome and opto data');

%% Generate trial and event table

if (~exist('trials','var') || options.redo)

    disp('Ongoing: making trial table');
    events{1} = allTrials;      events{2} = airpuffIdx;
    events{3} = waterIdx;       events{4} = rightLickON;
    events{5} = find(leftTone); events{6} = optoCue;

    trials = getTrialTable(options.task,events,rightSolenoid_rounded,airpuff_rounded,...
                pavlovian=options.pavlovian,reactionTime=options.reactionTime);
    % Calculate performance cutoff
    if contains(options.task,'reward'); [trials,cutoff_sample] = getSessionCutoff(trials,"->reward");
    elseif contains(options.task,'punish'); [trials,cutoff_sample] = getSessionCutoff(trials,"->punish");
    else; [trials,cutoff_sample] = getSessionCutoff(trials,"random"); 
    end
    params.analysis.cutoff_sample = cutoff_sample;
    disp(['     Session cutoff calculated: ',num2str(cutoff_sample)]);


    % For converting to datajoint
    disp('Ongoing: create tables for datajoint pipeline');
    % trialTable
    trialTable = trials(:,1:end-6);
    trialTable.block = ones(size(trials,1),1);
    trialTable.session_position = (1:size(trials,1))';
    trialTable = replaceNaN(trialTable,-1);
    parquetwrite(strcat(sessionpath,filesep,'trialTable.parquet'),trialTable);
    % eventTable
    eventTable = getEventTable(events,params);
    parquetwrite(strcat(sessionpath,filesep,'eventTable.parquet'),eventTable);
    % blockTable
    firstTrial = 1; lastTrial = size(trials,1);
    blockTable = table(firstTrial,lastTrial);
    parquetwrite(strcat(sessionpath,filesep,'blockTable.parquet'),blockTable);
    

    % Save to behavior_.mat
    save(strcat(sessionpath,filesep,'behavior_',params.session.name),'trials',...
        'eventTable','trialTable','blockTable','-append');
    disp('Finished: trial table saved');
end


%% Task specific params

if strcmp(options.task,'random')
    
    % Select event idx
    toneIdx = find(leftTone); stimIdx = optoCue;

    % Select baseline idx (align to each baseline lick)
    % baselineLicks = cell2mat(trials{~cellfun(@isempty,trials{1:end-1,'BaselineLicks'}),'BaselineLicks'});
    % if isempty(baselineLicks); baselineIdx = [];
    % else; baselineIdx = baselineLicks(:,1); end
    randomMinSample = 15*params.sync.behaviorFs;
    randomMaxSample = length(params.sync.timeNI) - (15*params.sync.behaviorFs);
    baselineIdx = randi([randomMinSample,randomMaxSample],100,1);
    
    % Create task legend
    stageTime = [-2,0;0,2];
    analysisEvents = {waterIdx,waterLickIdx,airpuffIdx,toneIdx,stimIdx,baselineIdx};
    eventTrialNum = {findTrials(waterIdx,trials),findTrials(waterLickIdx,trials),...
                    findTrials(airpuffIdx,trials),findTrials(toneIdx,trials),...
                    findTrials(stimIdx,trials),findTrials(baselineIdx,trials)};
    analysisLabels = {'Water','Rewarded licks','Airpuff','Tone','Stim','Baseline'};
    taskLegend = getLegend(analysisEvents,analysisLabels);

    stageColors = {[.75 .75 .75],bluePurpleRed(1,:)};
    stageLegend = {'Baseline','US'};
    
    eventTrialNum = eventTrialNum(~cellfun('isempty',analysisEvents));
    analysisLabels = analysisLabels(~cellfun('isempty',analysisEvents));
    analysisEvents = analysisEvents(~cellfun('isempty',analysisEvents));
    
    for i = 1:length(analysisEvents)
        disp(['Total ',analysisLabels{i},': ',num2str(length(analysisEvents{i}))]);
    end

else
    % Select baseline idx (align to each baseline lick)
    % baselineLicks = cell2mat(trials{~cellfun(@isempty,trials{1:end-1,'BaselineLicks'}),'BaselineLicks'});
    % if isempty(baselineLicks); baselineIdx = round((trials{2:end,"CueTime"} - trials{2:end,"ENL"}) - 5*params.sync.behaviorFs);
    % else; baselineIdx = baselineLicks(:,1); end
    randomMinSample = 15*params.sync.behaviorFs;
    randomMaxSample = length(params.sync.timeNI) - (15*params.sync.behaviorFs);
    baselineIdx = randi([randomMinSample,randomMaxSample],100,1);

    if contains(options.task,'reward')
        if options.performing
            stimTrials = trials{trials.isTone == 0 & trials.isStim == 1 ...
                & trials.isReward == 1 & trials.performing == 1,["TrialNumber","CueTime","OutcomeTime","ENL"]};
            stimOmissionTrials = trials{trials.isTone == 0 & trials.isStim == 1 ...
                & trials.isReward == 0 & trials.performing == 1,["TrialNumber","CueTime","OutcomeTime","ENL"]};
            pairTrials = trials{trials.isTone == 1 & trials.isStim == 1 ...
                & trials.isReward == 1 & trials.performing == 1,["TrialNumber","CueTime","OutcomeTime","ENL"]};
            pairOmissionTrials = trials{trials.isTone == 1 & trials.isStim == 1 ...
                & trials.isReward == 0 & trials.performing == 1,["TrialNumber","CueTime","OutcomeTime","ENL"]};
            toneTrials = trials{trials.isTone == 1 & trials.isStim == 0 ...
                & trials.performing == 1,["TrialNumber","CueTime","OutcomeTime","ENL"]};
            stimIdx = stimTrials(:,2);
            stimOmissionIdx = stimOmissionTrials(:,2);
            pairIdx = pairTrials(:,2);
            pairOmissionIdx = pairOmissionTrials(:,2);
            toneIdx = toneTrials(:,2);
        else
            stimTrials = trials{trials.isTone == 0 & trials.isStim == 1 & trials.isReward == 1,["TrialNumber","CueTime","OutcomeTime","ENL"]};
            stimOmissionTrials = trials{trials.isTone == 0 & trials.isStim == 1 & trials.isReward == 0,["TrialNumber","CueTime","OutcomeTime","ENL"]};
            pairTrials = trials{trials.isTone == 1 & trials.isStim == 1 & trials.isReward == 1,["TrialNumber","CueTime","OutcomeTime","ENL"]};
            pairOmissionTrials = trials{trials.isTone == 1 & trials.isStim == 1 & trials.isReward == 0,["TrialNumber","CueTime","OutcomeTime","ENL"]};
            toneTrials = trials{trials.isTone == 1 & trials.isStim == 0 & trials.isReward == 1,["TrialNumber","CueTime","OutcomeTime","ENL"]};
            toneOmissionTrials = trials{trials.isTone == 1 & trials.isStim == 0 & trials.isReward == 0,["TrialNumber","CueTime","OutcomeTime","ENL"]};
            if isempty(toneTrials) && ~isempty(toneOmissionTrials)
                toneTrials = toneOmissionTrials;
            end
            stimIdx = stimTrials(:,2);
            stimOmissionIdx = stimOmissionTrials(:,2);
            pairIdx = pairTrials(:,2);
            pairOmissionIdx = pairOmissionTrials(:,2);
            toneIdx = toneTrials(:,2);
            toneOmissionIdx = toneOmissionTrials(:,2);
        end
    elseif contains(options.task,'punish')
        if options.performing
            stimTrials = trials{trials.isTone == 0 & trials.isStim == 1 ...
                & trials.isPunishment == 1 & trials.performing == 1, ["TrialNumber","CueTime","OutcomeTime","ENL"]};
            stimOmissionTrials = trials{trials.isTone == 0 & trials.isStim == 1 ...
                & trials.isPunishment == 0 & trials.performing == 1, ["TrialNumber","CueTime","OutcomeTime","ENL"]};
            pairTrials = trials{trials.isTone == 1 & trials.isStim == 1 ...
                & trials.isPunishment == 1 & trials.performing == 1, ["TrialNumber","CueTime","OutcomeTime","ENL"]};
            pairOmissionTrials = trials{trials.isTone == 1 & trials.isStim == 1 ...
                & trials.isPunishment == 0 & trials.performing == 1, ["TrialNumber","CueTime","OutcomeTime","ENL"]};
            toneTrials = trials{trials.isTone == 1 & trials.isStim == 0 ...
                & trials.performing == 1,["TrialNumber","CueTime","OutcomeTime","ENL"]};
            stimIdx = stimTrials(:,2);
            stimOmissionIdx = stimOmissionTrials(:,2);
            pairIdx = pairTrials(:,2);
            pairOmissionIdx = pairOmissionTrials(:,2);
            toneIdx = toneTrials(:,2);
        else
            stimTrials = trials{trials.isTone == 0 & trials.isStim == 1 & trials.isPunishment == 1,["TrialNumber","CueTime","OutcomeTime","ENL"]};
            stimOmissionTrials = trials{trials.isTone == 0 & trials.isStim == 1 & trials.isPunishment == 0,["TrialNumber","CueTime","OutcomeTime","ENL"]};
            pairTrials = trials{trials.isTone == 1 & trials.isStim == 1 & trials.isPunishment == 1,["TrialNumber","CueTime","OutcomeTime","ENL"]};
            pairOmissionTrials = trials{trials.isTone == 1 & trials.isStim == 1 & trials.isPunishment == 0,["TrialNumber","CueTime","OutcomeTime","ENL"]};
            toneTrials = trials{trials.isTone == 1 & trials.isStim == 0 & trials.isPunishment == 1,["TrialNumber","CueTime","OutcomeTime","ENL"]};
            toneOmissionTrials = trials{trials.isTone == 1 & trials.isStim == 0 & trials.isPunishment == 0,["TrialNumber","CueTime","OutcomeTime","ENL"]};
            if isempty(toneTrials) && ~isempty(toneOmissionTrials)
                toneTrials = toneOmissionTrials;
            end
            stimIdx = stimTrials(:,2);
            stimOmissionIdx = stimOmissionTrials(:,2);
            pairIdx = pairTrials(:,2);
            pairOmissionIdx = pairOmissionTrials(:,2);
            toneIdx = toneTrials(:,2);
            toneOmissionIdx = toneOmissionTrials(:,2);
        end
    end

    stageTime = [-2,0;0,0.5;0.5,5];
    analysisEvents = {waterIdx,waterLickIdx,toneIdx,stimIdx,pairIdx,airpuffIdx,baselineIdx};
    eventTrialNum = {findTrials(waterIdx,trials),findTrials(waterLickIdx,trials),...
                    toneTrials(:,1),stimTrials(:,1),pairTrials(:,1),...
                    findTrials(airpuffIdx,trials),findTrials(baselineIdx,trials)};
    analysisLabels = {'Water','Rewarded licks','Tone only','Stim only','Pair','Airpuff','Baseline'};
    taskLegend = getLegend(analysisEvents,analysisLabels);

    stageColors = {[.75 .75 .75],bluePurpleRed(end,:),bluePurpleRed(1,:)};
    stageLegend = {'Baseline','CS','US'};

    eventTrialNum = eventTrialNum(~cellfun('isempty',analysisEvents));
    analysisLabels = analysisLabels(~cellfun('isempty',analysisEvents));
    analysisEvents = analysisEvents(~cellfun('isempty',analysisEvents));

    for i = 1:length(analysisEvents)
        disp(['Total ',analysisLabels{i},': ',num2str(length(analysisEvents{i}))]);
    end

    % Save idx to behavior.mat
    save(strcat(sessionpath,filesep,'behavior_',params.session.name),'allTrials',...
        'waterIdx','waterLickIdx','airpuffIdx','toneIdx','stimIdx','pairIdx',...
        'baselineIdx',...
        '-append');
end

%% Save photometry/lick/eye PSTHs

if options.analyzeTraces
    analysis = analyzeTraces(timeSeries,rightLick,analysisEvents,analysisLabels,params,...
                         stageTime=stageTime,...
                         trialNumber=eventTrialNum,trialTable=trials);
end

%% Plot photometry summary plots

% Find the number of photometry channels
photometryIdx = find(cellfun(@(x) contains(x,["NI","LJ"],"IgnoreCase",true), {timeSeries.system}));
photometryName = cellfun(@(x) unique(x,'rows'), {timeSeries(photometryIdx).name},'UniformOutput',false);
nSignals = length(photometryIdx);
disp(['Finished: found ', num2str(nSignals),' photometry signals']);

if options.plotPhotometry
    %% Plot pre-processing steps
    initializeFig(.67,.67); tiledlayout('flow');
    
    for i = 1:nSignals
        path = photometryIdx(i);
        nexttile;
        histogram(normrnd(0,1,size(timeSeries(path).data)),200); hold on
        histogram(timeSeries(path).data,200); hold on
        box off

        skew_lj = skewness(timeSeries(path).data); 
        kur_lj = kurtosis(timeSeries(path).data);
        xlabel('z-score'); ylabel('Count'); legend({'Normal distribution',timeSeries(path).name});
        
        title(timeSeries(path).name);
        subtitle(strcat("Skewness: ",num2str(skew_lj),", Kurtosis: ",num2str(kur_lj)));
    end
    
    % Save figure
    saveas(gcf,strcat(sessionpath,filesep,'Summary_photometry_distribution.png'));

    %% Loop through timeSeries struct
    for photometry = 1:nSignals

        % Load signal of interest
        path = photometryIdx(photometry);
        signal = timeSeries(path).data;
        finalFs = timeSeries(path).finalFs;
        system = timeSeries(path).system;
        name = timeSeries(path).name;
        
        %% Plot combined PSTH
        timeRange = [-1,5];
        % 2. Plot traces
        % Determine figure layout 
        if params.session.withCamera && params.session.withEyeTracking
            initializeFig(0.5, 1); tiledlayout(4,1);
        else
            initializeFig(0.5,0.5); tiledlayout(2,1);
        end

        if strcmp(options.task,'random')
            % 2.1 Plot photometry traces
            nexttile
            [~,~] = plotTraces(waterLickIdx,timeRange,signal,bluePurpleRed(1,:),params,...
                        signalFs=finalFs,signalSystem=system);
            [~,~] = plotTraces(airpuffIdx,timeRange,signal,[0.2, 0.2, 0.2],params,...
                        signalFs=finalFs,signalSystem=system);
            [~,~] = plotTraces(toneIdx,timeRange,signal,bluePurpleRed(350,:),params,...
                        signalFs=finalFs,signalSystem=system);
            [~,~] = plotTraces(stimIdx,timeRange,signal,bluePurpleRed(end,:),params,...
                        signalFs=finalFs,signalSystem=system);
            [~,~] = plotTraces(baselineIdx,timeRange,signal,[.75 .75 .75],params,...
                        signalFs=finalFs,signalSystem=system);
            plotEvent('',0);
            xlabel('Time (s)'); ylabel([name,' z-score']);
            legend(taskLegend(2:end),'Location','northeast');
            
            % 2.2 Plot lick traces
            nexttile
            plotLicks(waterLickIdx,timeRange,options.lick_binSize,bluePurpleRed(1,:),[],rightLick,params);
            plotLicks(airpuffIdx,timeRange,options.lick_binSize,[0.2, 0.2, 0.2],[],rightLick,params);
            plotLicks(toneIdx,timeRange,options.lick_binSize,bluePurpleRed(350,:),[],rightLick,params);
            plotLicks(stimIdx,timeRange,options.lick_binSize,bluePurpleRed(end,:),[],rightLick,params);
            plotLicks(baselineIdx,timeRange,options.lick_binSize,[.75 .75 .75],[],rightLick,params);
            plotEvent('',0);
            xlabel('Time (s)'); ylabel('Licks/s'); 
            legend(taskLegend(2:end),'Location','best');

            if params.session.withCamera && params.session.withEyeTracking
                % 2.3 Plot eye area traces
                nexttile
                [~,~] = plotTraces(waterLickIdx,timeRange,eyeArea_detrend,bluePurpleRed(1,:),params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(airpuffIdx,timeRange,eyeArea_detrend,[0.2, 0.2, 0.2],params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(toneIdx,timeRange,eyeArea_detrend,bluePurpleRed(350,:),params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(stimIdx,timeRange,eyeArea_detrend,bluePurpleRed(end,:),params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(baselineIdx,timeRange,eyeArea_detrend,[.75 .75 .75],params,signalSystem='camera',smooth=15);
                plotEvent('',0);
                xlabel('Time (s)'); ylabel('Eye area (z-score)');
                legend(taskLegend(2:end),'Location','northeast');
    
                % 2.4 Plot pupil area traces
                nexttile
                [~,~] = plotTraces(waterLickIdx,timeRange,pupilArea_detrend,bluePurpleRed(1,:),params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(airpuffIdx,timeRange,pupilArea_detrend,[0.2, 0.2, 0.2],params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(toneIdx,timeRange,pupilArea_detrend,bluePurpleRed(350,:),params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(stimIdx,timeRange,pupilArea_detrend,bluePurpleRed(end,:),params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(baselineIdx,timeRange,pupilArea_detrend,[.75 .75 .75],params,signalSystem='camera',smooth=15);
                plotEvent('',0);
                xlabel('Time (s)'); ylabel('Pupil area (z-score)');
                legend(taskLegend(2:end),'Location','northeast');
            end
    
        else
            % 2.1 Plot photometry traces
            nexttile
            [~,~] = plotTraces(waterLickIdx,timeRange,signal,bluePurpleRed(1,:),params,...
                        signalFs=finalFs,signalSystem=system);
            [~,~] = plotTraces(toneIdx,timeRange,signal,bluePurpleRed(350,:),params,...
                        signalFs=finalFs,signalSystem=system);
            [~,~] = plotTraces(stimIdx,timeRange,signal,bluePurpleRed(end,:),params,...
                        signalFs=finalFs,signalSystem=system);
            [~,~] = plotTraces(pairIdx,timeRange,signal,bluePurpleRed(150,:),params,...
                        signalFs=finalFs,signalSystem=system);
            [~,~] = plotTraces(airpuffIdx,timeRange,signal,[0.2, 0.2, 0.2],params,...
                        signalFs=finalFs,signalSystem=system);
            [~,~] = plotTraces(baselineIdx,timeRange,signal,[.75 .75 .75],params,...
                        signalFs=finalFs,signalSystem=system);
            plotEvent('',0);
            xlabel('Time (s)'); ylabel([name,' z-score']);
            legend(taskLegend(2:end),'Location','northeast');
            
            % 2.2 Plot lick traces
            nexttile
            plotLicks(waterLickIdx,timeRange,options.lick_binSize,bluePurpleRed(1,:),[],rightLick,params);
            plotLicks(toneIdx,timeRange,options.lick_binSize,bluePurpleRed(350,:),[],rightLick,params);
            plotLicks(stimIdx,timeRange,options.lick_binSize,bluePurpleRed(end,:),[],rightLick,params);
            plotLicks(pairIdx,timeRange,options.lick_binSize,bluePurpleRed(150,:),[],rightLick,params);
            plotLicks(airpuffIdx,timeRange,options.lick_binSize,[0.2, 0.2, 0.2],[],rightLick,params);
            plotLicks(baselineIdx,timeRange,options.lick_binSize,[.75 .75 .75],[],rightLick,params);
            plotEvent('',0);
            xlabel('Time (s)'); ylabel('Licks/s'); 
            legend(taskLegend(2:end),'Location','best');

            if params.session.withCamera && params.session.withEyeTracking
                % 2.3 Plot eye area traces
                nexttile
                [~,~] = plotTraces(waterLickIdx,timeRange,eyeArea_detrend,bluePurpleRed(1,:),params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(toneIdx,timeRange,eyeArea_detrend,bluePurpleRed(350,:),params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(stimIdx,timeRange,eyeArea_detrend,bluePurpleRed(end,:),params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(pairIdx,timeRange,eyeArea_detrend,bluePurpleRed(150,:),params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(airpuffIdx,timeRange,eyeArea_detrend,[0.2, 0.2, 0.2],params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(baselineIdx,timeRange,eyeArea_detrend,[.75 .75 .75],params,signalSystem='camera',smooth=15);
                plotEvent('',0);
                xlabel('Time (s)'); ylabel('Eye area (z-score)');
                legend(taskLegend(2:end),'Location','northeast');
    
                % 2.4 Plot pupil area traces
                nexttile
                [~,~] = plotTraces(waterLickIdx,timeRange,pupilArea_detrend,bluePurpleRed(1,:),params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(toneIdx,timeRange,pupilArea_detrend,bluePurpleRed(350,:),params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(stimIdx,timeRange,pupilArea_detrend,bluePurpleRed(end,:),params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(pairIdx,timeRange,pupilArea_detrend,bluePurpleRed(150,:),params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(airpuffIdx,timeRange,pupilArea_detrend,[0.2, 0.2, 0.2],params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(baselineIdx,timeRange,pupilArea_detrend,[.75 .75 .75],params,signalSystem='camera',smooth=15);
                plotEvent('',0);
                xlabel('Time (s)'); ylabel('Pupil area (z-score)');
                legend(taskLegend(2:end),'Location','northeast');
            end
        end 
        saveas(gcf,strcat(sessionpath,filesep,'Summary_events_',timeSeries(path).name,'.png'));

        %% Plot single stimulus PSTH
        if strcmp(options.task,'random')
            eventIdxes = {stimIdx,waterLickIdx,toneIdx,airpuffIdx};
            labels = {'Stim','Water','Tone','Airpuff'};
            eventDurations = [0.5,0,0.5,0.02];
            groupSizes = [20,30,10,30];
            longTimeRange = [-5,10];
            shortTimeRange = [-1,5]; 
            
            for event = 1:length(eventIdxes)
                eventIdx = eventIdxes{event};
                if isempty(eventIdx); continue; end
                label = labels{event}; eventDuration = eventDurations(event); 
                groupSize = groupSizes(event); % num of trials to plot in one line
                
                initializeFig(0.67,1); tiledlayout(4,4);
                
                % Plot short timescale
                nexttile(3,[1 2]);
                [traces,t] = plotTraces(eventIdx,shortTimeRange,signal,bluePurpleRed(1,:),params,...
                                signalFs=finalFs,signalSystem=system,...
                                plotIndividual=true);
                plotEvent(label,eventDuration); 
                xlabel('Time (s)'); ylabel([name,' z-score']);
                legend({['Shuffled (n=',num2str(length(eventIdx)),')'],...
                        [label,' (n=',num2str(length(eventIdx)),')']},...
                        'Location','northeast');
                
                nexttile(7,[1 2]);
                legendList = plotGroupTraces(traces,t,bluePurpleRed,groupSize=groupSize);
                plotEvent(label,eventDuration);
                xlabel('Time (s)'); ylabel([name,' z-score']);
                legend(legendList);

                % Plot heatmap
                nexttile(1,[4 2]);
                plotHeatmap(traces,t);
                plotEvent(label,eventDuration);
                xlabel('Time (s)'); ylabel('Trials');


                % Plot long timescale
                nexttile(11,[1 2]);
                [traces,t] = plotTraces(eventIdx,longTimeRange,signal,bluePurpleRed(1,:),params,...
                                signalFs=finalFs,signalSystem=system,...
                                plotIndividual=true);
                plotEvent(label,eventDuration);
                xlabel('Time (s)'); ylabel([name,' z-score']);
                legend({['Shuffled (n=',num2str(length(eventIdx)),')'],...
                        [label,' (n=',num2str(length(eventIdx)),')']},...
                        'Location','northeast');
                
                nexttile(15,[1 2]);
                legendList = plotGroupTraces(traces,t,bluePurpleRed,groupSize=groupSize);
                plotEvent(label,eventDuration);
                xlabel('Time (s)'); ylabel([name,' z-score']);
                legend(legendList);
                
                saveas(gcf,strcat(sessionpath,filesep,'Events_',timeSeries(path).name,'_',label,'.png'));
            end
        else
            eventIdxes = {stimIdx,pairIdx,toneIdx,waterLickIdx,airpuffIdx};
            omissionIdxes = {stimOmissionIdx, pairOmissionIdx,toneOmissionIdx,[],[]};
            labels = {'Stim','Pair','Tone','Water','Airpuff'};
            eventDurations = [0.5,0.5,0.5,0,0.02];
            groupSizes = [10,10,10,30,30];
            longTimeRange = [-5,10];
            shortTimeRange = [-1,5]; 
            
            for event = 1:length(eventIdxes)
                eventIdx = eventIdxes{event};
                omissionIdx = omissionIdxes{event};
                if isempty(eventIdx); continue; end
                label = labels{event}; eventDuration = eventDurations(event); 
                groupSize = groupSizes(event); % num of trials to plot in one line
                
                initializeFig(0.67,1);
                tiledlayout(4,4);

                % Plot short time scale
                nexttile(3,[1 2]);
                [traces,t] = plotTraces(eventIdx,shortTimeRange,signal,bluePurpleRed(1,:),params,...
                                signalFs=finalFs,signalSystem=system,...
                                plotIndividual=true);
                [~,~] = plotTraces(omissionIdx,shortTimeRange,signal,[0.3, 0.3, 0.3],params,...
                                signalFs=finalFs,signalSystem=system);
                plotEvent(label,eventDuration);
                xlabel('Time (s)'); ylabel([name,' z-score']);
                legend({['Shuffled (n=',num2str(length(eventIdx)),')'],...
                        [label,' (n=',num2str(length(eventIdx)),')'],...
                        [label,' omission (n=',num2str(length(omissionIdxes)),')']},...
                        'Location','northeast');
                
                nexttile(7,[1 2]);
                legendList = plotGroupTraces(traces,t,bluePurpleRed,groupSize=groupSize);
                plotEvent(label,eventDuration);
                xlabel('Time (s)'); ylabel([name,' z-score']);
                legend(legendList);

                % Plot heatmap
                nexttile(1,[4 2]);
                plotHeatmap(traces,t);
                plotEvent(label,eventDuration);
                xlabel('Time (s)'); ylabel('Trials');

                % Plot long time scale
                nexttile(11,[1 2]);
                [traces,t] = plotTraces(eventIdx,longTimeRange,signal,bluePurpleRed(1,:),params,...
                                signalFs=finalFs,signalSystem=system,...
                                plotIndividual=true);
                [~,~] = plotTraces(omissionIdx,longTimeRange,signal,[0.2, 0.2, 0.2],params,...
                                signalFs=finalFs,signalSystem=system);
                plotEvent(label,eventDuration);
                xlabel('Time (s)'); ylabel([name,' z-score']);
                legend({['Shuffled (n=',num2str(length(eventIdx)),')'],...
                        [label,' (n=',num2str(length(eventIdx)),')'],...
                        [label,' omission (n=',num2str(length(omissionIdxes)),')']},...
                        'Location','northeast');
                
                nexttile(15,[1 2]);
                legendList = plotGroupTraces(traces,t,bluePurpleRed,groupSize=groupSize);
                plotEvent(label,eventDuration);
                xlabel('Time (s)'); ylabel([name,' z-score']);
                legend(legendList);

                saveas(gcf,strcat(sessionpath,filesep,'Events_',timeSeries(path).name,'_',label,'.png'));
            end
        end
    end

    if options.analyzeTraces
        %% Plot subtrial average trend for baseline, CS, US
        nBins = 500;
        for signal = 1:nSignals
            % Load current signal name and row in analysis
            cur_signal = photometryName{signal};
            signalIdx = find(cellfun(@(x) strcmpi(x,cur_signal), {analysis.name}));
            eventPlotted = 0;
            % Plot analysis result for this signal
            initializeFig(0.8,0.8); tiledlayout(2,length(analysisEvents));
            for i = 1:length(signalIdx)
                row = signalIdx(i); eventPlotted = eventPlotted + 1;
                nexttile;
                for stage = 1:size(stageTime,1)
                    data = analysis(row).stageAvg.data(:,stage);
                    p = analysis(row).stageAvg.fit(stage,:);
                    x = 1:length(data);
                    scatter(x,data',100,stageColors{stage},'filled',MarkerFaceAlpha=0.5,HandleVisibility='off'); hold on
                    plot(x,polyval(p,x),Color=stageColors{stage},lineWidth=5);
                end
                title(analysis(row).event);
                xlabel('Trials'); ylabel([analysis(row).name,' signal average (z-score)']);
                legend(stageLegend); box off
    
                nexttile(eventPlotted+length(analysisEvents));
                for stage = 1:size(stageTime,1)
                    h = histogram(analysis(row).stageAvg.stats.bs(stage,:,1),nBins); hold on
                    h.FaceColor = stageColors{stage}; h.EdgeColor = stageColors{stage};
                    h.FaceAlpha = 0.25; h.EdgeAlpha = 0.25;
                    pval = min(analysis(row).stageAvg.stats.pval_slope(stage,:));
                    xline(analysis(row).stageAvg.fit(stage,1),'-',...
                        {['Slope (',stageLegend{stage},')'],['p=',num2str(pval)]},...
                        'Color',stageColors{stage},...
                        'LineWidth',3,...
                        'LabelOrientation','horizontal');
                end
                xlabel('Slope distribution (bootstrapped)'); ylabel('Count'); box off
            end
            % Save
            saveas(gcf,strcat(sessionpath,filesep,'Analysis_',cur_signal,'_subtrial_average.png'));
        end
    
        %% Plot subtrial peak trend for baseline, CS, US
        for signal = 1:nSignals
            % Load current signal name and row in analysis
            cur_signal = photometryName{signal};
            signalIdx = find(cellfun(@(x) strcmpi(x,cur_signal), {analysis.name}));
            eventPlotted = 0;
            % Plot analysis result for this signal
            initializeFig(0.8,0.8); tiledlayout(2,length(analysisEvents));
            for i = 1:length(signalIdx)
                row = signalIdx(i); eventPlotted = eventPlotted + 1;
                nexttile;
                for stage = 1:size(stageTime,1)
                    data = analysis(row).stageMax.data(:,stage);
                    p = analysis(row).stageMax.fit(stage,:);
                    x = 1:length(data);
                    scatter(x,data',100,stageColors{stage},'filled',MarkerFaceAlpha=0.5,HandleVisibility='off'); hold on
                    plot(x,polyval(p,x),Color=stageColors{stage},lineWidth=5);
                end
                title(analysis(row).event);
                xlabel('Trials'); ylabel([analysis(row).name,' signal peak (z-score)']);
                legend(stageLegend); box off
    
                nexttile(eventPlotted+length(analysisEvents));
                for stage = 1:size(stageTime,1)
                    h = histogram(analysis(row).stageMax.stats.bs(stage,:,1),nBins); hold on
                    h.FaceColor = stageColors{stage}; h.EdgeColor = stageColors{stage};
                    h.FaceAlpha = 0.25; h.EdgeAlpha = 0.25;
                    pval = min(analysis(row).stageMax.stats.pval_slope(stage,:));
                    xline(analysis(row).stageMax.fit(stage,1),'-',...
                        {['Slope (',stageLegend{stage},')'],['p=',num2str(pval)]},...
                        'Color',stageColors{stage},...
                        'LineWidth',3,...
                        'LabelOrientation','horizontal');
                end
                xlabel('Slope distribution (bootstrapped)'); ylabel('Count'); box off
            end
            % Save
            saveas(gcf,strcat(sessionpath,filesep,'Analysis_',cur_signal,'_subtrial_peak.png'));
        end
    
        %% Plot subtrial trough trend for baseline, CS, US
        for signal = 1:nSignals
            % Load current signal name and row in analysis
            cur_signal = photometryName{signal};
            signalIdx = find(cellfun(@(x) strcmpi(x,cur_signal), {analysis.name}));
            eventPlotted = 0;
            % Plot analysis result for this signal
            initializeFig(0.8,0.8); tiledlayout(2,length(analysisEvents));
            for i = 1:length(signalIdx)
                row = signalIdx(i); eventPlotted = eventPlotted + 1;
                nexttile;
                for stage = 1:size(stageTime,1)
                    data = analysis(row).stageMin.data(:,stage);
                    p = analysis(row).stageMin.fit(stage,:);
                    x = 1:length(data);
                    scatter(x,data',100,stageColors{stage},'filled',MarkerFaceAlpha=0.5,HandleVisibility='off'); hold on
                    plot(x,polyval(p,x),Color=stageColors{stage},lineWidth=5);
                end
                title(analysis(row).event);
                xlabel('Trials'); ylabel([analysis(row).name,' signal trough (z-score)']);
                legend(stageLegend); box off
    
                nexttile(eventPlotted+length(analysisEvents));
                for stage = 1:size(stageTime,1)
                    h = histogram(analysis(row).stageMin.stats.bs(stage,:,1),nBins); hold on
                    h.FaceColor = stageColors{stage}; h.EdgeColor = stageColors{stage};
                    h.FaceAlpha = 0.25; h.EdgeAlpha = 0.25;
                    pval = min(analysis(row).stageMin.stats.pval_slope(stage,:));
                    xline(analysis(row).stageMin.fit(stage,1),'-',...
                        {['Slope (',stageLegend{stage},')'],['p=',num2str(pval)]},...
                        'Color',stageColors{stage},...
                        'LineWidth',3,...
                        'LabelOrientation','horizontal');
                end
                xlabel('Slope distribution (bootstrapped)'); ylabel('Count'); box off
            end
            % Save
            saveas(gcf,strcat(sessionpath,filesep,'Analysis_',cur_signal,'_subtrial_trough.png'));
        end
    end
end

%% Plot behavior related plots

% Plot lick bout distribution
if options.plotBehavior
    initializeFig(0.5,0.5); tiledlayout(2,2);

    ENLinSec = trials{1:end-1,"ENL"} / params.sync.behaviorFs;
    ITIextra = trials{1:end-1,'ITI'} - ENLinSec;

    % Plot lick bout count distribution
    nexttile;
    histogram(lickBout(:,2),30); 
    xlabel('Licks per lick bout'); ylabel('Count'); box off

    % Plot ITI-ENL distribution
    nexttile;
    histogram(trials{1:end-1,'ITI'},30); 
    xlabel('ITI (s)'); ylabel('Count'); box off

    % Plot lick per trial vs trials
    nexttile;
    plot(trials{1:end-1,'TrialNumber'},trials{1:end-1,'nLicks'},'Color',bluePurpleRed(1,:),LineWidth=2);
    xlabel('Trials'); ylabel('Licks per trial'); box off

    % Plot ITI-ENL vs trials
    nexttile;
    plot(trials{1:end-1,'TrialNumber'},ENLinSec,'Color',bluePurpleRed(500,:),LineWidth=2); hold on
    plot(trials{1:end-1,'TrialNumber'},ITIextra,'Color',bluePurpleRed(1,:),LineWidth=2);
    xlabel('Trials'); ylabel('Time (s)'); ylim([0,Inf]); box off
    legend({'ENL','ITI-ENL'});

    saveas(gcf,strcat(sessionpath,filesep,'Behavior_ITI&LickBout.png'));
end


if options.plotBehavior && contains(options.task,'pairing')     
    %% Plot session overview for licking
    timeRange = [-5,10]; cameraTimeRange = [-1,5];
    markerSize = 20;

    % Get event time and number by trial type
    stimTrials(:,3) = stimTrials(:,3)./params.sync.behaviorFs;
    toneTrials(:,3) = toneTrials(:,3)./params.sync.behaviorFs;
    pairTrials(:,3) = pairTrials(:,3)./params.sync.behaviorFs;

    % getLicks by trial type
    [stimLickRate,~,stimLicks] = getLicks(timeRange,stimIdx,options.lick_binSize,[],rightLick,...
                                params.sync.behaviorFs,params.sync.timeNI);
    [pairLickRate,~,pairLicks] = getLicks(timeRange,pairIdx,options.lick_binSize,[],rightLick,...
                                params.sync.behaviorFs,params.sync.timeNI);
    [toneLickRate,~,toneLicks] = getLicks(timeRange,toneIdx,options.lick_binSize,[],rightLick,...
                                params.sync.behaviorFs,params.sync.timeNI);

    % 1. Plot lick raster and trace
    initializeFig(0.67,0.67); tiledlayout(3,2);
    % 1.1 Plot lick raster plot for session
    nexttile([3,1]);
    for i = 1:size(stimLicks,1)
        scatter(stimLicks{i},stimTrials(i,1),markerSize,'filled','MarkerFaceColor',bluePurpleRed(end,:)); hold on
        scatter(stimTrials(i,3),stimTrials(i,1),markerSize+10,bluePurpleRed(1,:),'pentagram','filled'); hold on
    end
    for i = 1:size(toneLicks,1)
        scatter(toneLicks{i},toneTrials(i,1),markerSize,'filled','MarkerFaceColor',bluePurpleRed(350,:)); hold on
        scatter(toneTrials(i,3),toneTrials(i,1),markerSize+10,bluePurpleRed(1,:),'pentagram','filled'); hold on
    end
    for i = 1:size(pairLicks,1)
        scatter(pairLicks{i},pairTrials(i,1),markerSize,'filled','MarkerFaceColor',bluePurpleRed(150,:)); hold on
        scatter(pairTrials(i,3),pairTrials(i,1),markerSize+10,bluePurpleRed(1,:),'pentagram','filled'); hold on
    end
    xlabel('Time (s)'); xlim([timeRange(1),timeRange(2)]);
    ylabel('Trial'); ylim([0,size(trials,1)]);
    plotEvent("",0.5);
    % 1.2 Plot lick traces across session
    traces = {stimLickRate,pairLickRate,toneLickRate};
    labels = {'Stim','Pair','Tone'};
    groupSizes = [20, 20, 10];
    for event = 1:length(traces)
        trace = traces{event};
        label = labels{event};
        groupSize = groupSizes(event);

        nexttile;
        t = linspace(timeRange(1),timeRange(2),size(trace,2));
        nLines = ceil(size(trace,1)/groupSize); legendList = cell(nLines,1);
        nColors = round(linspace(1,size(bluePurpleRed,1),nLines));
        for i = 1:nLines
            startTrial = (i-1)*groupSize+1; 
            if i == nLines; endTrial = size(trace,1);
            else; endTrial = i*groupSize; end
            plotSEM(t,trace(startTrial:endTrial,:),bluePurpleRed(nColors(i),:));
            legendList{i} = ['Trial ', num2str(startTrial),'-',num2str(endTrial)];
        end
        plotEvent(label,0.5);
        xlabel('Time (s)'); ylabel('Licks/s');
        legend(legendList);
    end
    saveas(gcf,strcat(sessionpath,filesep,'Behavior_LickOverview.png'));

    
    %% Plot session overview for eye
    if params.session.withCamera && params.session.withEyeTracking
        initializeFig(0.67,0.67); tiledlayout(3,4); 

        % Plot eye trace heatmap
        nexttile([3 2]);
        [allEyeTraces,t_cam] = plotTraces(trials{1:end-1,"CueTime"},cameraTimeRange,eyeArea,[0,0,0],params,plot=false,signalSystem='camera'); 
        imagesc(t_cam,1:size(trials,1),allEyeTraces);
        set(gca,'YDir','normal');
        colorbar; box off
        plotEvent('Cue',0.5);
        xlabel('Time (s)'); ylabel('Trials');

        % get camera trace by trial type
        [stimEyeArea,t_cam] = plotTraces(stimIdx,cameraTimeRange,eyeArea_detrend,[0,0,0],params,plot=false,signalSystem='camera');
        [pairEyeArea,~] = plotTraces(pairIdx,cameraTimeRange,eyeArea_detrend,[0,0,0],params,plot=false,signalSystem='camera');
        [toneEyeArea,~] = plotTraces(toneIdx,cameraTimeRange,eyeArea_detrend,[0,0,0],params,plot=false,signalSystem='camera');
        [stimPupilArea,~] = plotTraces(stimIdx,cameraTimeRange,pupilArea_detrend,[0,0,0],params,plot=false,signalSystem='camera');
        [pairPupilArea,~] = plotTraces(pairIdx,cameraTimeRange,pupilArea_detrend,[0,0,0],params,plot=false,signalSystem='camera');
        [tonePupilArea,~] = plotTraces(toneIdx,cameraTimeRange,pupilArea_detrend,[0,0,0],params,plot=false,signalSystem='camera');

        % Plot eye area across session
        traces = {stimEyeArea,pairEyeArea,toneEyeArea};
        labels = {'Stim','Pair','Tone'};
        groupSizes = [20, 20, 10];
        for event = 1:length(traces)
            trace = traces{event};
            label = labels{event};
            groupSize = groupSizes(event);
    
            nexttile(3+ 4*(event-1));
            nLines = ceil(size(trace,1)/groupSize); legendList = cell(nLines,1);
            nColors = round(linspace(1,size(bluePurpleRed,1),nLines));
            for i = 1:nLines
                startTrial = (i-1)*groupSize+1; 
                if i == nLines; endTrial = size(trace,1);
                else; endTrial = i*groupSize; end
                plotSEM(t_cam,trace(startTrial:endTrial,:),bluePurpleRed(nColors(i),:),smooth=15);
                legendList{i} = ['Trial ', num2str(startTrial),'-',num2str(endTrial)];
            end
            plotEvent(label,0.5);
            xlabel('Time (s)'); ylabel('Eye area (z-score)');
            legend(legendList);
        end
    
        % Plot pupil area across session
        traces = {stimPupilArea,pairPupilArea,tonePupilArea};
        labels = {'Stim','Pair','Tone'};
        groupSizes = [20, 20, 10];
        for event = 1:length(traces)
            trace = traces{event};
            label = labels{event};
            groupSize = groupSizes(event);
    
            nexttile(4+ 4*(event-1));
            nLines = ceil(size(trace,1)/groupSize); legendList = cell(nLines,1);
            nColors = round(linspace(1,size(bluePurpleRed,1),nLines));
            for i = 1:nLines
                startTrial = (i-1)*groupSize+1; 
                if i == nLines; endTrial = size(trace,1);
                else; endTrial = i*groupSize; end
                plotSEM(t_cam,trace(startTrial:endTrial,:),bluePurpleRed(nColors(i),:),smooth=15);
                legendList{i} = ['Trial ', num2str(startTrial),'-',num2str(endTrial)];
            end
            plotEvent(label,0.5);
            xlabel('Time (s)'); ylabel('Pupil area (z-score)');
            legend(legendList);
        end
        saveas(gcf,strcat(sessionpath,filesep,'Behavior_EyeOverview.png'));
    end

    %% Plot ENL aligned lick trace
    % Bin trials based on ENL length, plot lick trace aligning to ENL start
    % should not smear if animal lick specifically to the cue

    initializeFig(0.67,0.67); tiledlayout(3,3);
    timeRangeENL = [0,10];
    [stimLickRateENL,~,~] = getLicks(timeRangeENL,stimIdx-stimTrials(:,4),options.lick_binSize,[],rightLick,...
                                    params.sync.behaviorFs,params.sync.timeNI);
    [pairLickRateENL,~,~] = getLicks(timeRangeENL,pairIdx-pairTrials(:,4),options.lick_binSize,[],rightLick,...
                                params.sync.behaviorFs,params.sync.timeNI);
    [toneLickRateENL,~,~] = getLicks(timeRangeENL,toneIdx-toneTrials(:,4),options.lick_binSize,[],rightLick,...
                                params.sync.behaviorFs,params.sync.timeNI);

    nENLBins = [2,3,4];
    for i = 1:length(nENLBins)

        % 1.1 Sort trials based on ENL duration
        % For stim only
        nexttile;
        trialIdx = stimTrials; traces = stimLickRateENL;
        [~,sortIdx] = sortrows(trialIdx,4); 
        nModTrials = mod(height(trialIdx),nENLBins(i));
        ENLsortIdx = mat2cell(reshape(sortIdx(1:end-nModTrials),[],nENLBins(i))',ones(1,nENLBins(i)));
        ENLsortIdx{end} = [ENLsortIdx{end}, sortIdx(end-nModTrials+1:end)'];

        t = linspace(timeRangeENL(1),timeRangeENL(2),size(traces,2));
        legendList = cell(nENLBins(i),1);
        nColors = round(linspace(1,size(bluePurpleRed,1),nENLBins(i)));
        % 1.2 Plot lick trace
        for j = 1:nENLBins(i)
            plotSEM(t,traces(ENLsortIdx{j},:),bluePurpleRed(nColors(j),:));
            legendList{j} = ['ENL bin ', num2str(j),' (n=',num2str(length(ENLsortIdx{j})),')'];
        end
        xlabel('Time from ENL start (s)'); ylabel('Licks/s');
        legend(legendList);

        % For pair
        nexttile;
        trialIdx = pairTrials; traces = pairLickRateENL;
        [~,sortIdx] = sortrows(trialIdx,4); 
        nModTrials = mod(height(trialIdx),nENLBins(i));
        ENLsortIdx = mat2cell(reshape(sortIdx(1:end-nModTrials),[],nENLBins(i))',ones(1,nENLBins(i)));
        ENLsortIdx{end} = [ENLsortIdx{end}, sortIdx(end-nModTrials+1:end)'];

        t = linspace(timeRangeENL(1),timeRangeENL(2),size(traces,2));
        legendList = cell(nENLBins(i),1);
        nColors = round(linspace(1,size(bluePurpleRed,1),nENLBins(i)));
        % 1.2 Plot lick trace
        for j = 1:nENLBins(i)
            plotSEM(t,traces(ENLsortIdx{j},:),bluePurpleRed(nColors(j),:));
            legendList{j} = ['ENL bin ', num2str(j),' (n=',num2str(length(ENLsortIdx{j})),')'];
        end
        xlabel('Time from ENL start (s)'); ylabel('Licks/s');
        legend(legendList);

        % For tone only
        nexttile;
        trialIdx = toneTrials; traces = toneLickRateENL;
        [~,sortIdx] = sortrows(trialIdx,4); 
        nModTrials = mod(height(trialIdx),nENLBins(i));
        ENLsortIdx = mat2cell(reshape(sortIdx(1:end-nModTrials),[],nENLBins(i))',ones(1,nENLBins(i)));
        ENLsortIdx{end} = [ENLsortIdx{end}, sortIdx(end-nModTrials+1:end)'];

        t = linspace(timeRangeENL(1),timeRangeENL(2),size(traces,2));
        legendList = cell(nENLBins(i),1);
        nColors = round(linspace(1,size(bluePurpleRed,1),nENLBins(i)));
        % 1.2 Plot lick trace
        for j = 1:nENLBins(i)
            plotSEM(t,traces(ENLsortIdx{j},:),bluePurpleRed(nColors(j),:));
            legendList{j} = ['ENL bin ', num2str(j),' (n=',num2str(length(ENLsortIdx{j})),')'];
        end
        xlabel('Time from ENL start (s)'); ylabel('Licks/s');
        legend(legendList);
    end
    saveas(gcf,strcat(sessionpath,filesep,'Behavior_ENLAlignedLick.png'));

    %% Plot distributions
    initializeFig(1,1);
    tiledlayout(4,4);

    % 1. Distribution of baseline licks
    % Find out baseline period (eg 5s-10s after reward), bootstrap and
    % compare with stim onset
    % (i.e.) how likely I get 2+ licks within a 2s window randomly selected
    % from baseline (a window of 10 secs)?
    timeRangeBaseline = [-10,0]; timeRangeReaction = [0,options.reactionTime];
    nboot = 100000; nBins = 50;
    [~,allLicksBaseline,~] = getLicks(timeRangeBaseline,trials{1:end-1,'CueTime'}-trials{1:end-1,"ENL"},...
        options.lick_binSize,[],rightLick,params.sync.behaviorFs,params.sync.timeNI);
    [~,reactionLicks,~] = getLicks(timeRangeReaction,pairIdx,...
        options.lick_binSize,[],rightLick,params.sync.behaviorFs,params.sync.timeNI);
    al_pair = sum(reactionLicks,2);
    [~,reactionLicks,~] = getLicks(timeRangeReaction,stimIdx,...
        options.lick_binSize,[],rightLick,params.sync.behaviorFs,params.sync.timeNI);
    al_stim = sum(reactionLicks,2);
    [~,reactionLicks,~] = getLicks(timeRangeReaction,toneIdx,...
        options.lick_binSize,[],rightLick,params.sync.behaviorFs,params.sync.timeNI);
    al_tone = sum(reactionLicks,2);
    if options.pavlovian
        % Bootstrap baseline licking
        lickCountBaseline = zeros(size(trials,1),nboot);
        for i = 1:nboot
            % Randomly select timepoints (columns) correspond to
            % reactionTime window
            windowInBin = round(options.reactionTime / options.lick_binSize);
            baselineSamples = datasample(allLicksBaseline,windowInBin,2,'Replace',true); 
            lickCountBaseline(:,i) = sum(baselineSamples,2);
        end
        % Plot pair trials anticipatory distribution
        nexttile; events = al_pair; nboot_event = round(size(trials,1)*nboot/size(events,1));
        h = histogram(lickCountBaseline,nBins); 
        h.FaceColor = [0.75,0.75,0.75]; h.EdgeColor = [0.75,0.75,0.75]; hold on
        [~,bootsam] = bootstrp(nboot_event,[],events);
        h = histogram(events(bootsam),nBins);
        h.FaceColor = bluePurpleRed(150,:); h.EdgeColor = bluePurpleRed(150,:);
        % xline(options.minLicks,'-.','Big reward cutoff','LineWidth',3,'LabelOrientation','horizontal');
        xline(mean(events),'-r',{'Averge anticipatory licks','(pair trials)'},'LineWidth',3,'LabelOrientation','horizontal');
        xlabel('Anticipatory licks'); ylabel ('Count'); box off
        title("Baseline licks vs anticipatory licks (pair trials)");
        % Plot stim trials anticipatory distribution
        nexttile; events = al_stim; nboot_event = round(size(trials,1)*nboot/size(events,1));
        h = histogram(lickCountBaseline,nBins); 
        h.FaceColor = [0.75,0.75,0.75]; h.EdgeColor = [0.75,0.75,0.75]; hold on
        [~,bootsam] = bootstrp(nboot_event,[],events);
        h = histogram(events(bootsam),nBins);
        h.FaceColor = bluePurpleRed(end,:); h.EdgeColor = bluePurpleRed(end,:);
        % xline(options.minLicks,'-.','Big reward cutoff','LineWidth',3,'LabelOrientation','horizontal');
        xline(mean(events),'-r',{'Averge anticipatory licks','(stim only trials)'},'LineWidth',3,'LabelOrientation','horizontal');
        xlabel('Anticipatory licks'); ylabel ('Count'); box off
        title("Baseline licks vs anticipatory licks (stim trials)");
        % Plot tone trials anticipatory distribution
        nexttile; events = al_tone; nboot_event = round(size(trials,1)*nboot/size(events,1));
        h = histogram(lickCountBaseline,nBins); 
        h.FaceColor = [0.75,0.75,0.75]; h.EdgeColor = [0.75,0.75,0.75]; hold on
        [~,bootsam] = bootstrp(nboot_event,[],events);
        h = histogram(events(bootsam),nBins);
        h.FaceColor = bluePurpleRed(350,:); h.EdgeColor = bluePurpleRed(350,:);
        % xline(options.minLicks,'-.','Big reward cutoff','LineWidth',3,'LabelOrientation','horizontal');
        xline(mean(events),'-r',{'Averge anticipatory licks','(tone only trials)'},'LineWidth',3,'LabelOrientation','horizontal');
        xlabel('Anticipatory licks'); ylabel ('Count'); box off
        title("Baseline licks vs anticipatory licks (tone trials)");

    else
        % Bootstrap baseline licking success rate
        hitPercentBaseline = zeros(nboot,1);
        for i = 1:nboot
            % Randomly select timepoints (columns) correspond to
            % reactionTime window
            windowInBin = round(options.reactionTime / options.lick_binSize);
            baselineSamples = datasample(allLicksBaseline,windowInBin,2,'Replace',true); 
            hitPercentBaseline(i) = sum(sum(baselineSamples,2) >= options.minLicks)/size(baselineSamples,1);
        end
    
        % For pair
        % Calculate cue window success rate
        hitPercentCue = length(find(trials.isTone == 1 & trials.isStim == 1 & trials.nAnticipatoryLicks >= options.minLicks))/size(pairIdx,1);
        pval = sum(hitPercentBaseline>=hitPercentCue)/nboot;
        % Plot distribution
        nexttile; 
        h = histogram(hitPercentBaseline,nBins); 
        h.FaceColor = bluePurpleRed(150,:); h.EdgeColor = bluePurpleRed(150,:);
        hold on
        xline(hitPercentCue,'-r',{'Hit rate',['p=',num2str(pval)]},...
            'LineWidth',3,...
            'LabelOrientation','horizontal');
        xlim([0,1]); xlabel('Hit rate'); ylabel ('Count'); box off
        title("Baseline licks (pair trials)");
    
        % For stim only
        % Calculate cue window success rate
        hitPercentCue = length(find(trials.isTone == 0 & trials.isStim == 1 & trials.nAnticipatoryLicks >= options.minLicks))/size(stimIdx,1);
        pval = sum(hitPercentBaseline>=hitPercentCue)/nboot;
        % Plot distribution
        nexttile; 
        h = histogram(hitPercentBaseline,nBins);
        h.FaceColor = bluePurpleRed(end,:); h.EdgeColor = bluePurpleRed(end,:);
        hold on
        xline(hitPercentCue,'-r',{'Hit rate',['p=',num2str(pval)]},...
            'LineWidth',3,...
            'LabelOrientation','horizontal');
        xlim([0,1]); xlabel('Hit rate'); ylabel('Count'); box off
        title("Baseline licks (stim only trials)");
    
        % For tone only
        % Calculate cue window success rate
        hitPercentCue = length(find(trials.isTone == 1 & trials.isStim == 0 & trials.nAnticipatoryLicks >= options.minLicks))/size(toneIdx,1);
        pval = sum(hitPercentBaseline>=hitPercentCue)/nboot;
        % Plot distribution
        nexttile; 
        h = histogram(hitPercentBaseline,nBins);
        h.FaceColor = bluePurpleRed(350,:); h.EdgeColor = bluePurpleRed(350,:);
        hold on
        xline(hitPercentCue,'-r',{'Hit rate',['p=',num2str(pval)]},...
            'LineWidth',3,...
            'LabelOrientation','horizontal');
        xlim([0,1]); xlabel('Hit rate'); ylabel('Count'); box off
        title("Baseline licks (tone only trials)");
    end

    % Plot trend
    nexttile;
    plot(trials{1:end-1,"TrialNumber"},trials{1:end-1,"nBaselineLicks"},Color=[0.75,0.75,0.75],LineWidth=2); hold on
    scatter(pairTrials(:,1),trials{pairTrials(:,1),"nBaselineLicks"},100,bluePurpleRed(150,:),'filled'); hold on
    scatter(stimTrials(:,1),trials{stimTrials(:,1),"nBaselineLicks"},100,bluePurpleRed(end,:),'filled'); hold on
    scatter(toneTrials(:,1),trials{toneTrials(:,1),"nBaselineLicks"},100,bluePurpleRed(350,:),'filled'); hold on
    xlabel("Trials"); ylabel("Baseline licks"); box off
    title("Baseline licks (all trials)");


    stageCutoff = linspace(1,size(trials,1),4);
    stageCutoff = stageCutoff(2:3);
    
    % 2. Distribution of first lick reaction time
    rt_pair = trials{pairTrials(:,1),["TrialNumber","ReactionTime"]};
    rt_stim = trials{stimTrials(:,1),["TrialNumber","ReactionTime"]};
    rt_tone = trials{toneTrials(:,1),["TrialNumber","ReactionTime"]};
    rt_pair_early = rt_pair(rt_pair(:,1)<stageCutoff(1),2);
    rt_stim_early = rt_stim(rt_stim(:,1)<stageCutoff(1),2);
    rt_tone_early = rt_tone(rt_tone(:,1)<stageCutoff(1),2);
    rt_pair_mid = rt_pair(rt_pair(:,1)>=stageCutoff(1) & rt_pair(:,1)<stageCutoff(2),2);
    rt_stim_mid = rt_stim(rt_stim(:,1)>=stageCutoff(1) & rt_stim(:,1)<stageCutoff(2),2);
    rt_tone_mid = rt_tone(rt_tone(:,1)>=stageCutoff(1) & rt_tone(:,1)<stageCutoff(2),2);
    rt_pair_late = rt_pair(rt_pair(:,1)<=stageCutoff(2),2);
    rt_stim_late = rt_stim(rt_stim(:,1)<=stageCutoff(2),2);
    rt_tone_late = rt_tone(rt_tone(:,1)<=stageCutoff(2),2);
    

    % For early stage of session
    % Calculate bootstrap distribution
    stages = {rt_pair_early;rt_stim_early;rt_tone_early};
    colors = [bluePurpleRed(150,:);bluePurpleRed(end,:);bluePurpleRed(350,:)];
    nexttile;
    for s = 1:size(stages,1)
        if isempty(stages{s}) || numel(stages{s})==1; continue; end
        [~,bootsam] = bootstrp(nboot,[],stages{s});
        bs = stages{s}(bootsam) / params.sync.behaviorFs;

        h = histogram(bs,nBins); 
        h.FaceColor = colors(s,:); h.EdgeColor = colors(s,:); hold on
    end
    xlabel('First lick reaction time (s)'); ylabel('Count'); box off
    title("First lick reaction time (early stage)");

    % For middle stage of session
    % Calculate bootstrap distribution
    stages = {rt_pair_mid;rt_stim_mid;rt_tone_mid};
    colors = [bluePurpleRed(150,:);bluePurpleRed(end,:);bluePurpleRed(350,:)];
    nexttile;
    for s = 1:size(stages,1)
        if isempty(stages{s}) || numel(stages{s})==1; continue; end
        [~,bootsam] = bootstrp(nboot,[],stages{s});
        bs = stages{s}(bootsam) / params.sync.behaviorFs;

        h = histogram(bs,nBins); 
        h.FaceColor = colors(s,:); h.EdgeColor = colors(s,:); hold on
    end
    xlabel('First lick reaction time (s)'); ylabel('Count'); box off
    title("First lick reaction time (middle stage)");

    % For late stage of session
    % Calculate bootstrap distribution
    stages = {rt_pair_late;rt_stim_late;rt_tone_late};
    colors = [bluePurpleRed(150,:);bluePurpleRed(end,:);bluePurpleRed(350,:)];
    nexttile;
    for s = 1:size(stages,1)
        if isempty(stages{s}) || numel(stages{s})==1; continue; end
        [~,bootsam] = bootstrp(nboot,[],stages{s});
        bs = stages{s}(bootsam) / params.sync.behaviorFs;

        h = histogram(bs,nBins); 
        h.FaceColor = colors(s,:); h.EdgeColor = colors(s,:); hold on
    end
    xlabel('First lick reaction time (s)'); ylabel('Count'); box off
    title("First lick reaction time (late stage)");

    % Plot trend
    nexttile;
    plot(trials{1:end-1,"TrialNumber"},trials{1:end-1,"ReactionTime"}/params.sync.behaviorFs,Color=[0.75,0.75,0.75],LineWidth=2); hold on
    scatter(pairTrials(:,1),trials{pairTrials(:,1),"ReactionTime"}/params.sync.behaviorFs,100,bluePurpleRed(150,:),'filled'); hold on
    scatter(stimTrials(:,1),trials{stimTrials(:,1),"ReactionTime"}/params.sync.behaviorFs,100,bluePurpleRed(end,:),'filled'); hold on
    scatter(toneTrials(:,1),trials{toneTrials(:,1),"ReactionTime"}/params.sync.behaviorFs,100,bluePurpleRed(350,:),'filled'); hold on
    xlabel("Trials"); ylabel("Reaction time (s)"); box off
    title("First lick reaction time (all trials)");


    % 3. Distribution of anticipatory licks
    al_pair = trials{pairTrials(:,1),["TrialNumber","nAnticipatoryLicks"]};
    al_stim = trials{stimTrials(:,1),["TrialNumber","nAnticipatoryLicks"]};
    al_tone = trials{toneTrials(:,1),["TrialNumber","nAnticipatoryLicks"]};
    al_pair_early = al_pair(al_pair(:,1)<stageCutoff(1),2);
    al_stim_early = al_stim(al_stim(:,1)<stageCutoff(1),2);
    al_tone_early = al_tone(al_tone(:,1)<stageCutoff(1),2);
    al_pair_mid = al_pair(al_pair(:,1)>=stageCutoff(1) & al_pair(:,1)<stageCutoff(2),2);
    al_stim_mid = al_stim(al_stim(:,1)>=stageCutoff(1) & al_stim(:,1)<stageCutoff(2),2);
    al_tone_mid = al_tone(al_tone(:,1)>=stageCutoff(1) & al_tone(:,1)<stageCutoff(2),2);
    al_pair_late = al_pair(al_pair(:,1)<=stageCutoff(2),2);
    al_stim_late = al_stim(al_stim(:,1)<=stageCutoff(2),2);
    al_tone_late = al_tone(al_tone(:,1)<=stageCutoff(2),2);

    % For early stage of session
    % Calculate bootstrap distribution
    stages = {al_pair_early;al_stim_early;al_tone_early};
    colors = [bluePurpleRed(150,:);bluePurpleRed(end,:);bluePurpleRed(350,:)];
    nexttile;
    for s = 1:size(stages,1)
        if isempty(stages{s}) || numel(stages{s})==1; continue; end
        [~,bootsam] = bootstrp(nboot,[],stages{s});
        bs = stages{s}(bootsam);

        h = histogram(bs,nBins); 
        h.FaceColor = colors(s,:); h.EdgeColor = colors(s,:); hold on
    end
    % xline(options.minLicks,'-.','Big reward cutoff','LineWidth',3,'LabelOrientation','horizontal');
    xlabel('Anticipatory licks'); ylabel('Count'); box off
    title("Anticipatory licks (early stage)");

    % For middle stage of session
    % Calculate bootstrap distribution
    stages = {al_pair_mid;al_stim_mid;al_tone_mid};
    colors = [bluePurpleRed(150,:);bluePurpleRed(end,:);bluePurpleRed(350,:)];
    nexttile;
    for s = 1:size(stages,1)
        if isempty(stages{s}) || numel(stages{s})==1; continue; end
        [~,bootsam] = bootstrp(nboot,[],stages{s});
        bs = stages{s}(bootsam);

        h = histogram(bs,nBins); 
        h.FaceColor = colors(s,:); h.EdgeColor = colors(s,:); hold on
    end
    % xline(options.minLicks,'-.','Big reward cutoff','LineWidth',3,'LabelOrientation','horizontal');
    xlabel('Anticipatory licks'); ylabel('Count'); box off
    title("Anticipatory licks (middle stage)");

    % For late stage of session
    % Calculate bootstrap distribution
    stages = {al_pair_late;al_stim_late;al_tone_late};
    colors = [bluePurpleRed(150,:);bluePurpleRed(end,:);bluePurpleRed(350,:)];
    nexttile;
    for s = 1:size(stages,1)
        if isempty(stages{s}) || numel(stages{s})==1; continue; end
        [~,bootsam] = bootstrp(nboot,[],stages{s});
        bs = stages{s}(bootsam);

        h = histogram(bs,nBins); 
        h.FaceColor = colors(s,:); h.EdgeColor = colors(s,:); hold on
    end
    % xline(options.minLicks,'-.','Big reward cutoff','LineWidth',3,'LabelOrientation','horizontal');
    xlabel('Anticipatory licks'); ylabel('Count'); box off
    title("Anticipatory licks (late stage)");

    % Plot trend
    nexttile;
    yline(options.minLicks,'-.','Big reward cutoff','LineWidth',3,'LabelOrientation','horizontal'); hold on
    plot(trials{1:end-1,"TrialNumber"},trials{1:end-1,"nAnticipatoryLicks"},Color=[0.75,0.75,0.75],LineWidth=2); hold on
    scatter(pairTrials(:,1),trials{pairTrials(:,1),"nAnticipatoryLicks"},100,bluePurpleRed(150,:),'filled'); hold on
    scatter(stimTrials(:,1),trials{stimTrials(:,1),"nAnticipatoryLicks"},100,bluePurpleRed(end,:),'filled'); hold on
    scatter(toneTrials(:,1),trials{toneTrials(:,1),"nAnticipatoryLicks"},100,bluePurpleRed(350,:),'filled'); hold on
    xlabel("Trials"); ylabel("Anticipatory licks"); box off
    title("Anticipatory licks (all trials)");


    % 4. Distribution of reward rection time
    ort_pair = trials{pairTrials(:,1),["TrialNumber","OutcomeReactionTime"]};
    ort_stim = trials{stimTrials(:,1),["TrialNumber","OutcomeReactionTime"]};
    ort_tone = trials{toneTrials(:,1),["TrialNumber","OutcomeReactionTime"]};
    ort_pair_early = ort_pair(ort_pair(:,1)<stageCutoff(1),2);
    ort_stim_early = ort_stim(ort_stim(:,1)<stageCutoff(1),2);
    ort_tone_early = ort_tone(ort_tone(:,1)<stageCutoff(1),2);
    ort_pair_mid = ort_pair(ort_pair(:,1)>=stageCutoff(1) & ort_pair(:,1)<stageCutoff(2),2);
    ort_stim_mid = ort_stim(ort_stim(:,1)>=stageCutoff(1) & ort_stim(:,1)<stageCutoff(2),2);
    ort_tone_mid = ort_tone(ort_tone(:,1)>=stageCutoff(1) & ort_tone(:,1)<stageCutoff(2),2);
    ort_pair_late = ort_pair(ort_pair(:,1)<=stageCutoff(2),2);
    ort_stim_late = ort_stim(ort_stim(:,1)<=stageCutoff(2),2);
    ort_tone_late = ort_tone(ort_tone(:,1)<=stageCutoff(2),2);

    % For early stage of session
    % Calculate bootstrap distribution
    stages = {ort_pair_early;ort_stim_early;ort_tone_early};
    colors = [bluePurpleRed(150,:);bluePurpleRed(end,:);bluePurpleRed(350,:)];
    nexttile;
    for s = 1:size(stages,1)
        if isempty(stages{s}) || numel(stages{s})==1; continue; end
        [~,bootsam] = bootstrp(nboot,[],stages{s});
        bs = stages{s}(bootsam) / params.sync.behaviorFs;

        h = histogram(bs,nBins); 
        h.FaceColor = colors(s,:); h.EdgeColor = colors(s,:); hold on
    end
    xlabel('Outcome reaction time (s)'); ylabel('Count'); box off
    title("Outcome reaction time (early stage)");

    % For middle stage of session
    % Calculate bootstrap distribution
    stages = {ort_pair_mid;ort_stim_mid;ort_tone_mid};
    colors = [bluePurpleRed(150,:);bluePurpleRed(end,:);bluePurpleRed(350,:)];
    nexttile;
    for s = 1:size(stages,1)
        if isempty(stages{s}) || numel(stages{s})==1; continue; end
        [~,bootsam] = bootstrp(nboot,[],stages{s});
        bs = stages{s}(bootsam) / params.sync.behaviorFs;

        h = histogram(bs,nBins); 
        h.FaceColor = colors(s,:); h.EdgeColor = colors(s,:); hold on
    end
    xlabel('Outcome reaction time (s)'); ylabel('Count'); box off
    title("Outcome reaction time (middle stage)");

    % For late stage of session
    % Calculate bootstrap distribution
    stages = {ort_pair_late;ort_stim_late;ort_tone_late};
    colors = [bluePurpleRed(150,:);bluePurpleRed(end,:);bluePurpleRed(350,:)];
    nexttile;
    for s = 1:size(stages,1)
        if isempty(stages{s}) || numel(stages{s})==1; continue; end
        [~,bootsam] = bootstrp(nboot,[],stages{s});
        bs = stages{s}(bootsam) / params.sync.behaviorFs;

        h = histogram(bs,nBins); 
        h.FaceColor = colors(s,:); h.EdgeColor = colors(s,:); hold on
    end
    xlabel('Outcome reaction time (s)'); ylabel('Count'); box off
    title("Outcome reaction time (late stage)");

    % Plot trend
    nexttile;
    plot(trials{1:end-1,"TrialNumber"},trials{1:end-1,"OutcomeReactionTime"}/params.sync.behaviorFs,Color=[0.75,0.75,0.75],LineWidth=2); hold on
    scatter(pairTrials(:,1),trials{pairTrials(:,1),"OutcomeReactionTime"}/params.sync.behaviorFs,100,bluePurpleRed(150,:),'filled'); hold on
    scatter(stimTrials(:,1),trials{stimTrials(:,1),"OutcomeReactionTime"}/params.sync.behaviorFs,100,bluePurpleRed(end,:),'filled'); hold on
    scatter(toneTrials(:,1),trials{toneTrials(:,1),"OutcomeReactionTime"}/params.sync.behaviorFs,100,bluePurpleRed(350,:),'filled'); hold on
    xlabel("Trials"); ylabel("Outcome reaction time (s)"); box off
    title("Outcome reaction time (all trials)");
    
    saveas(gcf,strcat(sessionpath,filesep,'Summary_distributions.png'));

end
disp('Finished: all plots and struct are plotted and saved!');
return

end