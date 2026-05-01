function analyzeSessions_clamp(sessionpath,options)

arguments
    sessionpath string
    
    options.task
    options.outputName string %name of the output file in subfolder of the session
    
    options.laserSource string = 'clamp'
    options.intervalThreshold double = 100 % max interval of both clamps are off to be consider as clamp ON, in ms

    options.pavlovian logical = false
    options.reactionTime double = 1
    options.minLicks double = 2 % min licks to get reward
    options.combineOmission logical = false % combine omission trials

    options.analyzeTraces logical = true

    options.redo logical = true % Recalculate trial table and all preprocessing
    options.round logical = false % Round reward/airpuff/tone to get duration data
    options.performing logical = false % Only plot traces where the animal performs
    
    options.plotPhotometry logical = true % Plot photometry summary plot
    options.plotBehavior logical = true % Plot lick raster summary plot

    options.lick_binSize double = 0.1
    options.shortTimeRange double = [-1,5]
    options.longTimeRange double = [-5,10]
end

%% Notes
% Shun_analyzeBehavior_clamp
% Shun Li, 4/15/2026

%% Load data

[~,~,~,~,~,~,bluePurpleRed] = loadColors;
clampColor = [.232 .76 .58];
unclampColor = [217, 237, 223]./255;
             
% 1. Select session via uigetdir
dirsplit = strsplit(sessionpath,filesep); 
if ~isfield(options,'outputName'); options.outputName = dirsplit{end}; end
if ispc; projectPath = strcat('\\',fullfile(dirsplit{2:end-1}));
elseif isunix; projectPath = strcat('/',fullfile(dirsplit{2:end-1}));
end

% Get animal name and session date
if ~contains(options.outputName,{'-','_'})
    sessionName = dirsplit{end-1};
    dirsplit = strsplit(sessionName,{'-','_'});
else
    sessionName = options.outputName;
    dirsplit = strsplit(options.outputName,{'-','_'}); 
end
date = dirsplit{1}; animal = dirsplit{2}; sessionTask = dirsplit{3};
clear dirsplit

disp(strcat('**********',options.outputName,'**********'));
load(strcat(sessionpath,filesep,'timeseries_',options.outputName,'.mat'));
load(strcat(sessionpath,filesep,'data_',options.outputName,'.mat'));
load(strcat(sessionpath,filesep,'behavior_',options.outputName,'.mat'));
load(strcat(sessionpath,filesep,'sync_',options.outputName,'.mat'));

if ~isfield(params.session,'name'); params.session.name = options.outputName; end
if ~isfield(params.session,'date'); params.session.date = date; end
if ~isfield(params.session,'animal'); params.session.animal = animal; end
if ~isfield(params.session,'projectPath'); params.session.projectPath = projectPath; end

% Create analysis.mat
if ~isempty(dir(fullfile(sessionpath,"analysis_*.mat")))
    load(strcat(sessionpath,filesep,'analysis_',options.outputName,'.mat'));
else
    save(strcat(sessionpath,filesep,'analysis_',options.outputName),'sessionName','-v7.3');
    disp('Finished: analysis_.mat not found, created a new one');
end
disp(['Finished: Session ',options.outputName,' loaded']);

%% Load behaivor params

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
    if contains(sessionTask,["Random","H"],IgnoreCase=false)
        options.task = 'random';
        params.session.task = options.task;
    elseif contains(sessionTask,["P","PP","Punish"],IgnoreCase=false)
        options.task = 'punish pairing';
        params.session.task = options.task;
        if contains(sessionTask,"RP",IgnoreCase=false)
            options.task = 'reward pairing';
            params.session.task = options.task;
        end
    elseif contains(sessionTask,["R","RP","Reward"],IgnoreCase=false)
        options.task = 'reward pairing';
        params.session.task = options.task;
    else
        warning('Unrecognize session name pattern, set to random');
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
    save(strcat(sessionpath,filesep,'sync_',options.outputName),'params','-append');
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
if params.session.withCamera
    eyeAreaIdx = find(cellfun(@(x) strcmpi(x,'eyeArea'), {timeSeries.name}));
    if ~isempty(eyeAreaIdx); eyeArea_detrend = timeSeries(eyeAreaIdx).data; end
    pupilAreaIdx = find(cellfun(@(x) strcmpi(x,'pupilArea'), {timeSeries.name}));
    if ~isempty(pupilAreaIdx); pupilArea_detrend = timeSeries(pupilAreaIdx).data; end
end

save(strcat(sessionpath,filesep,'sync_',options.outputName),'params','-append');

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
        save(strcat(sessionpath,filesep,'timeseries_',options.outputName),'rightSolenoid_rounded','airpuff_rounded','-append');
        disp('Finished: rounding cue/outcome data');
    end
else
    rightSolenoid_rounded = rightSolenoid;
    airpuff_rounded = airpuff;
end


% Find start of lick bout
rightLickON = find(rightLick);
if ~exist('lickBout','var') || options.redo
    % Get lick bout start time (ILI < 0.5s)
    lickBout = getLickBout(rightLickON);
    save(strcat(sessionpath,filesep,'timeseries_',options.outputName),"lickBout",'-append');
end

% Define trial start
waterIdx = find(rightSolenoid_rounded);  
airpuffIdx = find(airpuff_rounded);

if ~exist('trials','var') || options.redo
    if strcmp(options.task,'random')
        [allTrials,~] = getTrials(find(leftTone),...
                                 waterIdx,airpuffIdx); % add rightLickON if only baseline licks
    elseif contains(options.task,'punish')
        [allTrials,~] = getTrials(find(leftTone),waterIdx);
    else
        [allTrials,~] = getTrials(find(leftTone));
    end
    save(strcat(sessionpath,filesep,'timeseries_',options.outputName),"allTrials",'-append');
end

% Find water lick (first lick in response to water)
waterLickIdx = nan(size(waterIdx));
for e = 1:length(waterIdx)
    nextLick = rightLickON(find(rightLickON>=waterIdx(e),1));
    if ~isempty(nextLick); waterLickIdx(e) = nextLick; end
end
waterLickIdx = rmmissing(waterLickIdx);


%% Find clamp ON time

if ~exist('clampON','var') || options.redo
    gap_ms = options.intervalThreshold;   % fill OFF gaps shorter than this
    noiseThresholdPct = 20;
    smooth_ms = 3;               % small smoothing window for noise resistance
    onPct = 99;
    offPct = 5;

    blueOn = analogToLogical(blueClamp, params.sync.behaviorFs, ...
                                onPct=onPct, offPct=offPct,...
                                noiseThresholdPct=noiseThresholdPct, smooth_ms=smooth_ms);
    redOn  = analogToLogical(redClamp,  params.sync.behaviorFs, ...
                                onPct=onPct, offPct=offPct,...
                                noiseThresholdPct=noiseThresholdPct, smooth_ms=smooth_ms);
    clampON = blueOn | redOn; % Combine channels
    
    % Fill short OFF gaps
    maxGapSamples = round(gap_ms/1000 * params.sync.behaviorFs);
    offRuns = ~clampON;
    d = diff([false; offRuns(:); false]);
    runStarts = find(d == 1);
    runEnds   = find(d == -1) - 1;
    for k = 1:numel(runStarts)
        runLen = runEnds(k) - runStarts(k) + 1;
        % fill only internal OFF gaps shorter than threshold
        if runLen <= maxGapSamples && runStarts(k) > 1 && runEnds(k) < numel(clampON)
            clampON(runStarts(k):runEnds(k)) = true;
        end
    end
    clampON = reshape(clampON, size(blueClamp));
    % Find rising edge of clampON
    clampOnset = find(diff(clampON > 0.5) > 0) + 1;

    %% Code to check 
    % sampleStart_sec = 96;
    % sampleWindow_sec = 4;
    % 
    % sampleStart = sampleStart_sec*params.sync.behaviorFs;
    % sampleWindow = sampleStart : sampleStart+sampleWindow_sec*params.sync.behaviorFs;
    % sampleTime = sampleWindow ./ params.sync.behaviorFs;
    % close all; initializeFig(0.8,0.4); tiledlayout(2,1);
    % 
    % nexttile;
    % plot(sampleTime, blueClamp(sampleWindow)); hold on;
    % plot(sampleTime, redClamp(sampleWindow)); hold on;
    % nexttile;
    % % plot(sampleTime, blueOn(sampleWindow)); hold on;
    % % plot(sampleTime, redOn(sampleWindow)); hold on;
    % plot(sampleTime, clampON(sampleWindow));
    % xlabel('time (s)');
    
    %%
    save(strcat(sessionpath,filesep,'data_',options.outputName),'clampON','clampOnset','-append');
end
disp('Finished: preprocess outcome and opto data');

%% Generate trial and event table

if (~exist('trials','var') || options.redo)

    disp('Ongoing: making trial table');
    events{1} = allTrials;      events{2} = airpuffIdx;
    events{3} = waterIdx;       events{4} = rightLickON;
    events{5} = find(leftTone); events{6} = clampON;

    trials = getTrialTable(options.task,events,rightSolenoid_rounded,airpuff_rounded,...
                pavlovian=options.pavlovian,reactionTime=options.reactionTime,...
                withClamp=true);

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
    % parquetwrite(strcat(sessionpath,filesep,'trialTable.parquet'),trialTable);
    % % eventTable
    % eventTable = getEventTable(events,params);
    % parquetwrite(strcat(sessionpath,filesep,'eventTable.parquet'),eventTable);
    % % blockTable
    % firstTrial = 1; lastTrial = size(trials,1);
    % blockTable = table(firstTrial,lastTrial);
    % parquetwrite(strcat(sessionpath,filesep,'blockTable.parquet'),blockTable);
    

    % Save to behavior_.mat
    % save(strcat(sessionpath,filesep,'behavior_',options.outputName),'trials',...
    %     'eventTable','trialTable','blockTable','-append');
    save(strcat(sessionpath,filesep,'behavior_',options.outputName),'trials','trialTable','-append');
    disp('Finished: trial table saved');
end


%% Task specific params

if strcmp(options.task,'random')
    
    % Select event idxs
    toneIdx = find(leftTone);
    clampIdx = clampOnset;

    % Select baseline idx (align to each baseline lick)
    % baselineLicks = cell2mat(trials{~cellfun(@isempty,trials{1:end-1,'BaselineLicks'}),'BaselineLicks'});
    % if isempty(baselineLicks); baselineIdx = [];
    % else; baselineIdx = baselineLicks(:,1); end
    randomMinSample = 15*params.sync.behaviorFs;
    randomMaxSample = length(params.sync.timeNI) - (15*params.sync.behaviorFs);
    baselineIdx = randi([randomMinSample,randomMaxSample],100,1); % changed to rightLickON if only baseline licks
    
    % Create task legend
    stageTime = [-2,0;0,2];
    analysisEvents = {waterIdx,waterLickIdx,airpuffIdx,toneIdx, clampIdx, baselineIdx};
    eventTrialNum = {findTrials(waterIdx,trials),findTrials(waterLickIdx,trials),...
                        findTrials(airpuffIdx,trials),findTrials(toneIdx,trials),...
                        findTrials(clampIdx,trials),findTrials(baselineIdx,trials)};
    analysisLabels = {'Water','Rewarded licks','Airpuff','Tone','Clamp','Baseline'};
    taskLegend = getLegend(analysisEvents,analysisLabels);

    stageColors = {[.75 .75 .75],bluePurpleRed(1,:)};
    stageLegend = {'Baseline','US'};
    
    eventTrialNum = eventTrialNum(~cellfun('isempty',analysisEvents));
    analysisLabels = analysisLabels(~cellfun('isempty',analysisEvents));
    analysisEvents = analysisEvents(~cellfun('isempty',analysisEvents));
    
    for i = 1:length(analysisEvents)
        disp(['Total ',analysisLabels{i},': ',num2str(length(analysisEvents{i}))]);
    end

    save(strcat(sessionpath,filesep,'behavior_',options.outputName),'allTrials',...
        'waterIdx','waterLickIdx','airpuffIdx','clampIdx',...
        '-append');

else
    % Select baseline idx (align to each baseline lick)
    % baselineLicks = cell2mat(trials{~cellfun(@isempty,trials{1:end-1,'BaselineLicks'}),'BaselineLicks'});
    % if isempty(baselineLicks); baselineIdx = round((trials{2:end,"CueTime"} - trials{2:end,"ENL"}) - 5*params.sync.behaviorFs);
    % else; baselineIdx = baselineLicks(:,1); end
    randomMinSample = floor(10*params.sync.behaviorFs);
    randomMaxSample = length(params.sync.timeNI) - randomMinSample;
    baselineIdx = randi([randomMinSample,randomMaxSample],100,1);

    if options.combineOmission
        if options.performing
            clampTrials = trials{trials.isTone == 1 & trials.isClamp == 1 & trials.performing == 1,["TrialNumber","CueTime","OutcomeTime","ENL"]};
            unclampTrials = trials{trials.isTone == 1 & trials.isClamp == 0 & trials.performing == 1,["TrialNumber","CueTime","OutcomeTime","ENL"]};
            clampIdx = clampTrials(:,2);
            unclampIdx = unclampTrials(:,2);
            clampOmissionIdx = []; unclampOmissionIdx = [];
        else
            clampTrials = trials{trials.isTone == 1 & trials.isClamp == 1,["TrialNumber","CueTime","OutcomeTime","ENL"]};
            unclampTrials = trials{trials.isTone == 1 & trials.isClamp == 0,["TrialNumber","CueTime","OutcomeTime","ENL"]};
            clampIdx = clampTrials(:,2);
            unclampIdx = unclampTrials(:,2);
            clampOmissionIdx = []; unclampOmissionIdx = [];
        end
    else
        if contains(options.task,'reward')
            if options.performing
                clampTrials = trials{trials.isTone == 1 & trials.isClamp == 1 ...
                    & trials.isReward == 1 & trials.performing == 1,["TrialNumber","CueTime","OutcomeTime","ENL"]};
                clampOmissionTrials = trials{trials.isTone == 0 & trials.isClamp == 1 ...
                    & trials.isReward == 0 & trials.performing == 1,["TrialNumber","CueTime","OutcomeTime","ENL"]};
                unclampTrials = trials{trials.isTone == 1 & trials.isClamp == 0 ...
                    & trials.isReward == 1 & trials.performing == 1,["TrialNumber","CueTime","OutcomeTime","ENL"]};
                unclampOmissionTrials = trials{trials.isTone == 1 & trials.isClamp == 0 ...
                    & trials.isReward == 0 & trials.performing == 1,["TrialNumber","CueTime","OutcomeTime","ENL"]};
                clampIdx = clampTrials(:,2);
                clampOmissionIdx = clampOmissionTrials(:,2);
                unclampIdx = unclampTrials(:,2);
                unclampOmissionIdx = unclampOmissionTrials(:,2);
            else
                clampTrials = trials{trials.isTone == 1 & trials.isClamp == 1 & trials.isReward == 1,["TrialNumber","CueTime","OutcomeTime","ENL"]};
                clampOmissionTrials = trials{trials.isTone == 1 & trials.isClamp == 1 & trials.isReward == 0,["TrialNumber","CueTime","OutcomeTime","ENL"]};
                unclampTrials = trials{trials.isTone == 1 & trials.isClamp == 0 & trials.isReward == 1,["TrialNumber","CueTime","OutcomeTime","ENL"]};
                unclampOmissionTrials = trials{trials.isTone == 1 & trials.isClamp == 0 & trials.isReward == 0,["TrialNumber","CueTime","OutcomeTime","ENL"]};
                clampIdx = clampTrials(:,2);
                clampOmissionIdx = clampOmissionTrials(:,2);
                unclampIdx = unclampTrials(:,2);
                unclampOmissionIdx = unclampOmissionTrials(:,2);
            end
        elseif contains(options.task,'punish')
            if options.performing
                clampTrials = trials{trials.isTone == 1 & trials.isClamp == 1 ...
                    & trials.isPunishment == 1 & trials.performing == 1,["TrialNumber","CueTime","OutcomeTime","ENL"]};
                clampOmissionTrials = trials{trials.isTone == 1 & trials.isClamp == 1 ...
                    & trials.isPunishment == 0 & trials.performing == 1,["TrialNumber","CueTime","OutcomeTime","ENL"]};
                unclampTrials = trials{trials.isTone == 1 & trials.isClamp == 0 ...
                    & trials.isPunishment == 1 & trials.performing == 1,["TrialNumber","CueTime","OutcomeTime","ENL"]};
                unclampOmissionTrials = trials{trials.isTone == 1 & trials.isClamp == 0 ...
                    & trials.isPunishment == 0 & trials.performing == 1,["TrialNumber","CueTime","OutcomeTime","ENL"]};
                clampIdx = clampTrials(:,2);
                clampOmissionIdx = clampOmissionTrials(:,2);
                unclampIdx = unclampTrials(:,2);
                unclampOmissionIdx = unclampOmissionTrials(:,2);
            else
                clampTrials = trials{trials.isTone == 1 & trials.isClamp == 1 & trials.isPunishment == 1,["TrialNumber","CueTime","OutcomeTime","ENL"]};
                clampOmissionTrials = trials{trials.isTone == 1 & trials.isClamp == 1 & trials.isPunishment == 0,["TrialNumber","CueTime","OutcomeTime","ENL"]};
                unclampTrials = trials{trials.isTone == 1 & trials.isClamp == 0 & trials.isPunishment == 1,["TrialNumber","CueTime","OutcomeTime","ENL"]};
                unclampOmissionTrials = trials{trials.isTone == 1 & trials.isClamp == 0 & trials.isPunishment == 0,["TrialNumber","CueTime","OutcomeTime","ENL"]};
                clampIdx = clampTrials(:,2);
                clampOmissionIdx = clampOmissionTrials(:,2);
                unclampIdx = unclampTrials(:,2);
                unclampOmissionIdx = unclampOmissionTrials(:,2);
            end
        end
    end

    stageTime = [-2,0;0,0.5;0.5,5];
    analysisEvents = {waterIdx,waterLickIdx,clampIdx,unclampIdx,airpuffIdx,baselineIdx};
    eventTrialNum = {findTrials(waterIdx,trials),findTrials(waterLickIdx,trials),...
                    clampTrials(:,1),unclampTrials(:,1),...
                    findTrials(airpuffIdx,trials),findTrials(baselineIdx,trials)};
    analysisLabels = {'Water','Rewarded licks','Tone (clamp)','Tone (unclamp)','Airpuff','Baseline'};
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
    save(strcat(sessionpath,filesep,'behavior_',options.outputName),'allTrials',...
        'waterIdx','waterLickIdx','airpuffIdx','clampIdx','unclampIdx',...
        'baselineIdx','clampTrials','unclampTrials',...
        '-append');
end

%% Save photometry/lick/eye PSTHs

if options.analyzeTraces
    if ~exist('timeSeries','var')
        timeSeries = struct([]);
    end
    analysis = analyzeTraces(timeSeries,rightLick,analysisEvents,analysisLabels,params,...
                         stageTime=stageTime,...
                         trialNumber=eventTrialNum,trialTable=trials);
end

%% Plot photometry summary plots

if options.plotPhotometry && exist('timeSeries','var')
    %% Find the number of photometry channels
    photometryIdx = find(cellfun(@(x) contains(x,["NI","LJ"],"IgnoreCase",true), {timeSeries.system}));
    photometryName = cellfun(@(x) unique(x,'rows'), {timeSeries(photometryIdx).name},'UniformOutput',false);
    nSignals = length(photometryIdx);
    disp(['Finished: found ', num2str(nSignals),' photometry signals']);

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
        xlabel('\DeltaF/F'); ylabel('Count'); legend({'Normal distribution',timeSeries(path).name});
        
        title(timeSeries(path).name);
        subtitle(strcat("Skewness: ",num2str(skew_lj),", Kurtosis: ",num2str(kur_lj)));
    end
    
    % Save figure
    saveFigures(gcf,'Summary_photometry_distribution',sessionpath,savePDF=true);

    %% Loop through timeSeries struct
    for photometry = 1:nSignals

        % Load signal of interest
        path = photometryIdx(photometry);
        signal = timeSeries(path).data;
        finalFs = timeSeries(path).finalFs;
        system = timeSeries(path).system;
        name = timeSeries(path).name;
        
        %% Plot combined PSTH
        timeRange = options.shortTimeRange;
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
            [~,~] = plotTraces(clampIdx,timeRange,signal,clampColor,params,...
                        signalFs=finalFs,signalSystem=system);
            [~,~] = plotTraces(baselineIdx,timeRange,signal,[.75 .75 .75],params,...
                        signalFs=finalFs,signalSystem=system);
            plotEvent('',0);
            xlabel('Time (s)'); 
            if contains(name,'clamp',IgnoreCase=true)
                ylabel([name,' (%)']);
            else
                ylabel([name,' \DeltaF/F']);
            end
            legend(taskLegend(2:end),'Location','northeast');
            
            % 2.2 Plot lick traces
            nexttile
            plotLicks(waterLickIdx,timeRange,options.lick_binSize,bluePurpleRed(1,:),[],rightLick,params);
            plotLicks(airpuffIdx,timeRange,options.lick_binSize,[0.2, 0.2, 0.2],[],rightLick,params);
            plotLicks(toneIdx,timeRange,options.lick_binSize,bluePurpleRed(350,:),[],rightLick,params);
            plotLicks(clampIdx,timeRange,options.lick_binSize,clampColor,[],rightLick,params);
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
                [~,~] = plotTraces(clampIdx,timeRange,eyeArea_detrend,clampColor,params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(baselineIdx,timeRange,eyeArea_detrend,[.75 .75 .75],params,signalSystem='camera',smooth=15);
                plotEvent('',0);
                xlabel('Time (s)'); ylabel('Eye area (a.u.)');
                legend(taskLegend(2:end),'Location','northeast');
    
                % 2.4 Plot pupil area traces
                nexttile
                [~,~] = plotTraces(waterLickIdx,timeRange,pupilArea_detrend,bluePurpleRed(1,:),params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(airpuffIdx,timeRange,pupilArea_detrend,[0.2, 0.2, 0.2],params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(toneIdx,timeRange,pupilArea_detrend,bluePurpleRed(350,:),params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(clampIdx,timeRange,pupilArea_detrend,clampColor,params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(baselineIdx,timeRange,pupilArea_detrend,[.75 .75 .75],params,signalSystem='camera',smooth=15);
                plotEvent('',0);
                xlabel('Time (s)'); ylabel('Pupil area (a.u.)');
                legend(taskLegend(2:end),'Location','northeast');
            end
    
        else
            % 2.1 Plot photometry traces
            nexttile
            [~,~] = plotTraces(waterLickIdx,timeRange,signal,bluePurpleRed(1,:),params,...
                        signalFs=finalFs,signalSystem=system);
            [~,~] = plotTraces(clampIdx,timeRange,signal,clampColor,params,...
                        signalFs=finalFs,signalSystem=system);
            [~,~] = plotTraces(unclampIdx,timeRange,signal,unclampColor,params,...
                        signalFs=finalFs,signalSystem=system);
            [~,~] = plotTraces(airpuffIdx,timeRange,signal,[0.2, 0.2, 0.2],params,...
                        signalFs=finalFs,signalSystem=system);
            [~,~] = plotTraces(baselineIdx,timeRange,signal,[.75 .75 .75],params,...
                        signalFs=finalFs,signalSystem=system);
            plotEvent('',0);
            xlabel('Time (s)'); 
            if contains(name,'clamp',IgnoreCase=true)
                ylabel([name,' (%)']);
            else
                ylabel([name,' \DeltaF/F']);
            end
            legend(taskLegend(2:end),'Location','northeast');
            
            % 2.2 Plot lick traces
            nexttile
            plotLicks(waterLickIdx,timeRange,options.lick_binSize,bluePurpleRed(1,:),[],rightLick,params);
            plotLicks(clampIdx,timeRange,options.lick_binSize,clampColor,[],rightLick,params);
            plotLicks(unclampIdx,timeRange,options.lick_binSize,unclampColor,[],rightLick,params);
            plotLicks(airpuffIdx,timeRange,options.lick_binSize,[0.2, 0.2, 0.2],[],rightLick,params);
            plotLicks(baselineIdx,timeRange,options.lick_binSize,[.75 .75 .75],[],rightLick,params);
            plotEvent('',0);
            xlabel('Time (s)'); ylabel('Licks/s'); 
            legend(taskLegend(2:end),'Location','best');

            if params.session.withCamera && params.session.withEyeTracking
                % 2.3 Plot eye area traces
                nexttile
                [~,~] = plotTraces(waterLickIdx,timeRange,eyeArea_detrend,bluePurpleRed(1,:),params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(clampIdx,timeRange,eyeArea_detrend,clampColor,params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(unclampIdx,timeRange,eyeArea_detrend,unclampColor,params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(airpuffIdx,timeRange,eyeArea_detrend,[0.2, 0.2, 0.2],params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(baselineIdx,timeRange,eyeArea_detrend,[.75 .75 .75],params,signalSystem='camera',smooth=15);
                plotEvent('',0);
                xlabel('Time (s)'); ylabel('Eye area (a.u.)');
                legend(taskLegend(2:end),'Location','northeast');
    
                % 2.4 Plot pupil area traces
                nexttile
                [~,~] = plotTraces(waterLickIdx,timeRange,pupilArea_detrend,bluePurpleRed(1,:),params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(clampIdx,timeRange,pupilArea_detrend,clampColor,params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(unclampIdx,timeRange,pupilArea_detrend,unclampColor,params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(airpuffIdx,timeRange,pupilArea_detrend,[0.2, 0.2, 0.2],params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(baselineIdx,timeRange,pupilArea_detrend,[.75 .75 .75],params,signalSystem='camera',smooth=15);
                plotEvent('',0);
                xlabel('Time (s)'); ylabel('Pupil area (a.u.)');
                legend(taskLegend(2:end),'Location','northeast');
            end
        end 
        saveFigures(gcf,strcat('Summary_events_',timeSeries(path).name),sessionpath,savePDF=true);

        %% Plot single stimulus PSTH
        if strcmp(options.task,'random')
            eventIdxes = {waterLickIdx,toneIdx,airpuffIdx,clampIdx};
            labels = {'Water','Tone','Airpuff','Clamp'};
            eventDurations = [0,0.25,0.02,5];
            groupSizes = [30,10,30,10];
            longTimeRange = options.longTimeRange;
            shortTimeRange = options.shortTimeRange; 
            
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
                xlabel('Time (s)');
                if contains(name,'clamp',IgnoreCase=true)
                    ylabel([name,' (%)']);
                else
                    ylabel([name,' \DeltaF/F']);
                end
                legend({[label,' (n=',num2str(length(eventIdx)),')']},...
                        'Location','northeast');
                
                nexttile(7,[1 2]);
                legendList = plotGroupTraces(traces,t,bluePurpleRed,groupSize=groupSize);
                plotEvent(label,eventDuration);
                xlabel('Time (s)'); 
                if contains(name,'clamp',IgnoreCase=true)
                    ylabel([name,' (%)']);
                else
                    ylabel([name,' \DeltaF/F']);
                end
                legend(legendList);

                % Plot heatmap
                nexttile(1,[4 2]);
                plotHeatmap(traces,t,centerColorMap=false);
                plotEvent(label,eventDuration);
                xlabel('Time (s)'); ylabel('Trials');


                % Plot long timescale
                nexttile(11,[1 2]);
                [traces,t] = plotTraces(eventIdx,longTimeRange,signal,bluePurpleRed(1,:),params,...
                                signalFs=finalFs,signalSystem=system,...
                                plotIndividual=true);
                plotEvent(label,eventDuration);
                xlabel('Time (s)');
                if contains(name,'clamp',IgnoreCase=true)
                    ylabel([name,' (%)']);
                else
                    ylabel([name,' \DeltaF/F']);
                end
                legend({[label,' (n=',num2str(length(eventIdx)),')']},...
                        'Location','northeast');
                
                nexttile(15,[1 2]);
                legendList = plotGroupTraces(traces,t,bluePurpleRed,groupSize=groupSize);
                plotEvent(label,eventDuration);
                xlabel('Time (s)');
                if contains(name,'clamp',IgnoreCase=true)
                    ylabel([name,' (%)']);
                else
                    ylabel([name,' \DeltaF/F']);
                end
                legend(legendList);
                
                saveFigures(gcf,strcat('Events_',timeSeries(path).name,'_',label),sessionpath,savePDF=true);
            end
        else
            eventIdxes = {clampIdx,unclampIdx,waterLickIdx,airpuffIdx};
            omissionIdxes = {clampOmissionIdx, unclampOmissionIdx,[],[],[]};
            labels = {'Tone (clamp)','Tone (unclamp)','Water','Airpuff'};
            eventDurations = [0.25,0.25,0,0.02];
            groupSizes = [10,10,30,30];
            longTimeRange = options.longTimeRange;
            shortTimeRange = options.shortTimeRange;
            
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
                [~,~] = plotTraces(omissionIdx,shortTimeRange,signal,[0.6, 0.6, 0.6],params,...
                                signalFs=finalFs,signalSystem=system);
                plotEvent(label,eventDuration);
                xlabel('Time (s)');
                if contains(name,'clamp',IgnoreCase=true)
                    ylabel([name,' (%)']);
                else
                    ylabel([name,' \DeltaF/F']);
                end
                legend({[label,' (n=',num2str(length(eventIdx)),')'],...
                        [label,' omission (n=',num2str(length(omissionIdx)),')']},...
                        'Location','northeast');
                
                nexttile(7,[1 2]);
                legendList = plotGroupTraces(traces,t,bluePurpleRed,groupSize=groupSize);
                plotEvent(label,eventDuration);
                xlabel('Time (s)'); 
                if contains(name,'clamp',IgnoreCase=true)
                    ylabel([name,' (%)']);
                else
                    ylabel([name,' \DeltaF/F']);
                end
                legend(legendList);

                % Plot heatmap
                nexttile(1,[4 2]);
                plotHeatmap(traces,t,centerColorMap=false);
                plotEvent(label,eventDuration);
                xlabel('Time (s)'); ylabel('Trials');

                % Plot long time scale
                nexttile(11,[1 2]);
                [traces,t] = plotTraces(eventIdx,longTimeRange,signal,bluePurpleRed(1,:),params,...
                                signalFs=finalFs,signalSystem=system,...
                                plotIndividual=true);
                [~,~] = plotTraces(omissionIdx,longTimeRange,signal,[0.6, 0.6, 0.6],params,...
                                signalFs=finalFs,signalSystem=system);
                plotEvent(label,eventDuration);
                xlabel('Time (s)');
                if contains(name,'clamp',IgnoreCase=true)
                    ylabel([name,' (%)']);
                else
                    ylabel([name,' \DeltaF/F']);
                end
                legend({[label,' (n=',num2str(length(eventIdx)),')'],...
                        [label,' omission (n=',num2str(length(omissionIdx)),')']},...
                        'Location','northeast');
                
                nexttile(15,[1 2]);
                legendList = plotGroupTraces(traces,t,bluePurpleRed,groupSize=groupSize);
                plotEvent(label,eventDuration);
                xlabel('Time (s)'); 
                if contains(name,'clamp',IgnoreCase=true)
                    ylabel([name,' (%)']);
                else
                    ylabel([name,' \DeltaF/F']);
                end
                legend(legendList);

                saveFigures(gcf,strcat('Events_',timeSeries(path).name,'_',label),sessionpath,savePDF=true);
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
                xlabel('Trials'); ylabel([analysis(row).name,' signal average (\DeltaF/F)']);
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
            saveFigures(gcf,strcat('Analysis_',cur_signal,'_subtrial_average'),sessionpath,savePDF=true);
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
                xlabel('Trials'); ylabel([analysis(row).name,' signal peak (\DeltaF/F)']);
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
            saveFigures(gcf,strcat('Analysis_',cur_signal,'_subtrial_peak'),sessionpath,savePDF=true);
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
                xlabel('Trials'); ylabel([analysis(row).name,' signal trough (\DeltaF/F)']);
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
            saveFigures(gcf,strcat('Analysis_',cur_signal,'_subtrial_trough'),sessionpath,savePDF=true);
        end

        %% Plot subtrial area trend for baseline, CS, US
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
                    data = analysis(row).stageArea.data(:,stage);
                    p = analysis(row).stageArea.fit(stage,:);
                    x = 1:length(data);
                    scatter(x,data',100,stageColors{stage},'filled',MarkerFaceAlpha=0.5,HandleVisibility='off'); hold on
                    plot(x,polyval(p,x),Color=stageColors{stage},lineWidth=5);
                end
                title(analysis(row).event);
                xlabel('Trials'); ylabel([analysis(row).name,' signal area (\DeltaF/F)']);
                legend(stageLegend); box off
    
                nexttile(eventPlotted+length(analysisEvents));
                for stage = 1:size(stageTime,1)
                    h = histogram(analysis(row).stageArea.stats.bs(stage,:,1),nBins); hold on
                    h.FaceColor = stageColors{stage}; h.EdgeColor = stageColors{stage};
                    h.FaceAlpha = 0.25; h.EdgeAlpha = 0.25;
                    pval = min(analysis(row).stageArea.stats.pval_slope(stage,:));
                    xline(analysis(row).stageArea.fit(stage,1),'-',...
                        {['Slope (',stageLegend{stage},')'],['p=',num2str(pval)]},...
                        'Color',stageColors{stage},...
                        'LineWidth',3,...
                        'LabelOrientation','horizontal');
                end
                xlabel('Slope distribution (bootstrapped)'); ylabel('Count'); box off
            end
            % Save
            saveFigures(gcf,strcat('Analysis_',cur_signal,'_subtrial_area'),sessionpath,savePDF=true);
        end
    end
end

%% Plot behavior related plots

% Plot lick bout distribution
if options.plotBehavior
    initializeFig(0.5,0.5); tiledlayout('flow');

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

    % Plot lick reaction time
    nexttile;
    histogram(trials{trials.isReward == 1,'OutcomeReactionTime'}/params.sync.behaviorFs,50);
    xlabel('First lick response time to water (s)'); ylabel('Count'); box off

    % Plot ITI-ENL vs trials
    nexttile;
    plot(trials{1:end-1,'TrialNumber'},ENLinSec,'Color',bluePurpleRed(500,:),LineWidth=2); hold on
    plot(trials{1:end-1,'TrialNumber'},ITIextra,'Color',bluePurpleRed(1,:),LineWidth=2);
    xlabel('Trials'); ylabel('Time (s)'); ylim([0,Inf]); box off
    legend({'ENL','ITI-ENL'});

    % Plot lick per trial vs trials
    % nexttile;
    % plot(trials{1:end-1,'TrialNumber'},trials{1:end-1,'nLicks'},'Color',bluePurpleRed(1,:),LineWidth=2);
    % xlabel('Trials'); ylabel('Licks per trial'); box off

    saveFigures(gcf,'Behavior_ITI&LickBout',sessionpath,savePDF=true);
end


if options.plotBehavior && contains(options.task,'pairing')     
    %% Plot session overview for licking
    timeRange = options.longTimeRange; cameraTimeRange = options.shortTimeRange;
    markerSize = 30;
    timeRange = [-3,5];

    % Get event time and number by trial type
    clampTrials(:,3) = clampTrials(:,3)./params.sync.behaviorFs;
    unclampTrials(:,3) = unclampTrials(:,3)./params.sync.behaviorFs;

    % getLicks by trial type
    [clampLickRate,~,clampLicks] = getLicks(timeRange,clampIdx,options.lick_binSize,[],rightLick,...
                                params.sync.behaviorFs,params.sync.timeNI);
    [unclampLickRate,~,unclampLicks] = getLicks(timeRange,unclampIdx,options.lick_binSize,[],rightLick,...
                                params.sync.behaviorFs,params.sync.timeNI);

    % 1. Plot lick raster and trace
    initializeFig(0.67,0.67); tiledlayout(3,2);
    % 1.1 Plot lick raster plot for session
    nexttile([2,1]);
    for i = 1:size(clampLicks,1)
        scatter(clampLicks{i},clampTrials(i,1),markerSize,'filled','MarkerFaceColor',clampColor); hold on
        scatter(clampTrials(i,3),clampTrials(i,1),markerSize+20,bluePurpleRed(1,:),'pentagram','filled'); hold on
    end
    for i = 1:size(unclampLicks,1)
        scatter(unclampLicks{i},unclampTrials(i,1),markerSize,'filled','MarkerFaceColor',unclampColor); hold on
        scatter(unclampTrials(i,3),unclampTrials(i,1),markerSize+20,bluePurpleRed(1,:),'pentagram','filled'); hold on
    end
    xlabel('Time (s)'); xlim([timeRange(1),timeRange(2)]);
    ylabel('Trial'); ylim([0,size(trials,1)]);
    plotEvent("",0.5);
    % 1.2 Plot lick traces across session
    traces = {clampLickRate,unclampLickRate};
    labels = {'Clamp','Unclamp'};
    groupSizes = [20, 20];
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
        plotEvent(label,0.25);
        xlabel('Time (s)'); ylabel('Licks/s');
        legend(legendList);
    end
    saveFigures(gcf,'Behavior_LickOverview',sessionpath,savePDF=true);
    
    %% Plot session overview for eye
    if params.session.withCamera && params.session.withEyeTracking
        initializeFig(0.67,0.67); tiledlayout(2,4); 

        % Plot eye trace heatmap
        nexttile([2 2]);
        [allEyeTraces,t_cam] = plotTraces(trials{1:end-1,"CueTime"},cameraTimeRange,eyeArea_detrend,[0,0,0],params,plot=false,signalSystem='camera'); 
        plotHeatmap(allEyeTraces,t_cam,centerColorMap=false);
        set(gca,'YDir','normal');
        colorbar; box off
        plotEvent('Cue',0.5);
        xlabel('Time (s)'); ylabel('Trials');

        % get camera trace by trial type
        [clampEyeArea,t_cam] = plotTraces(clampIdx,cameraTimeRange,eyeArea_detrend,[0,0,0],params,plot=false,signalSystem='camera');
        [unclampEyeArea,~] = plotTraces(unclampIdx,cameraTimeRange,eyeArea_detrend,[0,0,0],params,plot=false,signalSystem='camera');
        [clampPupilArea,~] = plotTraces(clampIdx,cameraTimeRange,pupilArea_detrend,[0,0,0],params,plot=false,signalSystem='camera');
        [unclampPupilArea,~] = plotTraces(unclampIdx,cameraTimeRange,pupilArea_detrend,[0,0,0],params,plot=false,signalSystem='camera');

        % Plot eye area across session
        traces = {clampEyeArea,unclampEyeArea};
        labels = {'Clamp','Unclamp'};
        groupSizes = [20, 20];
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
            plotEvent(label,0.25);
            xlabel('Time (s)'); ylabel('Eye area (a.u.)');
            legend(legendList);
        end
    
        % Plot pupil area across session
        traces = {clampPupilArea,unclampPupilArea};
        labels = {'Clamp','Unclamp'};
        groupSizes = [20, 20];
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
            plotEvent(label,0.25);
            xlabel('Time (s)'); ylabel('Pupil area (a.u.)');
            legend(legendList);
        end
        saveFigures(gcf,'Behavior_EyeOverview',sessionpath,savePDF=true);
    end

    %% Plot ENL aligned lick trace
    % Bin trials based on ENL length, plot lick trace aligning to ENL start
    % should not smear if animal lick specifically to the cue

    initializeFig(0.67,0.67); tiledlayout(2,3);
    timeRangeENL = [0,10];
    [clampLickRateENL,~,~] = getLicks(timeRangeENL,clampIdx-clampTrials(:,4),options.lick_binSize,[],rightLick,...
                                    params.sync.behaviorFs,params.sync.timeNI);
    [unclampLickRateENL,~,~] = getLicks(timeRangeENL,unclampIdx-unclampTrials(:,4),options.lick_binSize,[],rightLick,...
                                params.sync.behaviorFs,params.sync.timeNI);

    nENLBins = [2,3,4];
    for i = 1:length(nENLBins)

        % 1.1 Sort trials based on ENL duration
        % For clamp
        nexttile;
        trialIdx = clampTrials; traces = clampLickRateENL;
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

        % For unclamp
        nexttile;
        trialIdx = unclampTrials; traces = unclampLickRateENL;
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
    saveFigures(gcf,'Behavior_ENLAlignedLick',sessionpath,savePDF=true);

    %% Plot distributions
    initializeFig(1,1);
    tiledlayout(5,3);
    
    % 1. Distribution of baseline licks
    % Find out baseline period (eg 5s-10s after reward), bootstrap and
    % compare with stim onset
    % (i.e.) how likely I get 2+ licks within a 2s window randomly selected
    % from baseline (a window of 10 secs)?
    timeRangeBaseline = [-10,0]; timeRangeReaction = [0,options.reactionTime];
    nboot = 100000; nBins = 50;
    [~,allLicksBaseline,~] = getLicks(timeRangeBaseline,trials{1:end-1,'CueTime'}-trials{1:end-1,"ENL"},...
        options.lick_binSize,[],rightLick,params.sync.behaviorFs,params.sync.timeNI);
    [~,reactionLicks,~] = getLicks(timeRangeReaction,unclampIdx,...
        options.lick_binSize,[],rightLick,params.sync.behaviorFs,params.sync.timeNI);
    al_unclamp = sum(reactionLicks,2);
    [~,reactionLicks,~] = getLicks(timeRangeReaction,clampIdx,...
        options.lick_binSize,[],rightLick,params.sync.behaviorFs,params.sync.timeNI);
    al_clamp = sum(reactionLicks,2);
    if options.pavlovian
        % Bootstrap baseline licking
        nBaselineTrials = size(allLicksBaseline,1);
        lickCountBaseline = zeros(nBaselineTrials,nboot);
        
        for i = 1:nboot
            windowInBin = round(options.reactionTime / options.lick_binSize);
            baselineSamples = datasample(allLicksBaseline,windowInBin,2,'Replace',true); 
            lickCountBaseline(:,i) = sum(baselineSamples,2);
        end
        % Plot unclamp trials anticipatory distribution
        nexttile; 
        h = histogram(lickCountBaseline,nBins); h.FaceColor = [0.75,0.75,0.75]; h.EdgeColor = [0.75,0.75,0.75]; hold on
        events = al_unclamp(:);   % or al_clamp(:)
        if numel(events) < 2
            warning('Not enough events for bootstrap: n=%d', numel(events));
            h = histogram(events, nBins);   % or just skip plotting
        else
            nboot_event = round(size(trials,1) * nboot / numel(events));
            [~, bootsam] = bootstrp(nboot_event, [], events);
            h = histogram(events(bootsam), nBins);
        end
        h.FaceColor = unclampColor; h.EdgeColor = unclampColor;
        % xline(options.minLicks,'-.','Big reward cutoff','LineWidth',3,'LabelOrientation','horizontal');
        xline(mean(events),'-r',{'Averge anticipatory licks','(unclamp trials)'},'LineWidth',3,'LabelOrientation','horizontal');
        xlabel('Anticipatory licks'); ylabel ('Count'); box off
        title("Baseline licks vs anticipatory licks (unclamp trials)");
        % Plot clamp trials anticipatory distribution
        nexttile; 
        h = histogram(lickCountBaseline,nBins); h.FaceColor = [0.75,0.75,0.75]; h.EdgeColor = [0.75,0.75,0.75]; hold on
        events = al_clamp;
        if numel(events) < 2
            warning('Not enough events for bootstrap: n=%d', numel(events));
            h = histogram(events, nBins);   % or just skip plotting
        else
            nboot_event = round(size(trials,1) * nboot / numel(events));
            [~, bootsam] = bootstrp(nboot_event, [], events);
            h = histogram(events(bootsam), nBins);
        end
        h.FaceColor = clampColor; h.EdgeColor = clampColor;
        % xline(options.minLicks,'-.','Big reward cutoff','LineWidth',3,'LabelOrientation','horizontal');
        xline(mean(events),'-r',{'Averge anticipatory licks','(clamp only trials)'},'LineWidth',3,'LabelOrientation','horizontal');
        xlabel('Anticipatory licks'); ylabel ('Count'); box off
        title("Baseline licks vs anticipatory licks (clamp trials)");
    
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
    
        % For unclamp
        % Calculate cue window success rate
        hitPercentCue = length(find(trials.isTone == 1 & trials.isClamp == 0 & trials.nAnticipatoryLicks >= options.minLicks))/size(unclampIdx,1);
        pval = sum(hitPercentBaseline>=hitPercentCue)/nboot;
        % Plot distribution
        nexttile; 
        h = histogram(hitPercentBaseline,nBins); 
        h.FaceColor = unclampColor; h.EdgeColor = unclampColor;
        hold on
        xline(hitPercentCue,'-r',{'Hit rate',['p=',num2str(pval)]},...
            'LineWidth',3,...
            'LabelOrientation','horizontal');
        xlim([0,1]); xlabel('Hit rate'); ylabel ('Count'); box off
        title("Baseline licks (unclamp trials)");
    
        % For clamp only
        % Calculate cue window success rate
        hitPercentCue = length(find(trials.isTone == 1 & trials.isClamp == 1 & trials.nAnticipatoryLicks >= options.minLicks))/size(clampIdx,1);
        pval = sum(hitPercentBaseline>=hitPercentCue)/nboot;
        % Plot distribution
        nexttile; 
        h = histogram(hitPercentBaseline,nBins);
        h.FaceColor = clampColor; h.EdgeColor = clampColor;
        hold on
        xline(hitPercentCue,'-r',{'Hit rate',['p=',num2str(pval)]},...
            'LineWidth',3,...
            'LabelOrientation','horizontal');
        xlim([0,1]); xlabel('Hit rate'); ylabel('Count'); box off
        title("Baseline licks (clamp only trials)");
    end
    
    % Plot trend
    nexttile;
    plot(trials{1:end-1,"TrialNumber"},trials{1:end-1,"nBaselineLicks"},Color=[0.75,0.75,0.75],LineWidth=2); hold on
    scatter(unclampTrials(:,1),trials{unclampTrials(:,1),"nBaselineLicks"},100,unclampColor,'filled'); hold on
    scatter(clampTrials(:,1),trials{clampTrials(:,1),"nBaselineLicks"},100,clampColor,'filled'); hold on
    xlabel("Trials"); ylabel("Baseline licks"); box off
    title("Baseline licks (all trials)");
    
    
    stageCutoff = linspace(1,size(trials,1),4);
    stageCutoff = stageCutoff(2:3);
    
    % 2. Distribution of first lick reaction time
    rt_unclamp = trials{unclampTrials(:,1),["TrialNumber","ReactionTime"]};
    rt_clamp = trials{clampTrials(:,1),["TrialNumber","ReactionTime"]};
    rt_unclamp_early = rt_unclamp(rt_unclamp(:,1)<stageCutoff(1),2);
    rt_clamp_early = rt_clamp(rt_clamp(:,1)<stageCutoff(1),2);
    rt_unclamp_mid = rt_unclamp(rt_unclamp(:,1)>=stageCutoff(1) & rt_unclamp(:,1)<stageCutoff(2),2);
    rt_clamp_mid = rt_clamp(rt_clamp(:,1)>=stageCutoff(1) & rt_clamp(:,1)<stageCutoff(2),2);
    rt_unclamp_late = rt_unclamp(rt_unclamp(:,1)<=stageCutoff(2),2);
    rt_clamp_late = rt_clamp(rt_clamp(:,1)<=stageCutoff(2),2);
    
    
    % For early stage of session
    % Calculate bootstrap distribution
    stages = {rt_unclamp_early;rt_clamp_early};
    colors = [unclampColor;clampColor];
    nexttile;
    for s = 1:size(stages,1)
        if isempty(stages{s}) || isscalar(stages{s}); continue; end
        [~,bootsam] = bootstrp(nboot,[],stages{s});
        bs = stages{s}(bootsam) / params.sync.behaviorFs;
    
        h = histogram(bs,nBins); 
        h.FaceColor = colors(s,:); h.EdgeColor = colors(s,:); hold on
    end
    xlabel('First lick reaction time (s)'); ylabel('Count'); box off
    title("First lick reaction time (early stage)");
    
    % For middle stage of session
    % Calculate bootstrap distribution
    stages = {rt_unclamp_mid;rt_clamp_mid};
    colors = [unclampColor;clampColor];
    nexttile;
    for s = 1:size(stages,1)
        if isempty(stages{s}) || isscalar(stages{s}); continue; end
        [~,bootsam] = bootstrp(nboot,[],stages{s});
        bs = stages{s}(bootsam) / params.sync.behaviorFs;
    
        h = histogram(bs,nBins); 
        h.FaceColor = colors(s,:); h.EdgeColor = colors(s,:); hold on
    end
    xlabel('First lick reaction time (s)'); ylabel('Count'); box off
    title("First lick reaction time (middle stage)");
    
    % For late stage of session
    % Calculate bootstrap distribution
    stages = {rt_unclamp_late;rt_clamp_late};
    colors = [unclampColor;clampColor];
    nexttile;
    for s = 1:size(stages,1)
        if isempty(stages{s}) || isscalar(stages{s}); continue; end
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
    scatter(unclampTrials(:,1),trials{unclampTrials(:,1),"ReactionTime"}/params.sync.behaviorFs,100,unclampColor,'filled'); hold on
    scatter(clampTrials(:,1),trials{clampTrials(:,1),"ReactionTime"}/params.sync.behaviorFs,100,clampColor,'filled'); hold on
    xlabel("Trials"); ylabel("Reaction time (s)"); box off
    title("First lick reaction time (all trials)");
    
    
    % 3. Distribution of anticipatory licks
    al_unclamp = trials{unclampTrials(:,1),["TrialNumber","nAnticipatoryLicks"]};
    al_clamp = trials{clampTrials(:,1),["TrialNumber","nAnticipatoryLicks"]};
    al_unclamp_early = al_unclamp(al_unclamp(:,1)<stageCutoff(1),2);
    al_clamp_early = al_clamp(al_clamp(:,1)<stageCutoff(1),2);
    al_unclamp_mid = al_unclamp(al_unclamp(:,1)>=stageCutoff(1) & al_unclamp(:,1)<stageCutoff(2),2);
    al_clamp_mid = al_clamp(al_clamp(:,1)>=stageCutoff(1) & al_clamp(:,1)<stageCutoff(2),2);
    al_unclamp_late = al_unclamp(al_unclamp(:,1)<=stageCutoff(2),2);
    al_clamp_late = al_clamp(al_clamp(:,1)<=stageCutoff(2),2);
    
    % For early stage of session
    % Calculate bootstrap distribution
    stages = {al_unclamp_early;al_clamp_early};
    colors = [unclampColor;clampColor];
    nexttile;
    for s = 1:size(stages,1)
        if isempty(stages{s}) || isscalar(stages{s}); continue; end
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
    stages = {al_unclamp_mid;al_clamp_mid};
    colors = [unclampColor;clampColor];
    nexttile;
    for s = 1:size(stages,1)
        if isempty(stages{s}) || isscalar(stages{s}); continue; end
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
    stages = {al_unclamp_late;al_clamp_late};
    colors = [unclampColor;clampColor];
    nexttile;
    for s = 1:size(stages,1)
        if isempty(stages{s}) || isscalar(stages{s}); continue; end
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
    scatter(unclampTrials(:,1),trials{unclampTrials(:,1),"nAnticipatoryLicks"},100,unclampColor,'filled'); hold on
    scatter(clampTrials(:,1),trials{clampTrials(:,1),"nAnticipatoryLicks"},100,clampColor,'filled'); hold on
    xlabel("Trials"); ylabel("Anticipatory licks"); box off
    title("Anticipatory licks (all trials)");
    
    
    % 4. Distribution of reward rection time
    ort_unclamp = trials{unclampTrials(:,1),["TrialNumber","OutcomeReactionTime"]};
    ort_clamp = trials{clampTrials(:,1),["TrialNumber","OutcomeReactionTime"]};
    ort_unclamp_early = ort_unclamp(ort_unclamp(:,1)<stageCutoff(1),2);
    ort_clamp_early = ort_clamp(ort_clamp(:,1)<stageCutoff(1),2);
    ort_unclamp_mid = ort_unclamp(ort_unclamp(:,1)>=stageCutoff(1) & ort_unclamp(:,1)<stageCutoff(2),2);
    ort_clamp_mid = ort_clamp(ort_clamp(:,1)>=stageCutoff(1) & ort_clamp(:,1)<stageCutoff(2),2);
    ort_unclamp_late = ort_unclamp(ort_unclamp(:,1)<=stageCutoff(2),2);
    ort_clamp_late = ort_clamp(ort_clamp(:,1)<=stageCutoff(2),2);
    
    % For early stage of session
    % Calculate bootstrap distribution
    stages = {ort_unclamp_early;ort_clamp_early};
    colors = [unclampColor;clampColor];
    nexttile;
    for s = 1:size(stages,1)
        if isempty(stages{s}) || isscalar(stages{s}); continue; end
        [~,bootsam] = bootstrp(nboot,[],stages{s});
        bs = stages{s}(bootsam) / params.sync.behaviorFs;
    
        h = histogram(bs,nBins); 
        h.FaceColor = colors(s,:); h.EdgeColor = colors(s,:); hold on
    end
    xlabel('Outcome reaction time (s)'); ylabel('Count'); box off
    title("Outcome reaction time (early stage)");
    
    % For middle stage of session
    % Calculate bootstrap distribution
    stages = {ort_unclamp_mid;ort_clamp_mid};
    colors = [unclampColor;clampColor];
    nexttile;
    for s = 1:size(stages,1)
        if isempty(stages{s}) || isscalar(stages{s}); continue; end
        [~,bootsam] = bootstrp(nboot,[],stages{s});
        bs = stages{s}(bootsam) / params.sync.behaviorFs;
    
        h = histogram(bs,nBins); 
        h.FaceColor = colors(s,:); h.EdgeColor = colors(s,:); hold on
    end
    xlabel('Outcome reaction time (s)'); ylabel('Count'); box off
    title("Outcome reaction time (middle stage)");
    
    % For late stage of session
    % Calculate bootstrap distribution
    stages = {ort_unclamp_late;ort_clamp_late};
    colors = [unclampColor;clampColor];
    nexttile;
    for s = 1:size(stages,1)
        if isempty(stages{s}) || isscalar(stages{s}); continue; end
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
    scatter(unclampTrials(:,1),trials{unclampTrials(:,1),"OutcomeReactionTime"}/params.sync.behaviorFs,100,unclampColor,'filled'); hold on
    scatter(clampTrials(:,1),trials{clampTrials(:,1),"OutcomeReactionTime"}/params.sync.behaviorFs,100,clampColor,'filled'); hold on
    xlabel("Trials"); ylabel("Outcome reaction time (s)"); box off
    title("Outcome reaction time (all trials)");
    
    saveFigures(gcf,'Summary_distributions',sessionpath,savePDF=true);

end
disp('Finished: all plots and struct are plotted and saved!');
return

end