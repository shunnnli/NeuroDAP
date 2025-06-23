function analyzeSessions_toneRamp(sessionpath,options)

%% TODO:
% make sure lick stuff is correct

arguments
    sessionpath string
    
    options.task
    options.outputName string %name of the output file in subfolder of the session
    
    options.redStim logical = true
    options.redStimPattern cell
    options.blueStimPattern cell
    options.includeOtherStim logical = false % consider stim in other color as a single trials

    options.shortCueDuration double = 2
    options.longCueDuration double = 4

    options.pavlovian logical = false
    options.reactionTime double = 1
    options.minLicks double = 2 % min licks to get reward
    options.combineOmission logical = true % combine omission trials

    options.analyzeTraces logical = true

    options.redo logical = true % Recalculate trial table and all preprocessing
    options.round logical = false % Round reward/airpuff/tone to get duration data
    options.performing logical = false % Only plot traces where the animal performs
    
    options.plotPhotometry logical = true % Plot photometry summary plot
    options.plotBehavior logical = true % Plot lick raster summary plot

    options.lick_binSize double = 0.1
    options.shortTimeRange double = [-1,8]
    options.longTimeRange double = [-5,10]
end

%% Notes
% Shun_analyzeSessions_toneRamp
% Shun Li, 6/13/2025

%% Load data

[~,~,~,~,~,~,bluePurpleRed] = loadColors;
             
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

% Load stim patterns
if isfield(options,'redStimPattern') && (options.redo || ~isfield(params,'stim'))
    params.stim.pulseFreq_red = str2double(options.redStimPattern{1}); 
    params.stim.pulseDuration_red = str2double(options.redStimPattern{2}); 
    params.stim.stimDuration_red = str2double(options.redStimPattern{3});
    params.stim.nPulsesPerStim_red = (params.stim.stimDuration_red/1000) * params.stim.pulseFreq_red;
    params.stim.pulseInterval_red = (1/params.stim.pulseFreq_red) - params.stim.pulseDuration_red;
elseif ~isfield(options,'redStimPattern') && ~isfield(params,'stim') 
    params.stim.pulseFreq_red = 50; 
    params.stim.pulseDuration_red = 5; 
    params.stim.stimDuration_red = 500;
    params.stim.nPulsesPerStim_red = (params.stim.stimDuration_red/1000) * params.stim.pulseFreq_red;
    params.stim.pulseInterval_red = (1000/params.stim.pulseFreq_red) - params.stim.pulseDuration_red;
end

fn = fieldnames(params.stim); hasBlue = any(contains(fn, 'blue')); 
if isfield(options,'blueStimPattern') && (options.redo || ~hasBlue)
    params.stim.pulseFreq_blue = str2double(options.blueStimPattern{1}); 
    params.stim.pulseDuration_blue = str2double(options.blueStimPattern{2}); 
    params.stim.stimDuration_blue = str2double(options.blueStimPattern{3});
    params.stim.nPulsesPerStim_blue = (params.stim.stimDuration_blue/1000) * params.stim.pulseFreq_blue;
    params.stim.pulseInterval_blue = (1/params.stim.pulseFreq_blue) - params.stim.pulseDuration_blue;
elseif ~isfield(options,'blueStimPattern') && ~isfield(params,'stim') 
    params.stim.pulseFreq_blue = 30; 
    params.stim.pulseDuration_blue = 10; 
    params.stim.stimDuration_blue = 500;
    params.stim.nPulsesPerStim_blue = (params.stim.stimDuration_blue/1000) * params.stim.pulseFreq_blue;
    params.stim.pulseInterval_blue = (1000/params.stim.pulseFreq_blue) - params.stim.pulseDuration_blue;
end

if options.redStim
    params.stim.color = 'red';
    disp(['Stim params: color = ', params.stim.color]);
    disp(['Stim params: pulseFreq = ',num2str(params.stim.pulseFreq_red)]);
    disp(['Stim params: pulseDuration = ',num2str(params.stim.pulseDuration_red)]);
    disp(['Stim params: stimDuration = ',num2str(params.stim.stimDuration_red)]);
    disp(['Stim params: nPulsesPerStim = ',num2str(params.stim.nPulsesPerStim_red)]);
    disp(['Stim params: pulseInterval = ',num2str(params.stim.pulseInterval_red)]);
else
    params.stim.color = 'blue'; 
    disp(['Stim params: color = ', params.stim.color]);
    disp(['Stim params: pulseFreq = ',num2str(params.stim.pulseFreq_blue)]);
    disp(['Stim params: pulseDuration = ',num2str(params.stim.pulseDuration_blue)]);
    disp(['Stim params: stimDuration = ',num2str(params.stim.stimDuration_blue)]);
    disp(['Stim params: nPulsesPerStim = ',num2str(params.stim.nPulsesPerStim_blue)]);
    disp(['Stim params: pulseInterval = ',num2str(params.stim.pulseInterval_blue)]);
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
options.shortCueDuration = options.shortCueDuration + 1;
options.longCueDuration = options.longCueDuration + 1;
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
    options.shortCueDuration = params.analyze.shortCueDuration;
    options.longCueDuration = params.analyze.longCueDuration;
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


% Find start of opto cue
blueStimIdx = find(blueLaser);
redStimIdx = find(redLaser);

if strcmpi(params.stim.color,'red')
    cueStimPulses = redStimIdx; 
    nPulsesPerStim_cue = params.stim.nPulsesPerStim_red;
    stimDuration_cue = params.stim.stimDuration_red;
    otherStimPulses = blueStimIdx;
    nPulsesPerStim_other = params.stim.nPulsesPerStim_blue;
    stimDuration_other = params.stim.stimDuration_blue;
elseif strcmpi(params.stim.color,'blue')
    cueStimPulses = blueStimIdx; 
    nPulsesPerStim_cue = params.stim.nPulsesPerStim_blue;
    stimDuration_cue = params.stim.stimDuration_blue;
    otherStimPulses = redStimIdx;
    nPulsesPerStim_other = params.stim.nPulsesPerStim_red;
    stimDuration_other = params.stim.stimDuration_red;
end
    
if ~exist('optoStim1','var') || options.redo
    if  ~isempty(find(cueStimPulses, 1))
        if nPulsesPerStim_cue > 1
            % Save first pulse data for stim
            intervalThreshold = 10000;
            temp_interval = [100000,diff(cueStimPulses)];
            optoStim1 = cueStimPulses(temp_interval > intervalThreshold);
            save(strcat(sessionpath,filesep,'timeseries_',options.outputName),"optoStim1",'-append');
        else
            optoStim1 = cueStimPulses;
        end
    else
        optoStim1 = [];
    end

    if  ~isempty(find(otherStimPulses, 1))
        if nPulsesPerStim_other > 1
            % Save first pulse data for other color stim
            intervalThreshold = 10000;
            temp_interval = [100000,diff(otherStimPulses)];
            optoStim2 = otherStimPulses(temp_interval > intervalThreshold);
            save(strcat(sessionpath,filesep,'timeseries_',options.outputName),"optoStim2",'-append');
        else
            optoStim2 = otherStimPulses;
        end
    else
        optoStim2 = [];
    end
    save(strcat(sessionpath,filesep,'timeseries_',options.outputName),"optoStim1","optoStim2",'-append');
end


% Find start of lick bout
rightLickON = find(rightLick);
if ~exist('lickBout','var') || options.redo
    % Get lick bout start time (ILI < 0.5s)
    lickBout = getLickBout(rightLickON);
    save(strcat(sessionpath,filesep,'timeseries_',options.outputName),"lickBout",'-append');
end

%% Find different cue

leftToneON = find(leftTone); rightToneON = find(rightTone);
tol        = round(0.005 * params.sync.behaviorFs);            % 5 ms in samples
[lia,locb] = ismembertol(leftToneON, rightToneON, tol);
sineTone   = leftToneON(lia);              % those within ±tol of a right
rampTone   = leftToneON(~lia);
jumpTone   = setdiff(rightToneON, rightToneON(locb(lia)));

v = leftTone(rampTone);
D = [abs(v-options.shortCueDuration); abs(v-options.longCueDuration)];      
[~, cls] = min(D, [], 1);      % cls(i)==1 means v(i) is nearer 3, ==2 nearer 6
rampToneShort = rampTone(cls==1);
rampToneLong = rampTone(cls==2);

v = rightTone(jumpTone);
D = [abs(v-options.shortCueDuration); abs(v-options.longCueDuration)];      
[~, cls] = min(D, [], 1);      % cls(i)==1 means v(i) is nearer 3, ==2 nearer 6
jumpToneShort = jumpTone(cls==1);
jumpToneLong = jumpTone(cls==2);

%% Combine stim&tone to form trial start
waterIdx = find(rightSolenoid_rounded);  
airpuffIdx = find(airpuff_rounded);

if ~exist('trials','var') || options.redo
    if contains(options.task,'punish')
        if options.includeOtherStim
            [allTrials,~] = getTrials(rampToneShort,rampToneLong,...
                                 jumpToneShort,jumpToneLong,sineTone,...
                                 optoStim1,optoStim2,waterIdx);
        else
            [allTrials,~] = getTrials(rampToneShort,rampToneLong,...
                                 jumpToneShort,jumpToneLong,sineTone,...
                                 optoStim1,waterIdx);
        end
    elseif contains(options.task,'reward')
        if options.includeOtherStim
            [allTrials,~] = getTrials(rampToneShort,rampToneLong,...
                                 jumpToneShort,jumpToneLong,sineTone,...
                                 optoStim1,optoStim2);
        else
            [allTrials,~] = getTrials(rampToneShort,rampToneLong,...
                                 jumpToneShort,jumpToneLong,sineTone,...
                                 optoStim1);
        end
    end
    save(strcat(sessionpath,filesep,'timeseries_',options.outputName),"allTrials",'-append');
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
    events{5} = or(leftTone,rightTone); 
    events{6} = optoStim1;      events{7} = optoStim2;
    events{8} = {rampToneShort,rampToneLong,jumpToneShort,jumpToneLong,sineTone};

    trials = getTrialTable_toneRamp(options.task,events,rightSolenoid_rounded,airpuff_rounded,...
                pavlovian=options.pavlovian,reactionTime=options.reactionTime);
    % Calculate performance cutoff
    if contains(options.task,'reward'); [trials,cutoff_sample] = getSessionCutoff(trials,"->reward");
    elseif contains(options.task,'punish'); [trials,cutoff_sample] = getSessionCutoff(trials,"->punish");
    else; [trials,cutoff_sample] = getSessionCutoff(trials,"random"); 
    end
    params.analysis.cutoff_sample = cutoff_sample;
    disp(['     Session cutoff calculated: ',num2str(cutoff_sample)]);

    % Save to behavior_.mat
    save(strcat(sessionpath,filesep,'behavior_',options.outputName),'trials','-append');
    disp('Finished: trial table saved');
end


%% Task specific params

% Select baseline idx (align to each baseline lick)
% baselineLicks = cell2mat(trials{~cellfun(@isempty,trials{1:end-1,'BaselineLicks'}),'BaselineLicks'});
% if isempty(baselineLicks); baselineIdx = round((trials{2:end,"CueTime"} - trials{2:end,"ENL"}) - 5*params.sync.behaviorFs);
% else; baselineIdx = baselineLicks(:,1); end
randomMinSample = floor(10*params.sync.behaviorFs);
randomMaxSample = length(params.sync.timeNI) - randomMinSample;
baselineIdx = randi([randomMinSample,randomMaxSample],100,1);

rampShortTrials = trials{strcmpi(trials.toneType,'ramp_short'),["TrialNumber","CueTime","OutcomeTime","ENL"]};
rampLongTrials = trials{strcmpi(trials.toneType,'ramp_long'),["TrialNumber","CueTime","OutcomeTime","ENL"]};
jumpShortTrials = trials{strcmpi(trials.toneType,'jump_short'),["TrialNumber","CueTime","OutcomeTime","ENL"]};
jumpLongTrials = trials{strcmpi(trials.toneType,'jump_long'),["TrialNumber","CueTime","OutcomeTime","ENL"]};
sineTrials = trials{strcmpi(trials.toneType,'sine'),["TrialNumber","CueTime","OutcomeTime","ENL"]};

stageTime = [-2,0;0,6;6,10];
analysisEvents = {waterIdx,waterLickIdx,...
    rampToneShort,rampToneLong,jumpToneShort,jumpToneLong,sineTone,...
    airpuffIdx,baselineIdx,optoStim2};
eventTrialNum = {findTrials(waterIdx,trials),findTrials(waterLickIdx,trials),...
                rampShortTrials(:,1),rampLongTrials(:,1),...
                jumpShortTrials(:,1),jumpLongTrials(:,1),...
                sineTrials(:,1),...
                findTrials(airpuffIdx,trials),findTrials(baselineIdx,trials),findTrials(optoStim2,trials)};
if strcmp(params.stim.color,'red')
    analysisLabels = {'Water','Rewarded licks',...
        'Ramp (short)','Ramp (Long)','Jump (short)','Jump (Long)','Sine',...
        'Airpuff','Baseline','Blue stim'};
else
    analysisLabels = {'Water','Rewarded licks',...
        'Ramp (short)','Ramp (Long)','Jump (short)','Jump (Long)','Sine',...
        'Airpuff','Baseline','Red stim'};
end
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
    'waterIdx','waterLickIdx','airpuffIdx','baselineIdx',...
    'rampToneShort','rampToneLong','jumpToneShort','jumpToneLong','sineTone',...
    'rampShortTrials','rampLongTrials','jumpShortTrials','jumpLongTrials','sineTrials',...
    '-append');

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

        %% Plot single stimulus PSTH
        
        eventIdxes = {rampToneShort,rampToneLong,jumpToneShort,jumpToneLong,sineTone,...
                        waterLickIdx,airpuffIdx,optoStim2};
        omissionIdxes = {[],[],[],[],[],[],[],[]};
        eventDurations = [options.shortCueDuration,options.longCueDuration,...
                          options.shortCueDuration,options.longCueDuration,...
                          options.shortCueDuration,...
                          0,0.02,stimDuration_other/1000];
        longTimeRange = options.longTimeRange;
        shortTimeRange = options.shortTimeRange;
        if strcmp(params.stim.color,'red')
            labels = {'Ramp (short)','Ramp (Long)','Jump (short)','Jump (Long)','Sine',...
                        'Water','Airpuff','Blue stim'};
            groupSizes = [10,10,10,10,10,10,10,30,30,20];
        else
            labels = {'Ramp (short)','Ramp (Long)','Jump (short)','Jump (Long)','Sine',...
                        'Water','Airpuff','Red stim'};
            groupSizes = [10,10,10,10,10,30,30,10];
        end
        
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
            legend({[label,' (n=',num2str(length(eventIdx)),')'],...
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
            legend({[label,' (n=',num2str(length(eventIdx)),')'],...
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
                xlabel('Trials'); ylabel([analysis(row).name,' signal area (z-score)']);
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
            saveas(gcf,strcat(sessionpath,filesep,'Analysis_',cur_signal,'_subtrial_area.png'));
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

    saveas(gcf,strcat(sessionpath,filesep,'Behavior_ITI&LickBout.png'));
end


if options.plotBehavior && contains(options.task,'pairing')     
    %% Plot session overview for licking
    timeRange = options.longTimeRange; cameraTimeRange = options.shortTimeRange;
    markerSize = 30;
    timeRange = [-3,10];

    % Get event time and number by trial type
    rampShortTrials(:,3) = rampShortTrials(:,3)./params.sync.behaviorFs;
    rampLongTrials(:,3) = rampLongTrials(:,3)./params.sync.behaviorFs;
    jumpShortTrials(:,3) = jumpShortTrials(:,3)./params.sync.behaviorFs;
    jumpLongTrials(:,3) = jumpLongTrials(:,3)./params.sync.behaviorFs;
    sineTrials(:,3) = sineTrials(:,3)./params.sync.behaviorFs;

    % getLicks by trial type
    [rampShortLickRate,~,rampShortLicks] = getLicks(timeRange,rampToneShort,options.lick_binSize,[],rightLick,...
                                params.sync.behaviorFs,params.sync.timeNI);
    [rampLongLickRate,~,rampLongLicks] = getLicks(timeRange,rampToneLong,options.lick_binSize,[],rightLick,...
                                params.sync.behaviorFs,params.sync.timeNI);
    [jumpShortLickRate,~,jumpShortLicks] = getLicks(timeRange,jumpToneShort,options.lick_binSize,[],rightLick,...
                                params.sync.behaviorFs,params.sync.timeNI);
    [jumpLongLickRate,~,jumpLongLicks] = getLicks(timeRange,jumpToneLong,options.lick_binSize,[],rightLick,...
                                params.sync.behaviorFs,params.sync.timeNI);
    [sineLickRate,~,sineLicks] = getLicks(timeRange,sineTone,options.lick_binSize,[],rightLick,...
                                params.sync.behaviorFs,params.sync.timeNI);

    % 1. Plot lick raster and trace
    initializeFig(0.67,1); tiledlayout(5,2);
    rampLongColor = bluePurpleRed(end,:); rampShortColor = addOpacity(rampLongColor,0.5);
    jumpLongColor = bluePurpleRed(150,:); jumpShortColor = addOpacity(jumpLongColor,0.5);
    sineColor = bluePurpleRed(350,:);

    % 1.1 Plot lick raster plot for session
    nexttile([5,1]);
    plot_licks = rampShortLicks; plot_trials = rampShortTrials; plot_color = rampShortColor;
    for i = 1:size(plot_licks,1)
        scatter(plot_licks{i},plot_trials(i,1),markerSize,'filled','MarkerFaceColor',plot_color); hold on
        scatter(plot_trials(i,3),plot_trials(i,1),markerSize+20,bluePurpleRed(1,:),'pentagram','filled'); hold on
    end
    plot_licks = rampLongLicks; plot_trials = rampLongTrials; plot_color = rampLongColor;
    for i = 1:size(plot_licks,1)
        scatter(plot_licks{i},plot_trials(i,1),markerSize,'filled','MarkerFaceColor',plot_color); hold on
        scatter(plot_trials(i,3),plot_trials(i,1),markerSize+20,bluePurpleRed(1,:),'pentagram','filled'); hold on
    end
    plot_licks = jumpShortLicks; plot_trials = jumpShortTrials; plot_color = jumpShortColor;
    for i = 1:size(plot_licks,1)
        scatter(plot_licks{i},plot_trials(i,1),markerSize,'filled','MarkerFaceColor',plot_color); hold on
        scatter(plot_trials(i,3),plot_trials(i,1),markerSize+20,bluePurpleRed(1,:),'pentagram','filled'); hold on
    end
    plot_licks = jumpLongLicks; plot_trials = jumpLongTrials; plot_color = jumpLongColor;
    for i = 1:size(plot_licks,1)
        scatter(plot_licks{i},plot_trials(i,1),markerSize,'filled','MarkerFaceColor',plot_color); hold on
        scatter(plot_trials(i,3),plot_trials(i,1),markerSize+20,bluePurpleRed(1,:),'pentagram','filled'); hold on
    end
    plot_licks = sineLicks; plot_trials = sineTrials; plot_color = sineColor;
    for i = 1:size(plot_licks,1)
        scatter(plot_licks{i},plot_trials(i,1),markerSize,'filled','MarkerFaceColor',plot_color); hold on
        scatter(plot_trials(i,3),plot_trials(i,1),markerSize+20,bluePurpleRed(1,:),'pentagram','filled'); hold on
    end
    xlabel('Time (s)'); xlim([timeRange(1),timeRange(2)]);
    ylabel('Trial'); ylim([0,size(trials,1)]);
    plotEvent("",0.5);
    % 1.2 Plot lick traces across session
    traces = {rampShortLickRate,rampLongLickRate,jumpShortLickRate,jumpLongLickRate,sineLickRate};
    labels = {'Ramp (short)','Ramp (Long)','Jump (short)','Jump (Long)','Sine'};
    groupSizes = [10,10,10,10,10];
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
    saveFigures(gcf,'Behavior_LickOverview',sessionpath,savePNG=true,savePDF=false);

end
disp('Finished: all plots and struct are plotted and saved!');
return

end