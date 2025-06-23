function analyzeSessions_forSally(sessionpath,options)

arguments
    sessionpath string
    
    options.task
    options.outputName string %name of the output file in subfolder of the session

    options.analyzeTraces logical = true

    options.redo logical = true % Recalculate trial table and all preprocessing
    options.round logical = false % Round reward/airpuff/tone to get duration data
    
    options.plotPhotometry logical = true % Plot photometry summary plot
    options.plotBehavior logical = true % Plot lick raster summary plot

    options.shortTimeRange double = [-1,5]
    options.longTimeRange double = [-5,15]
end

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


% Define baselineSystem
if ~isfield(params.session,'baselineSystem')
    params.session.baselineSystem = 'NI';
end

save(strcat(sessionpath,filesep,'sync_',options.outputName),'params','-append');

%% Load behavior labels

load(strcat(sessionpath,filesep,'labeled_output.mat'));

labeledData.Duration = labeledData.EndFrame - labeledData.StartFrame;

% Clean up potential mistakes
for r = 1:height(labeledData)
    % IF negative duration, switch start and end frame
    if labeledData(r,:).Duration < 0
        temp = labeledData(r,:).StartFrame;
        labeledData(r,:).StartFrame = labeledData(r,:).EndFrame;
        labeledData(r,:).EndFrame = temp;
    end
end
labeledData.Duration = labeledData.EndFrame - labeledData.StartFrame;
behaviors = unique(labeledData.Label);


%% Find the number of photometry channels
photometryIdx = find(cellfun(@(x) contains(x,["NI","LJ"],"IgnoreCase",true), {timeSeries.system}));
photometryName = cellfun(@(x) unique(x,'rows'), {timeSeries(photometryIdx).name},'UniformOutput',false);
nSignals = length(photometryIdx);
disp(['Finished: found ', num2str(nSignals),' photometry signals']);

%% Plot behavior events as single trials

% Plotting configs
colorListIdx = round(linspace(1,500,numel(behaviors)));

for photometry = 1:nSignals
    path = photometryIdx(photometry);
    signal = timeSeries(path).data;
    finalFs = timeSeries(path).finalFs;
    system = timeSeries(path).system;

    initializeFig(1,1); tiledlayout(1,numel(behaviors));
    for b = 1:numel(behaviors)
        nexttile; axis off;
        % Calculate the longest instance in sec
        maxDuration = max(labeledData.Duration(strcmpi(labeledData.Label,behaviors{b})));
        maxDurationSec = maxDuration / params.sync.camFs;
        color = bluePurpleRed(colorListIdx(b),:);
        eventIdx = labeledData.StartFrame(strcmpi(labeledData.Label,behaviors{b}));

        [~,~] = plotTraces(eventIdx,[-2,maxDurationSec],...
                                signal,color,params,eventSystem='cam',...
                                signalFs=finalFs,signalSystem=system,...
                                plotIndividual=true);
        plotEvent(behaviors{b},0,color=color);
        xlim([-2,maxDurationSec]);
        legend([behaviors{b}, ' (n=',num2str(length(eventIdx)),')']);
    end
end

%% Find & plot optogenetic stim

blueLaserON = find(blueLaser); blueLaserInterval = [0,diff(blueLaserON)];
blueStimIdx = blueLaserON(blueLaserInterval > 10*params.sync.behaviorFs);

blueStimStartIdx = blueLaserInterval > 10*params.sync.behaviorFs;
blueStimDuration = 40 - round(blueLaserInterval(blueStimStartIdx)/params.sync.behaviorFs); % in sec
durationList = unique(blueStimDuration);
colorListIdx = round(linspace(1,500,numel(durationList)));

% Plot
for photometry = 1:nSignals
    path = photometryIdx(photometry);
    signal = timeSeries(path).data;
    finalFs = timeSeries(path).finalFs;
    system = timeSeries(path).system;

    initializeFig(1,1); tiledlayout(1,numel(durationList));
    for d = 1:numel(durationList)
        nexttile;
        color = bluePurpleRed(colorListIdx(d),:);
        cur_duration = durationList(d);
        blueStimIdx_duration = blueStimIdx(blueStimDuration == cur_duration);
        [~,~] = plotTraces(blueStimIdx_duration,[-2,40],...
                            signal,color,params,eventSystem='ni',...
                            signalFs=finalFs,signalSystem=system,...
                            plotIndividual=true);
        plotEvent(['Opto stim (', num2str(cur_duration),' sec)'],cur_duration,color=color);
        xlim([-2,40]);
        legend(['Opto stim (n=',num2str(length(blueStimIdx_duration)),')']);
    end
end

%% Scatter plot of grooming events aligned to start of opto-stim



end