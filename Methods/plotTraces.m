function [traces,timestamp] = plotTraces(eventIdx,timeRange,signal,color,params,options)

%% Notes
% plotTraces.m, Shun Li
% 2023/09/02: change to generalize to event or signal streams of choice
% 2023/09/05: fixed minor error to make it applicable for camera and
% downsampled LJ photometry

%% Function
arguments
    eventIdx double
    timeRange double
    signal double
    color
    params struct

    options.eventSystem string = 'ni'
    options.signalSystem string = 'lj'
    options.signalFs double = nan

    options.rmmissing logical = true % remove nan rows
    options.smooth double = 0; % 0: no smoothing, else is samples of smooth data
    options.smoothMethod string = 'movmean';
    options.baselineWindow double = 0 % subtracted by x secs before event time

    % Plot options
    options.plot logical = true % whether or not to plot traces
    options.LineStyle (1,1) string = "-"
    options.LineWidth (1,1) {mustBeNumeric} = 2
    options.plotShuffled logical = false
    options.shuffledColor double = [.75, .75, .75]
    options.nShuffle double = 10000

end

% 0. Define event system
if ~isfield(params.session,'baselineSystem'); options.eventSystem = 'ni';
else; options.eventSystem = params.session.baselineSystem; end

% 1. Initialize event and signal system time
% 1.1 Define event system time
if strcmpi(options.eventSystem,'ni')
    timeRef = params.sync.timeNI; 
elseif strcmpi(options.eventSystem,'lj') || strcmpi(options.eventSystem,'labjack')
    timeRef = params.sync.timePhotometry; 
elseif contains(options.eventSystem,'cam','IgnoreCase',true)
    timeRef = params.sync.timeCamera; 
elseif strcmpi(options.eventSystem,'imec')
    timeRef = params.sync.timeImec; 
end
% 1.2 Define signal system time
if strcmpi(options.signalSystem,'ni')
    timeTarget = params.sync.timeNI; 
    if ~isnan(options.signalFs); signalFs = options.signalFs;
    else
        if isfield(params.sync,'ni_photometryFs'); signalFs = params.sync.ni_photometryFs;
        else; signalFs = 50; end
        disp(['     plotTraces: signalFs not provided, set to ',num2str(signalFs)]);
    end
    syncFs = params.sync.behaviorFs;
elseif strcmpi(options.signalSystem,'lj') || strcmpi(options.signalSystem,'labjack')
    timeTarget = params.sync.timePhotometry;
    if ~isnan(options.signalFs); signalFs = options.signalFs;
    else
        if isfield(params.sync,'photometryFs'); signalFs = mode(params.sync.photometryFs);
        else; signalFs = 50; end
        disp(['     plotTraces: signalFs not provided, set to ',num2str(signalFs)]);
    end
    syncFs = params.sync.labjackFs;
elseif contains(options.signalSystem,'cam','IgnoreCase',true)
    timeTarget = params.sync.timeCamera;
    if ~isnan(options.signalFs); signalFs = options.signalFs;
    else
        if isfield(params.sync,'camFs'); signalFs = params.sync.camFs;
        else; signalFs = 100; end
    end
    syncFs = params.sync.camFs;
elseif strcmpi(options.signalSystem,'imec')
    timeTarget = params.sync.timeImec;
    if ~isnan(options.signalFs); signalFs = options.signalFs;
    else
        if isfield(params.sync,'apFs'); signalFs = params.sync.apFs;
        else; signalFs = 30000; end
    end
    syncFs = params.sync.apFs;
end

% 2. Get signal traces around a specific event
if ~strcmp(options.eventSystem,options.signalSystem)
    % If event system and signal system is different
    eventInSignal = findCorrespondingTime(eventIdx,timeRef,timeTarget);
    eventInSec = eventInSignal / syncFs;
else
    % If event system and signal system is the same
    eventInSec = eventIdx / syncFs;
end
timestamp = timeRange(1):(1/signalFs):timeRange(2);
traces = zeros(length(eventInSec),length(timestamp));

% 3. generate traces
for i = 1:length(eventInSec)
    % If first bin calculated is before recording start, reset to 1
    if (floor(eventInSec(i)+timeRange(1)) <= 0) || isnan(eventInSec(i))
        traces(i,:) = nan(1,length(timestamp));
        continue
    end

    firstBin = floor((eventInSec(i)+timestamp(1))*signalFs);
    lastBin = firstBin + length(timestamp) - 1;
    eventBin = floor(eventInSec(i)*signalFs);
    if options.baselineWindow ~= 0
        firstBin_baseline = floor((eventInSec(i) + options.baselineWindow(1))*signalFs);
        lastBin_baseline = floor((eventInSec(i) + options.baselineWindow(2))*signalFs);
        baseline = mean(signal(firstBin_baseline:lastBin_baseline));
    else
        baseline = 0;
    end
    
    if lastBin > length(signal)
        trace = nan(1,length(timestamp));
    else
        trace = signal(firstBin:lastBin) - baseline;
    end
    
    traces(i,:) = trace;
end
% 3.1 Remove missing if necessary
if options.rmmissing; traces = rmmissing(traces); end

% 4. plot traces
if options.plot
    % 4.1 Smooth data if neccessary
    if options.smooth ~= 0
        if options.smooth == 1; smoothWindow = 0.1*signalFs;
        else; smoothWindow = options.smooth;
        end
    else
        smoothWindow = 0;
    end

    % 4.2 Plot shuffled data
    if options.plotShuffled
        shuffled = shuffleTraces(traces);
        plotSEM(timestamp,shuffled,options.shuffledColor,smooth=smoothWindow,smoothMethod=options.smoothMethod,...
        LineStyle=options.LineStyle,LineWidth=options.LineWidth); hold on
    end

    % 4.3 Plot SEM
    plotSEM(timestamp,traces,color,smooth=smoothWindow,smoothMethod=options.smoothMethod,...
        LineStyle=options.LineStyle,LineWidth=options.LineWidth);
    xlabel('Time (s)'); ylabel('z-score'); hold on
end

end

