function [traces,timestamp] = getTraces(eventIdx,timeRange,signal,params,options)

arguments
    eventIdx double
    timeRange double
    signal double
    params struct

    % options.getTraces logical = true
    options.eventSystem string = 'ni'
    options.signalSystem string = 'lj'
    options.signalFs double = nan

    options.rmmissing logical = true % remove nan rows

    options.baseline double = 0
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
    % eventBin = floor(eventInSec(i)*signalFs);
    if options.baseline ~= 0
        firstBin_baseline = floor((eventInSec(i) + options.baseline(1))*signalFs);
        lastBin_baseline = floor((eventInSec(i) + options.baseline(2))*signalFs);
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


end