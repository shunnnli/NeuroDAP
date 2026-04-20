function [traces,timestamp] = getTraces(eventIdx,timeRange,signal,options)

arguments
    eventIdx double
    timeRange double
    signal double
    
    options.params struct

    options.sameSystem logical = false

    options.eventSystem string = 'ni'
    options.signalSystem string = 'lj'
    options.signalFs double = nan % Fs of extracted signal
    options.syncFs double = nan % Fs of event system / native sync axis
    
    options.rmmissing logical = true % remove nan rows
    options.baseline double = 0
end

%% Check inputs

if strcmpi(options.eventSystem, options.signalSystem)
    options.sameSystem = true;
end

params = options.params;

%% Initialize event and signal timing
signalFs = options.signalFs;
nativeSignalFs = nan;
timeRef = [];
timeTarget = [];

if options.sameSystem
    if ~strcmp(options.eventSystem,options.signalSystem)
        options.eventSystem = options.signalSystem;
    end

    nativeSignalFs = inferSystemFs(params, options.signalSystem);

    % If caller provides signalFs, always use it for extracted traces.
    % If not provided, fall back to the native Fs stored in params.
    if isnan(signalFs)
        if ~isnan(nativeSignalFs)
            signalFs = nativeSignalFs;
        else
            warning('getTraces: did not provide signalFs and could not infer it from params, set signalFs = 50');
            signalFs = 50;
        end
    end

    % syncFs describes the event index space, not the extracted signal space.
    if isnan(options.syncFs)
        if ~isnan(nativeSignalFs)
            options.syncFs = nativeSignalFs;
        else
            options.syncFs = signalFs;
        end
    end
else
    if ~isfield(options,'params')
        error('getTraces: params is required if sameSystem=false');
    end

    % 0. Define event system
    if ~isfield(options,'eventSystem')
        if ~isfield(params.session,'baselineSystem'); options.eventSystem = 'ni';
        else; options.eventSystem = params.session.baselineSystem; end
    end

    % added by Emily Ferenczi 7/29/25 to make compatible with news rig
    if isfield(options,'eventSystem')
        if isfield(params,'session') && isfield(params.session,'baselineSystem')
            if ~strcmp(params.session.baselineSystem,options.eventSystem)
                options.eventSystem = params.session.baselineSystem; 
            end
        end
    end
    
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
    else
        error('getTraces: unsupported eventSystem %s', options.eventSystem);
    end
    
    % 1.2 Define signal system time and native signal Fs
    if strcmpi(options.signalSystem,'ni')
        timeTarget = params.sync.timeNI; 
        nativeSignalFs = inferSystemFs(params, options.signalSystem);
        if isnan(signalFs)
            if ~isnan(nativeSignalFs); signalFs = nativeSignalFs;
            else; signalFs = 50; end
            disp(['     getTraces: signalFs not provided, set to ',num2str(signalFs)]);
        end
    elseif strcmpi(options.signalSystem,'lj') || strcmpi(options.signalSystem,'labjack')
        timeTarget = params.sync.timePhotometry;
        nativeSignalFs = inferSystemFs(params, options.signalSystem);
        if isnan(signalFs)
            if ~isnan(nativeSignalFs); signalFs = nativeSignalFs;
            else; signalFs = 50; end
            disp(['     getTraces: signalFs not provided, set to ',num2str(signalFs)]);
        end
    elseif contains(options.signalSystem,'cam','IgnoreCase',true)
        timeTarget = params.sync.timeCamera;
        nativeSignalFs = inferSystemFs(params, options.signalSystem);
        if isnan(signalFs)
            if ~isnan(nativeSignalFs); signalFs = nativeSignalFs;
            else; signalFs = 100; end
        end
    elseif strcmpi(options.signalSystem,'imec')
        timeTarget = params.sync.timeImec;
        nativeSignalFs = inferSystemFs(params, options.signalSystem);
        if isnan(signalFs)
            if ~isnan(nativeSignalFs); signalFs = nativeSignalFs;
            else; signalFs = 30000; end
        end
    else
        error('getTraces: unsupported signalSystem %s', options.signalSystem);
    end

    % syncFs describes the event/reference axis, not the extracted signal.
    if isnan(options.syncFs)
        eventFs = inferSystemFs(params, options.eventSystem);
        if ~isnan(eventFs)
            options.syncFs = eventFs;
        elseif ~isempty(timeRef) && numel(timeRef) > 1
            dt = median(diff(timeRef),'omitnan');
            if ~isnan(dt) && dt > 0; options.syncFs = 1/dt; end
        end
    end
end

%% Convert events into seconds relative to signal start
if ~strcmp(options.eventSystem,options.signalSystem) || ~options.sameSystem
    % Different systems: first map event indices to the common sync timebase,
    % then convert that time into seconds relative to the signal start.
    if isempty(timeRef) || isempty(timeTarget)
        error('getTraces: missing sync time vectors for different-system alignment');
    end

    if isintegerlike(eventIdx)
        eventSample = round(eventIdx);
        validEvent = eventSample >= 1 & eventSample <= numel(timeRef);
        eventSyncSec = nan(size(eventSample));
        eventSyncSec(validEvent) = timeRef(eventSample(validEvent));
    else
        % eventIdx already provided in seconds on the shared sync axis
        eventSyncSec = eventIdx;
    end

    signalStartSec = timeTarget(1);
    eventInSec = eventSyncSec - signalStartSec;
else
    % Same system: if the signal was downsampled, use the native Fs from params
    % for event indices, but still extract traces at the provided signalFs.
    if isnan(options.syncFs)
        options.syncFs = signalFs;
    end

    if isintegerlike(eventIdx) && ~isnan(nativeSignalFs) && ~isnan(signalFs) && abs(signalFs - nativeSignalFs) > 0.5
        eventInSec = eventIdx / nativeSignalFs;
    else
        eventInSec = eventIdx / options.syncFs;
    end
end

timestamp = timeRange(1):(1/signalFs):timeRange(2);
traces = zeros(length(eventInSec),length(timestamp));

%% Generate traces
for i = 1:length(eventInSec)
    % If first bin calculated is before recording start, reset to 1
    if eventInSec(i)+timeRange(1) < 0 || isnan(eventInSec(i))
        traces(i,:) = nan(1,length(timestamp));
        continue
    end

    firstBin = round((eventInSec(i)+timestamp(1))*signalFs);
    if firstBin == 0; firstBin = 1; end
    lastBin = firstBin + length(timestamp) - 1;

    if any(options.baseline ~= 0)
        firstBin_baseline = round((eventInSec(i) + options.baseline(1))*signalFs);
        lastBin_baseline = round((eventInSec(i) + options.baseline(2))*signalFs);
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

function tf = isintegerlike(x)
tf = isnumeric(x) && all(isfinite(x(:))) && all(abs(x(:) - round(x(:))) < 1e-9);
end

function fs = inferSystemFs(params, systemName)
fs = nan;
if ~isfield(params,'sync')
    return
end

s = lower(char(systemName));

if strcmp(s,'ni')
    if isfield(params.sync,'behaviorFs'); fs = scalarize(params.sync.behaviorFs); return; end
    if isfield(params.sync,'ni_photometryFs'); fs = scalarize(params.sync.ni_photometryFs); return; end
elseif strcmp(s,'lj') || strcmp(s,'labjack') || strcmp(s,'photometry')
    if isfield(params.sync,'photometryFs'); fs = scalarize(params.sync.photometryFs); return; end
    if isfield(params.sync,'labjackFs'); fs = scalarize(params.sync.labjackFs); return; end
elseif contains(s,'cam')
    if isfield(params.sync,'camFs'); fs = scalarize(params.sync.camFs); return; end
elseif strcmp(s,'imec')
    if isfield(params.sync,'apFs'); fs = scalarize(params.sync.apFs); return; end
elseif strcmp(s,'lfp')
    if isfield(params.sync,'lfpFs'); fs = scalarize(params.sync.lfpFs); return; end
end
end

function x = scalarize(v)
if isempty(v)
    x = nan;
elseif isscalar(v)
    x = double(v);
else
    x = double(mode(v(:)));
end
end
