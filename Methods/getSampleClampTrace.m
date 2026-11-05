function sample = getSampleClampTrace(timeSeries,params,referenceIdx,timeRange,options)
%GETSAMPLECLAMPTRACE Extract one DA / clampIndex example from loaded signals.
% referenceIdx is a native NI sample (e.g. trials.CueTime). timeRange is in
% seconds relative to that sample. Uses getTraces for NI/LJ synchronization
% and getClampIdx for red minus blue. No additional DA normalization is used.
%
% DA / timestampDA and clampIndex / timestampClamp are row vectors. DA and
% clamp commands retain their stored sampling rates. blueClamp and redClamp
% contain the original percentage commands in the requested window.
% sessionTimeDA_sec and sessionTimeClamp_sec are seconds since the first NI
% sample; referenceTime_sec uses the same origin. removeClampBaseline=true
% subtracts each command's median in [-0.2,0) s, even if that baseline is
% outside the displayed window. Windows outside a recording are rejected.

arguments
    timeSeries struct
    params (1,1) struct
    referenceIdx (1,1) double {mustBeFinite,mustBeInteger,mustBePositive}
    timeRange (1,2) double {mustBeFinite}
    options.daIdx (1,1) double {mustBeInteger,mustBePositive} = 1
    options.removeClampBaseline (1,1) logical = true
end

if timeRange(2) <= timeRange(1)
    error('getSampleClampTrace:InvalidWindow','Window end must be greater than window start.');
end
if ~all(isfield(timeSeries,{'name','data','finalFs','system'})) || options.daIdx > numel(timeSeries)
    error('getSampleClampTrace:InvalidTimeSeries','Select a valid channel with name, data, finalFs and system.');
end
if ~isfield(params,'sync') || ~isfield(params.sync,'timeNI') || ...
        ~isfield(params.sync,'behaviorFs') || referenceIdx > numel(params.sync.timeNI)
    error('getSampleClampTrace:InvalidSync','Reference must be a sample in params.sync.timeNI with behaviorFs defined.');
end
if isfield(params,'session') && isfield(params.session,'baselineSystem') && ...
        ~strcmpi(params.session.baselineSystem,'ni')
    error('getSampleClampTrace:InvalidSync','Clamp examples require trials referenced to the NI behavior system.');
end
% getTraces defaults to the session baseline system for cross-system events.
params.session.baselineSystem = 'NI';
names = string({timeSeries.name});
blueIdx = find(strcmpi(names,'blueClamp'));
redIdx = find(strcmpi(names,'redClamp'));
if numel(blueIdx) ~= 1 || numel(redIdx) ~= 1
    error('getSampleClampTrace:MissingClamp', ...
        'Expected one blueClamp and one redClamp in timeSeries. Rerun loadSessions with withClamp=true.');
end
da = timeSeries(options.daIdx);
blue = timeSeries(blueIdx);
red = timeSeries(redIdx);
if contains(lower(string(da.name)),'clamp') || ...
        ~any(strcmpi(da.system,{'ni','lj','labjack'}))
    error('getSampleClampTrace:InvalidDA','Select a photometry channel, rather than a clamp or camera channel.');
end
for idx = [options.daIdx,blueIdx,redIdx]
    validateattributes(timeSeries(idx).data,{'numeric'},{'vector','real','nonempty'},'','timeSeries.data');
    validateattributes(timeSeries(idx).finalFs,{'numeric'},{'scalar','finite','positive'},'','timeSeries.finalFs');
end
if ~strcmpi(blue.system,'ni') || ~strcmpi(red.system,'ni') || blue.finalFs ~= red.finalFs
    error('getSampleClampTrace:ClampSamplingMismatch', ...
        'blueClamp and redClamp must be NI signals with the same finalFs. Rerun loadSessions.');
end

% DA and clamp commands can have different finalFs and recording start times.
[daTrace,daTime] = extractChannel(da,params,referenceIdx,timeRange);
clampRange = timeRange;
if options.removeClampBaseline
    clampRange = [min(timeRange(1),-0.2),max(timeRange(2),0)];
end
[blueTrace,clampTime] = extractChannel(blue,params,referenceIdx,clampRange);
[redTrace,~] = extractChannel(red,params,referenceIdx,clampRange);
if options.removeClampBaseline
    % Exclude zero despite floating-point noise in colon-generated timestamps.
    baseline = find(clampTime >= -0.2-1e-9 & clampTime < -1e-9);
    if isempty(baseline) || any(~isfinite(blueTrace(baseline))) || any(~isfinite(redTrace(baseline)))
        error('getSampleClampTrace:MissingBaseline','A finite clamp baseline in [-0.2,0) s is required.');
    end
    index = getClampIdx(blueTrace,redTrace,type='diff',baselineWindow=baseline);
else
    index = getClampIdx(blueTrace,redTrace,type='diff',removeBaseline=false);
end
keep = clampTime >= timeRange(1)-1e-9 & clampTime <= timeRange(2)+1e-9;

sample.referenceIdx = referenceIdx;
sample.referenceTime_sec = double(params.sync.timeNI(referenceIdx)-params.sync.timeNI(1));
sample.timeRange = timeRange;
sample.DAName = char(da.name);
sample.DAUnits = signalUnits(da);
sample.DASystem = char(da.system);
sample.DAFs = double(da.finalFs);
sample.DA = daTrace;
sample.timestampDA = daTime;
sample.sessionTimeDA_sec = sample.referenceTime_sec+daTime;
sample.clampFs = double(blue.finalFs);
sample.blueClamp = blueTrace(keep);
sample.redClamp = redTrace(keep);
sample.clampIndex = index(keep);
sample.timestampClamp = clampTime(keep);
sample.sessionTimeClamp_sec = sample.referenceTime_sec+sample.timestampClamp;
sample.removeClampBaseline = options.removeClampBaseline;
end

function [trace,timestamp] = extractChannel(channel,params,referenceIdx,timeRange)
[trace,timestamp] = getTraces(referenceIdx,timeRange,double(channel.data(:))', ...
    params=params,signalFs=double(channel.finalFs),signalSystem=string(channel.system), ...
    eventSystem='ni',rmmissing=false);
if isempty(trace) || ~any(isfinite(trace))
    error('getSampleClampTrace:WindowOutsideRecording', ...
        ['Window [%g, %g] s (including the baseline, if enabled) is unavailable for %s. ', ...
         'Choose a reference/window inside both recordings.'],timeRange(1),timeRange(2),char(channel.name));
end
end

function units = signalUnits(channel)
units = 'Stored signal';
if isfield(channel,'detrend_type') && ~isempty(channel.detrend_type)
    switch lower(string(channel.detrend_type))
        case 'rolling-z'
            units = 'z-score';
        case 'dff'
            units = '\DeltaF/F';
        case 'none'
            units = 'V';
        otherwise
            units = char(channel.detrend_type);
    end
end
end
