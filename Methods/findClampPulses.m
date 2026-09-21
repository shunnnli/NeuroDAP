function [events,detection] = findClampPulses(blueClamp,redClamp,Fs,options)
%FINDCLAMPPULSES Find open-loop clamp power sweeps in the NI command traces.
% Detects the calibration pulses commanded on the blue and red clamp DAC
% channels and groups them into power levels. Command voltages are converted
% to percent of the configured DAC range, so power is percent of range and
% not physical optical power. Detection runs on a decimated copy of the
% commands (targetFs, 200 Hz by default); event indices are mapped back to
% the original NI samples (1-based, offset exclusive) so they can be used
% with plotTraces/getTraces. Long blocks are excluded and exactly maxDuration
% pulses are included. Nearby measured powers are grouped within
% powerTolerance percentage points, without merging sparse but distinct
% levels. Pulses clipped by the start or end of the recording are still
% returned, with complete = false.
%
% Example:
%   [events,detection] = findClampPulses(blueClamp,redClamp,params.sync.behaviorFs, ...
%       redClampRange=[25 500],blueClampRange=[800 1600],laserTime=params.sync.timeNI);

arguments
    blueClamp double {mustBeVector,mustBeNonempty,mustBeFinite}
    redClamp double {mustBeVector,mustBeNonempty,mustBeFinite}
    Fs (1,1) double {mustBePositive,mustBeFinite}
    options.blueClampRange (1,2) double {mustBeFinite} = [800,1600]
    options.redClampRange (1,2) double {mustBeFinite} = [25,500]
    options.targetFs (1,1) double {mustBePositive,mustBeFinite} = 200
    options.thresholdPct (1,1) double {mustBeNonnegative,mustBeLessThan(options.thresholdPct,100)} = 5
    options.minDuration (1,1) double {mustBePositive,mustBeFinite} = 0.5
    options.maxDuration (1,1) double {mustBePositive} = 5
    options.minGap (1,1) double {mustBeNonnegative,mustBeFinite} = 2
    options.powerTolerance (1,1) double {mustBeNonnegative,mustBeFinite} = 2
    options.laserTime double = [] % shared sync clock for the NI samples
end

if numel(blueClamp) ~= numel(redClamp)
    error('findClampPulses:LengthMismatch','Laser traces must have equal lengths.');
end
if options.blueClampRange(2) <= options.blueClampRange(1) || ...
        options.redClampRange(2) <= options.redClampRange(1)
    error('findClampPulses:InvalidRange','Each DAC range must have an increasing minimum and maximum.');
end
if options.maxDuration < options.minDuration
    error('findClampPulses:InvalidDuration','maxDuration must be at least minDuration.');
end
hasSync = ~isempty(options.laserTime);
if hasSync
    laserTime = options.laserTime;
    if ~isvector(laserTime) || numel(laserTime) ~= numel(blueClamp) || ...
            numel(laserTime) < 2 || any(~isfinite(laserTime(:))) || any(diff(laserTime(:)) <= 0)
        error('findClampPulses:InvalidSync', ...
            'laserTime must contain one finite, strictly increasing timestamp per NI sample.');
    end
end

stride = max(1,round(Fs/options.targetFs));
analysisFs = Fs/stride;
sampleIdx = 1:stride:numel(blueClamp);
red = voltage2percent(redClamp(sampleIdx),options.redClampRange);
blue = voltage2percent(blueClamp(sampleIdx),options.blueClampRange);
redEvents = extractPulses(red,analysisFs,"red","excite",options);
blueEvents = extractPulses(blue,analysisFs,"blue","inhibit",options);
events = sortrows([redEvents;blueEvents],'onset_analysis_idx');
events.onset_idx = (events.onset_analysis_idx-1)*stride+1;
events.offset_idx = min((events.offset_analysis_idx-1)*stride+1,numel(blueClamp)+1);
events.onset_sec = (events.onset_idx-1)/Fs;
if hasSync
    events.onset_sync_sec = reshape(laserTime(events.onset_idx),[],1);
end

% Group nearby measured powers into the commanded levels.
events.level_pct = nan(height(events),1);
channels = ["red","blue"];
for c = 1:numel(channels)
    rows = find(events.channel == channels(c));
    [powers,order] = sort(events.power_pct(rows));
    rows = rows(order);
    % Bound the full spread of each group to avoid chaining nearby levels.
    first = 1;
    while first <= numel(rows)
        last = first;
        while last < numel(rows) && powers(last+1)-powers(first) <= options.powerTolerance
            last = last+1;
        end
        eventRows = sort(rows(first:last));
        events.level_pct(eventRows) = round(median(events.power_pct(eventRows)));
        first = last+1;
    end
end

% Sync vectors live in sync_*.mat; do not duplicate them in the outputs.
options = rmfield(options,'laserTime');
detection = struct('originalFs',Fs,'analysisFs',analysisFs,'stride',stride, ...
    'nSamples',numel(blueClamp),'options',options);
if isempty(events)
    warning('findClampPulses:NoPulses','No calibration pulses met the power and duration criteria.');
end
end

function events = extractPulses(power,Fs,channel,type,options)
active = power(:) > options.thresholdPct;
edges = diff([false;active;false]);
onsets = find(edges == 1);
offsets = find(edges == -1);
starts = zeros(0,1); stops = zeros(0,1); powers = zeros(0,1);
k = 1;
while k <= numel(onsets)
    on = onsets(k); off = offsets(k);
    while k < numel(onsets) && onsets(k+1)-off <= round(options.minGap*Fs)
        k = k+1;
        off = offsets(k);
    end
    duration = (off-on)/Fs;
    if duration >= options.minDuration && duration <= options.maxDuration+1/Fs
        starts(end+1,1) = on; %#ok<AGROW>
        stops(end+1,1) = off; %#ok<AGROW>
        powers(end+1,1) = mean(power(on:off-1)); %#ok<AGROW>
    end
    k = k+1;
end
complete = starts > 1 & stops <= numel(power);
events = table(repmat(channel,numel(starts),1),repmat(type,numel(starts),1), ...
    starts,stops,(stops-starts)/Fs,powers,complete, ...
    'VariableNames',{'channel','type','onset_analysis_idx','offset_analysis_idx', ...
    'duration_sec','power_pct','complete'});
end
