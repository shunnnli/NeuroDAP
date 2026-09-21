function calibration = analyzeClampCalibration(photometry_raw,blueClamp,redClamp,Fs,options)
%ANALYZECLAMPCALIBRATION Open-loop power sweeps from NI voltage recordings.
% Follows PID-tuning.ipynb: 200 Hz sampling, rounded 10-bit photometry ADC,
% 20 ms EMA, a 1 s pre-onset F0, and median dF/F0 in the last 2 s of each pulse.
% Power is percent of the configured DAC range, not physical optical power.
% Event indices refer to the original NI samples (1-based, offset exclusive).
% Long blocks are excluded; exactly 5 s pulses are included. Nearby measured
% powers are grouped within powerTolerance percentage points, without merging
% sparse but distinct levels. Boundary pulses remain in events as invalid.

arguments
    photometry_raw double {mustBeVector,mustBeNonempty,mustBeFinite}
    blueClamp double {mustBeVector,mustBeNonempty,mustBeFinite}
    redClamp double {mustBeVector,mustBeNonempty,mustBeFinite}
    Fs (1,1) double {mustBePositive,mustBeFinite}
    options.blueClampRange (1,2) double {mustBeFinite} = [800,1600]
    options.redClampRange (1,2) double {mustBeFinite} = [25,500]
    options.targetFs (1,1) double {mustBePositive,mustBeFinite} = 200
    options.emaTau (1,1) double {mustBeNonnegative,mustBeFinite} = 0.020
    options.preTime (1,1) double {mustBePositive,mustBeFinite} = 1
    options.postTime (1,1) double {mustBeNonnegative,mustBeFinite} = 2
    options.lateTime (1,1) double {mustBePositive,mustBeFinite} = 2
    options.thresholdPct (1,1) double {mustBeNonnegative,mustBeLessThan(options.thresholdPct,100)} = 5
    options.minDuration (1,1) double {mustBePositive,mustBeFinite} = 0.5
    options.maxDuration (1,1) double {mustBePositive} = 5
    options.minGap (1,1) double {mustBeNonnegative,mustBeFinite} = 2
    options.powerTolerance (1,1) double {mustBeNonnegative,mustBeFinite} = 2
end

if numel(photometry_raw) ~= numel(blueClamp) || numel(photometry_raw) ~= numel(redClamp)
    error('analyzeClampCalibration:LengthMismatch','NI photometry and laser traces must have equal lengths.');
end
if options.blueClampRange(2) <= options.blueClampRange(1) || ...
        options.redClampRange(2) <= options.redClampRange(1)
    error('analyzeClampCalibration:InvalidRange','Each DAC range must have an increasing minimum and maximum.');
end
if options.maxDuration < options.minDuration
    error('analyzeClampCalibration:InvalidDuration','maxDuration must be at least minDuration.');
end

stride = max(1,round(Fs/options.targetFs));
analysisFs = Fs/stride;
sampleIdx = 1:stride:numel(photometry_raw);
signal = round(voltage2arduino(photometry_raw(sampleIdx)));
signal = signal(:);
if options.emaTau > 0
    alpha = 1-exp(-1/(analysisFs*options.emaTau));
    signal = filter(alpha,[1,-(1-alpha)],signal,(1-alpha)*signal(1));
end
red = voltage2percent(redClamp(sampleIdx),options.redClampRange);
blue = voltage2percent(blueClamp(sampleIdx),options.blueClampRange);
redEvents = extractPulses(red,analysisFs,"red","excite",options);
blueEvents = extractPulses(blue,analysisFs,"blue","inhibit",options);
events = sortrows([redEvents;blueEvents],'onset_analysis_idx');
events.onset_idx = (events.onset_analysis_idx-1)*stride+1;
events.offset_idx = min((events.offset_analysis_idx-1)*stride+1,numel(photometry_raw)+1);
events.onset_sec = (events.onset_idx-1)/Fs;
events.level_pct = nan(height(events),1);
events.validTrial = false(height(events),1);

preSamples = max(1,round(options.preTime*analysisFs));
postSamples = round(options.postTime*analysisFs);
lateSamples = max(1,round(options.lateTime*analysisFs));
groups = struct('channel',{},'type',{},'power_pct',{},'eventRows',{}, ...
    'nTrials',{},'time_sec',{},'duration_sec',{},'raw_adc',{},'baseline_adc',{}, ...
    'dff',{},'mean_dff',{},'sem_dff',{},'lateOnMedian_dff',{},'medianResponse_dff',{});

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
        level = round(median(events.power_pct(eventRows)));
        events.level_pct(eventRows) = level;
        on = events.onset_analysis_idx(eventRows);
        off = events.offset_analysis_idx(eventRows);
        onSamples = off-on;
        offsets = -preSamples:(max(onSamples)+postSamples);
        raw = nan(numel(on),numel(offsets));
        dff = raw;
        f0 = nan(numel(on),1);
        late = f0;
        for trial = 1:numel(on)
            % A full pre-window and complete stimulation are needed for F0
            % and the late-on statistic. Missing post-window data stays NaN.
            if ~events.complete(eventRows(trial)) || on(trial)-preSamples < 1
                continue
            end
            idx = on(trial)+offsets;
            available = idx >= 1 & idx <= numel(signal);
            raw(trial,available) = signal(idx(available));
            f0(trial) = mean(raw(trial,offsets < 0));
            if ~isfinite(f0(trial)) || abs(f0(trial)) < 1e-9
                continue
            end
            dff(trial,:) = (raw(trial,:)-f0(trial))/f0(trial);
            lateWindow = offsets >= max(0,onSamples(trial)-lateSamples) & offsets < onSamples(trial);
            late(trial) = median(dff(trial,lateWindow),'omitnan');
            events.validTrial(eventRows(trial)) = true;
        end
        counts = sum(isfinite(dff),1);
        sem = std(dff,0,1,'omitnan')./sqrt(counts);
        sem(counts == 0) = NaN;
        g = numel(groups)+1;
        groups(g) = struct('channel',channels(c),'type',events.type(eventRows(1)), ...
            'power_pct',level,'eventRows',eventRows, ...
            'nTrials',sum(events.validTrial(eventRows)),'time_sec',offsets/analysisFs, ...
            'duration_sec',onSamples/analysisFs,'raw_adc',raw,'baseline_adc',f0, ...
            'dff',dff,'mean_dff',mean(dff,1,'omitnan'),'sem_dff',sem, ...
            'lateOnMedian_dff',late,'medianResponse_dff',median(late,'omitnan'));
        first = last+1;
    end
end

fits = struct('channel',{},'slope',{},'intercept',{},'r2',{});
for c = 1:numel(channels)
    selected = arrayfun(@(group) group.channel == channels(c),groups);
    x = [groups(selected).power_pct];
    y = [groups(selected).medianResponse_dff];
    finite = isfinite(x) & isfinite(y);
    x = x(finite); y = y(finite);
    slope = NaN; intercept = NaN; r2 = NaN;
    if numel(unique(x)) >= 2
        p = polyfit(x,y,1);
        slope = p(1); intercept = p(2);
        total = sum((y-mean(y)).^2);
        if total > 0; r2 = 1-sum((y-polyval(p,x)).^2)/total; end
    end
    fits(c) = struct('channel',channels(c),'slope',slope,'intercept',intercept,'r2',r2);
end

calibration = struct('source','photometry_raw (NI)','originalFs',Fs, ...
    'analysisFs',analysisFs,'options',options,'events',events,'groups',groups,'fits',fits);
if isempty(events)
    warning('analyzeClampCalibration:NoPulses','No calibration pulses met the power and duration criteria.');
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
