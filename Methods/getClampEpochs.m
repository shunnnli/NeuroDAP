function [clampEpochs,clampOn,clampTime_sec] = getClampEpochs(timeSeries,options)
%GETCLAMPEPOCHS Stretches of session time when the clamp command was on.
%
% Mirrors analyzeSessions_clamp: either stored percentage command at or above
% threshold_pct counts as clamping, and OFF gaps shorter than gapFill_sec are
% bridged, because the command pulses within a single clamp epoch.
%
% clampEpochs is [start,end] in seconds since the first NI sample, one row per
% epoch. clampOn is the state of each stored clamp sample and clampTime_sec the
% matching session time, so a caller can resample the state onto another signal.

arguments
    timeSeries struct
    options.threshold_pct (1,1) double {mustBeFinite,mustBePositive} = 1
    options.gapFill_sec (1,1) double {mustBeFinite,mustBeNonnegative} = 0.1
end

names = string({timeSeries.name});
blueIdx = find(strcmpi(names,'blueClamp'),1);
redIdx = find(strcmpi(names,'redClamp'),1);
if isempty(blueIdx) || isempty(redIdx)
    error('getClampEpochs:MissingClamp', ...
        'Expected blueClamp and redClamp in timeSeries. Rerun loadSessions with withClamp=true.');
end

clampFs = double(timeSeries(blueIdx).finalFs);
clampOn = double(timeSeries(blueIdx).data(:)) >= options.threshold_pct | ...
          double(timeSeries(redIdx).data(:)) >= options.threshold_pct;

gapSamples = round(options.gapFill_sec*clampFs);
if gapSamples > 0
    offEdges = diff([false;~clampOn;false]);
    offStarts = find(offEdges == 1); offEnds = find(offEdges == -1)-1;
    for k = 1:numel(offStarts)
        % Only internal gaps, so an epoch running to either end stays bounded
        if offEnds(k)-offStarts(k)+1 <= gapSamples && offStarts(k) > 1 && offEnds(k) < numel(clampOn)
            clampOn(offStarts(k):offEnds(k)) = true;
        end
    end
end

clampTime_sec = (0:numel(clampOn)-1)'/clampFs + double(timeSeries(blueIdx).time_offset);
clampEdges = diff([false;clampOn;false]);
clampOnsets = find(clampEdges == 1);
clampOffsets = min(find(clampEdges == -1)-1,numel(clampOn));
clampEpochs = [clampTime_sec(clampOnsets),clampTime_sec(clampOffsets)+1/clampFs];
end
