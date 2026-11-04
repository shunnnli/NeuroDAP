function stats = dmdSweepMeanStats(traces,timeMs,windowMs,peakAverageMs,baselineMs)

% Peak, latency and area of one spot's baseline-corrected sweep MEAN.
%
% This is the notebook's definition, and it is not the same quantity as
% loadSlicesDMD's own stats.response.max/min: those average a peak window on
% each sweep and then average across repetitions, where this takes the
% extremes of the averaged trace. Averaging first cancels noise as sqrt(n),
% so the two differ by more than rounding on a quiet spot.
%
% Each extremum is averaged over peakAverageMs centred on its own sample, so
% one noisy sample cannot set an amplitude. Times are in ms from stimulus
% onset.
%
% baselineMs defaults to dmdDefaultCriteria().baselineMs, but loadSlicesDMD's
% control window follows controlWindowLength and baselineGuard, so it is
% passed explicitly from there.

arguments
    traces double
    timeMs double
    windowMs (1,2) double
    peakAverageMs (1,1) double
    baselineMs (1,2) double = dmdDefaultCriteria().baselineMs
end

stats = struct('max',NaN,'min',NaN,'maxTime',NaN,'minTime',NaN,'auc',NaN);

timeMs = timeMs(:)';
meanTrace = dmdBaselineCorrectedMean(traces,timeMs,baselineMs);

mask = timeMs >= windowMs(1) & timeMs <= windowMs(2);
values = meanTrace(mask);
windowTime = timeMs(mask);
if isempty(values) || all(isnan(values)); return; end

sampleIntervalMs = median(diff(timeMs),'omitnan');
stats.auc = sum(values,'omitnan') * sampleIntervalMs / 1000;

[~,maxIdx] = max(values);
[~,minIdx] = min(values);

stats.maxTime = windowTime(maxIdx);
stats.minTime = windowTime(minIdx);
% Locate extrema inside the response window, but average on the full trace,
% as the notebook does (including when a peak lies at a window boundary).
stats.max = localAverageAround(meanTrace,timeMs,stats.maxTime,peakAverageMs);
stats.min = localAverageAround(meanTrace,timeMs,stats.minTime,peakAverageMs);

end


function value = localAverageAround(values,timeMs,peakTime,averageMs)

half = averageMs/2;
mask = timeMs >= peakTime-half & timeMs <= peakTime+half;
value = mean(values(mask),'omitnan');

end
