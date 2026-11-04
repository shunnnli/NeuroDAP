function stats = dmdControlMeanStats(traces,outputFs,endMs)
% Extrema of the selected-sweep control mean, timed relative to stimulus.
arguments
    traces double
    outputFs (1,1) double
    endMs (1,1) double
end
if isempty(traces)
    stats = struct('max',NaN,'min',NaN,'maxTime',NaN,'minTime',NaN,'auc',NaN);
    return
end
timeMs = (-(size(traces,2)-1):0)*1000/outputFs + endMs;
windowMs = [timeMs(1) timeMs(end)];
stats = dmdSweepMeanStats(traces,timeMs,windowMs, ...
    dmdDefaultCriteria().peakAverageMs,windowMs);
end
