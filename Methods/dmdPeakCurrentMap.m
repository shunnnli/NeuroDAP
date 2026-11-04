function responseMap = dmdPeakCurrentMap(spotTraces,keptRows,spotLocations,mapSize,timeMs,windowMs,baselineMs,metric)
% Signed IPSC-EPSC current of each QC-filtered, baseline-corrected spot mean.
% max + min subtracts the magnitude of the negative EPSC from the IPSC.
% AUC mode integrates that same mean, with equal weight per retained sweep.
% Unlike a hotspot verdict, a measured amplitude needs only one usable sweep.
% Unsampled spots and spots with no usable sweeps remain NaN.

arguments
    spotTraces cell
    keptRows cell
    spotLocations
    mapSize (1,2) double
    timeMs double
    windowMs (1,2) double
    baselineMs (1,2) double
    metric (1,1) string {mustBeMember(metric,["netPeak","peakToPeak","auc"])} = "netPeak"
end

if iscell(spotLocations); spotLocations = spotLocations{1}; end
responseMap = nan(mapSize);
criteria = dmdDefaultCriteria();
for k = 1:min(numel(spotTraces),size(spotLocations,1))
    location = spotLocations(k,:);
    if any(~isfinite(location)) || isempty(spotTraces{k}); continue; end
    traces = spotTraces{k}(keptRows{k},:);
    if isempty(traces); continue; end
    stats = dmdSweepMeanStats(traces,timeMs,windowMs,criteria.peakAverageMs,baselineMs);
    cols = max(1,floor(location(1))):min(mapSize(2),ceil(location(2)));
    rows = max(1,floor(location(3))):min(mapSize(1),ceil(location(4)));
    if metric == "netPeak"
        responseMap(rows,cols) = stats.max + stats.min;
    elseif metric == "peakToPeak"
        responseMap(rows,cols) = stats.max - stats.min;
    else
        responseMap(rows,cols) = stats.auc;
    end
end
end
