function out = dmdHotspotVerdict(mode,traces,timeMs,criteria,reference)

% Score one spot's repeated sweeps the way pyNeuroDAP's Shun_dmd_analysis
% notebook does.
%
%   out = dmdHotspotVerdict('reanalysis',traces,timeMs,criteria)
%       The offline re-test (HOTSPOT_TYPE = "reanalysis"), judged on the
%       sweep average rather than sweep by sweep.
%
% traces is reps x samples; timeMs is the matching time base in ms, zero at
% stimulus onset. criteria comes from dmdDefaultCriteria.
%
% Note the reanalysis test runs on the stored trace as it is, with no further
% levelling, which is what the notebook does: reclassify_hotspots hands
% test_trace_for_hotspot the cached traces directly.

arguments
    mode (1,1) string
    traces double
    timeMs double
    criteria (1,1) struct
    reference (1,1) struct = struct()
end

timeMs = timeMs(:)';

switch mode
    case "reanalysis"
        out = localReanalysis(mean(traces,1,'omitnan'),timeMs,criteria);
    otherwise
        error('dmdHotspotVerdict:unknownMode','Unknown mode "%s".',mode);
end

end


function out = localReanalysis(meanTrace,timeMs,criteria)

out = struct('baselineSd',NaN,'threshold',NaN,'nOverThreshold',0, ...
             'peakAbs',NaN,'isHotspot',false);

meanTrace = meanTrace(:)';
baseline = meanTrace(timeMs >= criteria.baselineMs(1) & timeMs <= criteria.baselineMs(2));
window = meanTrace(timeMs >= criteria.reanalysisWindowMs(1) & timeMs <= criteria.reanalysisWindowMs(2));
baseline = baseline(isfinite(baseline));
window = window(isfinite(window));
if numel(baseline) < 2 || isempty(window); return; end

out.baselineSd = std(baseline,0);
out.threshold = criteria.thresholdFactor * out.baselineSd;
out.nOverThreshold = sum(abs(window) >= out.threshold);
out.peakAbs = max(abs(window));
out.isHotspot = out.baselineSd > 0 ...
    && out.nOverThreshold >= criteria.minSamples ...
    && out.peakAbs >= criteria.minAmplitudePa;

end
