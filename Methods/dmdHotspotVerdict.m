function out = dmdHotspotVerdict(mode,traces,timeMs,criteria,reference)

% Score one spot's repeated sweeps the way pyNeuroDAP's Shun_dmd_analysis
% notebook does.
%
%   out = dmdHotspotVerdict('reanalysis',traces,timeMs,criteria)
%       The offline re-test (HOTSPOT_TYPE = "reanalysis"), judged on the
%       sweep average rather than sweep by sweep.
%
%   ref = dmdHotspotVerdict('reference',spotTraces,timeMs,criteria)
%       The null level and spread every spot of one depth is scored against,
%       from a cell array holding each spot's reps x samples matrix. One spot
%       contributes a single control reading and so cannot give a spread at
%       all; taking the between-spot spread as a property of the depth is the
%       assumption that makes the z possible. It averages the statistic, never
%       the traces: averaging control traces across spots would cancel their
%       noise as sqrt(nSpots) and leave a reference near zero.
%
%   out = dmdHotspotVerdict('score',traces,timeMs,criteria,reference)
%       That spot's charge and path-length excesses, as pC/pA over the window
%       and as a z on the reference's spread.
%
% traces is reps x samples; timeMs is the matching time base in ms, zero at
% stimulus onset. criteria comes from dmdDefaultCriteria.
%
% Note the reanalysis test runs on the stored trace as it is, with no further
% levelling, which is what the notebook does: reclassify_hotspots hands
% test_trace_for_hotspot the cached traces directly. The detection metrics do
% level each sweep first, via dmdBaselineCorrectedMean.

arguments
    mode (1,1) string
    traces                                  % reps x samples, or a cell array of them
    timeMs double
    criteria (1,1) struct
    reference (1,1) struct = struct()
end

timeMs = timeMs(:)';

switch mode
    case "reanalysis"
        out = localReanalysis(mean(traces,1,'omitnan'),timeMs,criteria);
    case "reference"
        out = localReference(traces,timeMs,criteria);
    case "score"
        out = localScore(traces,timeMs,criteria,reference);
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


function ref = localReference(spotTraces,timeMs,criteria)

if ~iscell(spotTraces)
    error('dmdHotspotVerdict:referenceNeedsCell', ...
          'Mode "reference" takes a cell array holding each spot''s traces.');
end

ref = struct('chargeCentre',NaN,'chargeSpread',NaN, ...
             'derivativeCentre',NaN,'derivativeSpread',NaN,'nSpots',0);

[~,controlWindows] = dmdDetectionWindows(timeMs,criteria.detectionWindowMs,criteria.detectionGapMs);

nSpots = numel(spotTraces);
chargeReadings = nan(nSpots,1);
travelReadings = nan(nSpots,1);
for s = 1:nSpots
    meanTrace = dmdBaselineCorrectedMean(spotTraces{s},timeMs,criteria.baselineMs);
    [chargeReadings(s),travelReadings(s)] = ...
        localControlDensities(meanTrace,timeMs,controlWindows,criteria);
end

[ref.chargeCentre,ref.chargeSpread,ref.nSpots] = ...
    localCentreAndSpread(chargeReadings,criteria.minSpotsForReference);
[ref.derivativeCentre,ref.derivativeSpread] = ...
    localCentreAndSpread(travelReadings,criteria.minSpotsForReference);

end


function out = localScore(traces,timeMs,criteria,reference)

required = {'chargeCentre','chargeSpread','derivativeCentre','derivativeSpread'};
if ~all(isfield(reference,required))
    error('dmdHotspotVerdict:scoreNeedsReference', ...
          'Mode "score" needs a reference from dmdHotspotVerdict(''reference'',...).');
end

out = struct('chargeExcessPc',NaN,'chargeExcessZ',NaN, ...
             'derivativeExcessPa',NaN,'derivativeExcessZ',NaN, ...
             'nSweeps',size(traces,1));

windowMs = criteria.detectionWindowMs;
duration = windowMs(2) - windowMs(1);
meanTrace = dmdBaselineCorrectedMean(traces,timeMs,criteria.baselineMs);
[chargeDensity,travelDensity] = dmdWindowDensity(meanTrace,timeMs,windowMs,criteria.derivativeLagMs);

[out.chargeExcessPc,out.chargeExcessZ] = localExcess( ...
    chargeDensity,reference.chargeCentre,reference.chargeSpread,duration);
[out.derivativeExcessPa,out.derivativeExcessZ] = localExcess( ...
    travelDensity,reference.derivativeCentre,reference.derivativeSpread,duration);

end


function [q,g] = localControlDensities(meanTrace,timeMs,controlWindows,criteria)

nWindows = size(controlWindows,1);
charges = nan(nWindows,1);
travels = nan(nWindows,1);
for w = 1:nWindows
    [charges(w),travels(w)] = ...
        dmdWindowDensity(meanTrace,timeMs,controlWindows(w,:),criteria.derivativeLagMs);
end

q = median(charges,'omitnan');
g = median(travels,'omitnan');

end


function [centre,spread,nUsable] = localCentreAndSpread(readings,minSpots)

readings = readings(isfinite(readings));
nUsable = numel(readings);

centre = NaN; spread = NaN;
if nUsable == 0; return; end

centre = median(readings);
if nUsable >= minSpots
    spread = 1.4826 * median(abs(readings - centre));
end

end


function [excess,z] = localExcess(density,centre,spread,duration)

excess = (density - centre) * duration;

z = NaN;
if isfinite(spread) && spread > 0
    z = (density - centre) / spread;
end

end
