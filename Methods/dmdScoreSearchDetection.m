function detection = dmdScoreSearchDetection(searchCurrentMap,searchRepMap,searchQC,options)

% Score one search's spots with the notebook's detection metrics.
%
% Called from loadSlicesDMD once a search's depths have been aggregated, where
% a spot's repeated sweeps are still together with the QC rows of the sweeps
% they came from. Two passes per depth: build the null that depth's spots are
% measured against, then score each spot against it.
%
% A depth too sparse for a MAD borrows the whole search's spots. That fallback
% only ever fires on coarse search depths the notebook does not analyze, since
% the repeated final depth always carries far more than minSpotsForReference
% spots -- and where it could fire on a repeated search, that search holds a
% single depth, so the wider pool is the same set and the answer is unchanged.
%
% searchCurrentMap{depth}{spot} is reps x samples; searchRepMap{depth}{spot}
% names the sweep each row came from; searchQC{depth} is the struct of columns
% localSummarizeDMDQC returns. options carries timeRange, controlWindowLength,
% baselineGuard, analysisWindowLength and maxRs.

arguments
    searchCurrentMap cell
    searchRepMap cell
    searchQC cell
    options (1,1) struct
end

criteria = dmdDefaultCriteria();
criteria.maxRs = options.maxRs;
baselineMs = [-options.controlWindowLength, -options.baselineGuard];

nDepth = numel(searchCurrentMap);
emptyStats = struct('max',NaN,'min',NaN,'maxTime',NaN,'minTime',NaN,'auc',NaN);
detection = struct('chargeExcessZ',{cell(nDepth,1)},'chargeExcessPc',{cell(nDepth,1)}, ...
                   'derivativeExcessZ',{cell(nDepth,1)},'reanalysisHotspot',{cell(nDepth,1)}, ...
                   'nSweepsKept',{cell(nDepth,1)},'filtered',{cell(nDepth,1)});

nSamples = localSampleCount(searchCurrentMap);
if nSamples == 0; return; end
timeMs = linspace(options.timeRange(1),options.timeRange(2),nSamples);

% Drop sweeps failing the series-resistance limit before anything is averaged.
filteredByDepth = cell(nDepth,1);
for depthIdx = 1:nDepth
    spotTraces = searchCurrentMap{depthIdx};
    sweepNames = searchRepMap{depthIdx};
    nSpots = numel(spotTraces);

    filtered = cell(nSpots,1);
    kept = zeros(nSpots,1);
    for spotIdx = 1:nSpots
        traces = spotTraces{spotIdx};
        if isempty(traces)
            filtered{spotIdx} = zeros(0,nSamples);
        else
            keep = localKeepRowsByRs(sweepNames{spotIdx},searchQC{depthIdx},criteria.maxRs);
            filtered{spotIdx} = traces(keep,:);
        end
        kept(spotIdx) = size(filtered{spotIdx},1);
    end

    filteredByDepth{depthIdx} = filtered;
    detection.nSweepsKept{depthIdx} = kept;
end

allUsable = vertcat(filteredByDepth{:});
allUsable = allUsable(cellfun(@(x) size(x,1) > 0, allUsable));
searchRef = [];

for depthIdx = 1:nDepth
    spots = filteredByDepth{depthIdx};
    usable = spots(cellfun(@(x) size(x,1) > 0, spots));

    try
        if numel(usable) >= criteria.minSpotsForReference
            ref = dmdHotspotVerdict('reference',usable,timeMs,criteria);
        else
            if isempty(searchRef)
                searchRef = dmdHotspotVerdict('reference',allUsable,timeMs,criteria);
            end
            ref = searchRef;
        end
    catch ME
        if strcmp(ME.identifier,'dmdDetectionWindows:noControlWindow')
            error('dmdScoreSearchDetection:storedWindowTooShort',...
                  ['timeRange is [%g %g] ms, which leaves no room for a detection control '...
                   'window. Widen it (the notebook expects [-30 100]) and reload.'],...
                  options.timeRange(1),options.timeRange(2));
        end
        rethrow(ME);
    end

    nSpots = numel(spots);
    chargeZ = nan(nSpots,1);
    chargePc = nan(nSpots,1);
    derivativeZ = nan(nSpots,1);
    reanalysis = false(nSpots,1);
    filteredStats = repmat(emptyStats,nSpots,1);

    for spotIdx = 1:nSpots
        traces = spots{spotIdx};
        scored = dmdHotspotVerdict('score',traces,timeMs,criteria,ref);
        chargeZ(spotIdx) = scored.chargeExcessZ;
        chargePc(spotIdx) = scored.chargeExcessPc;
        derivativeZ(spotIdx) = scored.derivativeExcessZ;

        if size(traces,1) > 0
            reanalysis(spotIdx) = dmdHotspotVerdict('reanalysis',traces,timeMs,criteria).isHotspot;
            filteredStats(spotIdx) = dmdSweepMeanStats(traces,timeMs, ...
                [0 options.analysisWindowLength],criteria.peakAverageMs,baselineMs);
        end
    end

    detection.chargeExcessZ{depthIdx} = chargeZ;
    detection.chargeExcessPc{depthIdx} = chargePc;
    detection.derivativeExcessZ{depthIdx} = derivativeZ;
    detection.reanalysisHotspot{depthIdx} = reanalysis;
    detection.filtered{depthIdx} = filteredStats;
end

end


function nSamples = localSampleCount(searchCurrentMap)

% The stored window is timeRange, so the sample count fixes the time base.
nSamples = 0;
for depthIdx = 1:numel(searchCurrentMap)
    for spotIdx = 1:numel(searchCurrentMap{depthIdx})
        if ~isempty(searchCurrentMap{depthIdx}{spotIdx})
            nSamples = size(searchCurrentMap{depthIdx}{spotIdx},2);
            return
        end
    end
end

end


function keep = localKeepRowsByRs(sweepNames,qc,maxRs)

% Join each stored row to its sweep's QC row and apply the Rs limit. A row
% whose sweep has no QC entry is kept: no evidence is not evidence of a bad
% patch.

keep = true(numel(sweepNames),1);
if isempty(sweepNames) || isempty(qc) || ~isfield(qc,'Sweep'); return; end

qcSweeps = string(qc.Sweep);
for r = 1:numel(sweepNames)
    idx = find(qcSweeps == sweepNames(r),1);
    if isempty(idx); continue; end

    rsHeader = NaN;
    if isfield(qc,'Rs_headerString'); rsHeader = qc.Rs_headerString(idx); end
    keep(r) = dmdFilterSweepsByRs(qc.Rs(idx),rsHeader,maxRs);
end

end
