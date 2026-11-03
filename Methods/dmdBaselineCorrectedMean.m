function meanTrace = dmdBaselineCorrectedMean(traces,timeMs,baselineMs)

% The sweep average, after removing each sweep's own baseline median.
%
% Per sweep rather than per spot: a sweep that sat at a different holding
% level contributes its offset to the average otherwise, and the charge
% metrics are only invariant to an offset shared by every sweep.
%
% Mirrors baseline_corrected_mean() in pyNeuroDAP's Shun_dmd_analysis
% notebook.

arguments
    traces double
    timeMs double
    baselineMs (1,2) double
end

timeMs = timeMs(:)';

if isempty(traces)
    meanTrace = nan(1,numel(timeMs));
    return
end

mask = timeMs >= baselineMs(1) & timeMs <= baselineMs(2);
if any(mask)
    level = median(traces(:,mask),2,'omitnan');
    level(~isfinite(level)) = 0;
else
    level = zeros(size(traces,1),1);
end

meanTrace = mean(traces - level,1,'omitnan');

end
