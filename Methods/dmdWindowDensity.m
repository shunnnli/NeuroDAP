function [chargePcPerMs, travelPaPerMs] = dmdWindowDensity(trace,timeMs,windowMs,derivativeLagMs)

% Rectified charge density (pC/ms) and path-length density (pA/ms) of one
% trace over one window.
%
% Both are per ms of window so windows of unequal length stay comparable, and
% both are meant for the sweep average rather than single sweeps: averaging
% first cancels noise as sqrt(n) while rectifying first keeps all of it.
%
% The charge de-levels by the window's own median before rectifying, so only
% what moves inside the window counts. Without that, a few pA of drift is
% picked up on every sample, and because the control windows always sit later
% than the response window, drift would inflate the control side and bias the
% excess negative. The path length measures movement rather than position, so
% a constant offset contributes nothing to it either way.
%
% Mirrors charge_density_pc_per_ms() and derivative_density_pa_per_ms() in
% pyNeuroDAP's Shun_dmd_analysis notebook.

arguments
    trace double
    timeMs double
    windowMs (1,2) double
    derivativeLagMs (1,1) double
end

chargePcPerMs = NaN;
travelPaPerMs = NaN;

duration = windowMs(2) - windowMs(1);
if duration <= 0; return; end

timeMs = timeMs(:)';
trace = trace(:)';
values = trace(timeMs >= windowMs(1) & timeMs <= windowMs(2));
if isempty(values) || all(isnan(values)); return; end

sampleIntervalMs = median(diff(timeMs),'omitnan');

deLevelled = values - median(values,'omitnan');
chargePcPerMs = sum(abs(deLevelled),'omitnan') * sampleIntervalMs / 1000 / duration;

lag = max(1,round(derivativeLagMs/sampleIntervalMs));
if numel(values) > lag
    travelPaPerMs = sum(abs(values(1+lag:end)-values(1:end-lag)),'omitnan') / duration;
end

end
