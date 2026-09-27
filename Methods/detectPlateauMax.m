function level = detectPlateauMax(values,options)

% Largest level that is actually held by the signal, ie the highest plateau
% rather than the highest single sample. A level only qualifies if at least
% minSamples samples sit within tolerance of it, so isolated noise spikes are
% stepped over. Returns NaN if no level qualifies (or values is empty).
%
% Used to find the full scale of a clamp command, either in ADC counts
% (voltage2percent) or in percent (rescaleClampTraces).
%
% Example:
%   fullScale = detectPlateauMax(adc_counts,minSamples=5,tolerance=5.75);

arguments
    values double
    options.minSamples double = 5 % samples required at the detected level
    options.tolerance double = 0  % how close to the level still counts as at it
end

values = values(isfinite(values));
level = NaN;
if isempty(values); return; end

% Spike rejection needs a trace longer than the sample requirement
minSamples = max(1,round(options.minSamples));
if numel(values) < minSamples; minSamples = 1; end

level = max(values);
while sum(values >= level - options.tolerance) < minSamples
    remaining = values(values < level - options.tolerance);
    if isempty(remaining); level = NaN; return; end
    level = max(remaining);
end

end
