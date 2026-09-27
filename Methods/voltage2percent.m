function [pct,usedRange] = voltage2percent(voltage, range, options)

% Convert a clamp/laser command voltage into percentage of its command range.
%
% range = [min_count, max_count] in ADC counts. By default the nominal
% max_count is used. With autoMax=true, the actual full scale is detected
% from the data instead: the clamp command reaches the same maximum value at
% some point within every session (briefly or continuously), so the largest
% sustained command level in the trace IS 100%. This avoids the systematic
% error of normalizing by a nominal max that does not match the hardware
% (eg a command topping out at 500 counts read against [25,600] gives 83%).
%
% A level only counts as the full scale if at least autoMaxMinSamples samples
% sit within autoMaxTolerance counts of it, so single noise spikes are
% ignored. If the detected span is smaller than autoMaxMinSpan of the nominal
% span (eg the clamp was never turned on in this session), the nominal
% max_count is kept and a warning is thrown.
%
% Example:
%   redClamp_pct = voltage2percent(redClamp,[25,600],autoMax=true);

arguments
    voltage double
    range double
    options.Vref double = 5
    options.ADC_resolution double = 4095

    options.autoMax logical = false       % detect full scale from the data
    options.autoMaxMinSamples double = 5  % samples needed at the detected level
    options.autoMaxTolerance double       % counts; default is 1% of nominal span
    options.autoMaxMinSpan double = 0.5   % detected span / nominal span floor
end
    min_count = range(1);
    max_count = range(end);
    nominal_span = max_count - min_count;

    % Convert voltage to ADC counts
    adc_counts = (voltage / options.Vref) * options.ADC_resolution;

    % Detect the full scale command actually used in this session
    if options.autoMax
        if ~isfield(options,'autoMaxTolerance')
            options.autoMaxTolerance = 0.01 * nominal_span;
        end
        detected = detectPlateauMax(adc_counts,minSamples=options.autoMaxMinSamples,...
                                    tolerance=options.autoMaxTolerance);

        if isnan(detected) || (detected - min_count) < options.autoMaxMinSpan * nominal_span
            warning('voltage2percent:AutoMaxFailed',...
                ['Detected max (',num2str(detected),' counts) is too close to the range floor (',...
                 num2str(min_count),' counts), kept nominal max (',num2str(max_count),' counts). ',...
                 'Was the clamp on during this session?']);
        else
            max_count = detected;
        end
    end

    % Map to percentage
    pct = ((adc_counts - min_count) / (max_count - min_count)) * 100;

    % Clamp to 0-100%
    pct = max(0, min(100, pct));
    pct(isnan(adc_counts)) = NaN; % min/max ignore NaN, would turn them into 100%
    usedRange = [min_count, max_count];
end
