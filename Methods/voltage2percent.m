function pct = voltage2percent(voltage, range, options)

arguments
    voltage double
    range double
    options.Vref double = 5
    options.ADC_resolution double = 4095
end
    min_count = range(1);
    max_count = range(end);

    % Convert voltage to ADC counts
    adc_counts = (voltage / options.Vref) * options.ADC_resolution;

    % Map to percentage
    pct = ((adc_counts - min_count) / (max_count - min_count)) * 100;

    % Clamp to 0-100%
    pct = max(0, min(100, pct));
end