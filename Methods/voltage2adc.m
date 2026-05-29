function adc_counts = voltage2adc(voltage, options)

arguments
    voltage double
    options.Vref double = 5
    options.ADC_resolution double = 4095
end
    % Convert voltage to ADC counts
    adc_counts = (voltage / options.Vref) * options.ADC_resolution;
end