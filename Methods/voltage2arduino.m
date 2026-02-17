function output = voltage2arduino(voltage,options)

    % Convert raw voltage signal to arduino signal
    arguments
        voltage double
        options.Vref double = 5
        options.ADCresolution double = 1023
    end

    output = (voltage / options.Vref) * options.ADCresolution;
end