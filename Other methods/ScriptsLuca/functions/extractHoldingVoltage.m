function holdingVoltage = extractHoldingVoltage(holdingVoltageDataFile)

    holdingVoltageData = holdingVoltageDataFile.data(10001:end)*100;
    actualMeanVoltage = mean(holdingVoltageData);
    
    if actualMeanVoltage >= -75 && actualMeanVoltage < -65; holdingVoltage = -70;
    elseif actualMeanVoltage >= -40 && actualMeanVoltage < -30; holdingVoltage = -35; 
    elseif actualMeanVoltage >= -5 && actualMeanVoltage < 5; holdingVoltage = 0; 
    elseif actualMeanVoltage >= 5 && actualMeanVoltage < 15; holdingVoltage = 10;
    else holdingVoltage = round(actualMeanVoltage/5)*5; % rounding in steps of 5
    end

end