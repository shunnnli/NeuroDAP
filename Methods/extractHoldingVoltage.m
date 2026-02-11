function holdingVoltage = extractHoldingVoltage(holdingVoltageDataFile)

    holdingVoltageData = holdingVoltageDataFile.data * 100;
    nbins = 200;
    [counts, edges] = histcounts(holdingVoltageData, nbins);
    [~, i] = max(counts);
    holdingVoltage = mean(edges(i:i+1));   % bin center of the most-populated bin
    
%     if actualMeanVoltage >= -75 && actualMeanVoltage < -65; holdingVoltage = -70;
%     elseif actualMeanVoltage >= -40 && actualMeanVoltage < -30; holdingVoltage = -35; 
%     elseif actualMeanVoltage >= -5 && actualMeanVoltage < 5; holdingVoltage = 0; 
%     elseif actualMeanVoltage >= 5 && actualMeanVoltage < 15; holdingVoltage = 10;
%     else holdingVoltage = round(actualMeanVoltage/5)*5; % rounding in steps of 5
%     end

end