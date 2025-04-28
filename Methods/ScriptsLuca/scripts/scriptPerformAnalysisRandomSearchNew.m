function peakData = scriptPerformAnalysisRandomSearchNew(averageData, baselineSD, halfISI, responseSign, thresholdFactor)
    
    warning('off', 'signal:findpeaks:largeMinPeakHeight');
    
    peakWindowDuration = 500;
    peakWindowStart = halfISI;
    peakWindowEnd = min(halfISI + peakWindowDuration,length(averageData));
    
    controlWindowDuration = peakWindowDuration;
    controlWindowStart = max(1, halfISI - 1 - controlWindowDuration);
    controlWindowEnd = halfISI - 1;
    
    thresholdPeak =  thresholdFactor*baselineSD;
    
    peakData = [];
    peakData.thresholdValue = thresholdPeak;

    % Peak features
    [pulsePeaks, pulsePeakIndices] = findpeaks(responseSign*averageData(peakWindowStart:peakWindowEnd));
    [pulsePeakHeight, maxPulsePeakIndex] = max(pulsePeaks);
    pulsePeakTime = pulsePeakIndices(maxPulsePeakIndex);
    
    peakData.heightPulsePeak = responseSign*pulsePeakHeight;
    peakData.timePulsePeak = pulsePeakTime;
    peakData.areaPulse = trapz(peakWindowStart:peakWindowEnd, averageData(peakWindowStart:peakWindowEnd))/10000;  
    
    if pulsePeakHeight > thresholdPeak
        
        peakData.isPulsePeakAboveThreshold = 1;
        
    else
        
        peakData.isPulsePeakAboveThreshold = 0;

    end
    
    % Control trace features
    [controlPeaks, ~] = findpeaks(responseSign*averageData(controlWindowStart:controlWindowEnd));
    [controlPeakHeight, ~] = max(controlPeaks);
    
    peakData.heightControlPeak = responseSign*controlPeakHeight;
    peakData.areaControl = trapz(controlWindowStart:controlWindowEnd, averageData(controlWindowStart:controlWindowEnd))/10000;
    
    if controlPeakHeight > thresholdPeak
        
        peakData.isControlPeakAboveThreshold = 1;
        
    else
        
        peakData.isControlPeakAboveThreshold = 0;
        
    end
        
end