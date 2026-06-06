function peakData = scriptPerformAnalysisVC(data, baselineSD, pulseTime, responseSign)
    
    warningState = warning('off', 'signal:findpeaks:largeMinPeakHeight');
    
    preStimulusAverage = mean(data(max(1,(pulseTime-300)):(pulseTime)));
    averageData = data - preStimulusAverage;

    peakData = [];

    peakWindowDuration = 500;
    peakWindowStart = pulseTime;
    peakWindowEnd = min(pulseTime + peakWindowDuration, size(averageData,2));
    
    controlWindowDuration = peakWindowDuration;
    controlWindowStart = max(1, pulseTime - controlWindowDuration);
    controlWindowEnd = pulseTime;
    
    thresholdPulsePeak =  5*baselineSD;
    
    [pulsePeaks, pulsePeakIndices] = findpeaks(responseSign*averageData(peakWindowStart:peakWindowEnd));
    [pulsePeakHeight, maxPulsePeakIndex] = max(pulsePeaks);
    pulsePeakTime = pulsePeakIndices(maxPulsePeakIndex); %+ peakWindowStart;
    
    peakData.heightPulsePeak = responseSign*pulsePeakHeight;
    peakData.timePulsePeak = pulsePeakTime;
    
    if pulsePeakHeight > thresholdPulsePeak
        
        peakData.isPulsePeakAboveThreshold = 1;
        
    else
        
        peakData.isPulsePeakAboveThreshold = 0;
        
    end

    peakData.areaPulse = trapz(peakWindowStart:peakWindowEnd, averageData(peakWindowStart:peakWindowEnd))/10000;    

    [controlPeaks, ~] = findpeaks(responseSign*averageData(controlWindowStart:controlWindowEnd)); %, 'MinPeakHeight', thresholdPulsePeak
    [controlPeakHeight, ~] = max(controlPeaks);
    
    peakData.heightControlPeak = responseSign*controlPeakHeight;
    peakData.areaControl = trapz(controlWindowStart:controlWindowEnd, averageData(controlWindowStart:controlWindowEnd))/10000;
    
    peakData.trace = data;

end