function boxResponse = performAnalysisRandomSearch(data, baselineSD, halfISI)
    
    warning('off', 'signal:findpeaks:largeMinPeakHeight');

    preStimulusAverage = mean(data((halfISI-300):(halfISI)));
    preStimulusSD = baselineSD;
    peakWindowDuration = 500;
    peakWindowStart = halfISI;
    peakWindowEnd = halfISI + peakWindowDuration;
    thresholdPeak =  5*preStimulusSD;
    averageData = data - preStimulusAverage;
    [peaks, ~] = findpeaks(averageData(peakWindowStart:peakWindowEnd), 'MinPeakHeight', thresholdPeak);
    
    if ~isempty(peaks)
        
        boxResponse = 1;
        
    else
        
        boxResponse = 0;
        
    end
        
end