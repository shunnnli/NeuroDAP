function peakData = extractPeakData(data, responseThreshold, responseSign)
    
    warning('off', 'signal:findpeaks:largeMinPeakHeight');
    
    if responseSign == 1

        [peaks, ~] = findpeaks(data);
        peakAmplitude = max(peaks);

    elseif responseSign == -1

        [peaks, ~] = findpeaks(-data);
        peakAmplitude = -max(peaks);
    
    end
    
    if peakAmplitude > responseThreshold
        
        isAboveThreshold = 1;
        
    else
        
        isAboveThreshold = 0;       
        
    end
    
    peakData = [peakAmplitude,isAboveThreshold];
        
end