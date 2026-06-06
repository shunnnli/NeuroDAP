function peakData = scriptPerformAnalysisVCtwoSigns(data, baselineSD, pulseTime, responseSign)
    
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
    
    averageDataPeakWindow = averageData(peakWindowStart:peakWindowEnd);
    averageDataControlWindow = averageData(controlWindowStart:controlWindowEnd);
    
    peakAverageWindow = 30;
    thresholdPulsePeak =  5*baselineSD;
    
    if responseSign ~= 2
        
        [pulsePeaks, pulsePeakIndices] = findpeaks(responseSign*averageDataPeakWindow);
        [~, maxPulsePeakIndex] = max(pulsePeaks);
        pulsePeakTime = pulsePeakIndices(maxPulsePeakIndex);
        pulsePeakHeight = mean(responseSign*averageDataPeakWindow(max(1,(pulsePeakTime-peakAverageWindow/2)):min((pulsePeakTime+peakAverageWindow/2-1),peakWindowDuration)));

        peakData.heightPulsePeak = responseSign*pulsePeakHeight;
        peakData.timePulsePeak = pulsePeakTime/10;

        if pulsePeakHeight > thresholdPulsePeak

            peakData.isPulsePeakAboveThreshold = 1;

        else

            peakData.isPulsePeakAboveThreshold = 0;

        end

        peakData.areaPulse = trapz(peakWindowStart:peakWindowEnd, averageDataPeakWindow)/10000;    

        [controlPeaks, controlPeakIndices] = findpeaks(responseSign*averageDataControlWindow);
        [~, maxControlPeakIndex] = max(controlPeaks);
        controlPeakTime = controlPeakIndices(maxControlPeakIndex);
        controlPeakHeight = mean(responseSign*averageDataControlWindow(max(1,(controlPeakTime-peakAverageWindow/2)):min((controlPeakTime+peakAverageWindow/2-1),peakWindowDuration)));

        peakData.heightControlPeak = responseSign*controlPeakHeight;
        peakData.areaControl = trapz(controlWindowStart:controlWindowEnd, averageDataControlWindow)/10000;

    else
        
        positiveResponseSign = 1;
        negativeResponseSign = -1;
        
        % Positive response
        [pulsePeaks, pulsePeakIndices] = findpeaks(positiveResponseSign*averageDataPeakWindow);      
        [~, maxPulsePeakIndex] = max(pulsePeaks);
        pulsePeakTime = pulsePeakIndices(maxPulsePeakIndex);
        pulsePeakHeight = mean(positiveResponseSign*averageDataPeakWindow(max(1,(pulsePeakTime-peakAverageWindow/2)):min((pulsePeakTime+peakAverageWindow/2-1),peakWindowDuration)));

        heightPulsePeakPositive = positiveResponseSign*pulsePeakHeight;
        timePulsePeakPositive = pulsePeakTime;

        if pulsePeakHeight > thresholdPulsePeak

            isPulsePeakAboveThresholdPositive = 1;

        else

            isPulsePeakAboveThresholdPositive = 0;

        end

        %peakData.areaPulse = trapz(peakWindowStart:peakWindowEnd, averageData(peakWindowStart:peakWindowEnd))/10000;    

        [controlPeaks, controlPeakIndices] = findpeaks(positiveResponseSign*averageDataControlWindow);
        [~, maxControlPeakIndex] = max(controlPeaks);
        controlPeakTime = controlPeakIndices(maxControlPeakIndex);
        controlPeakHeight = mean(positiveResponseSign*averageDataControlWindow(max(1,(controlPeakTime-peakAverageWindow/2)):min((controlPeakTime+peakAverageWindow/2-1),peakWindowDuration)));

        heightControlPeakPositive = positiveResponseSign*controlPeakHeight;
        %peakData.areaControl = trapz(controlWindowStart:controlWindowEnd, averageData(controlWindowStart:controlWindowEnd))/10000;

        % Negative response
        [pulsePeaks, pulsePeakIndices] = findpeaks(negativeResponseSign*averageDataPeakWindow);        
        [~, maxPulsePeakIndex] = max(pulsePeaks);
        pulsePeakTime = pulsePeakIndices(maxPulsePeakIndex);
        pulsePeakHeight = mean(negativeResponseSign*averageDataPeakWindow(max(1,(pulsePeakTime-peakAverageWindow/2)):min((pulsePeakTime+peakAverageWindow/2-1),peakWindowDuration)));

        heightPulsePeakNegative = negativeResponseSign*pulsePeakHeight;
        timePulsePeakNegative = pulsePeakTime;

        if pulsePeakHeight > thresholdPulsePeak

            isPulsePeakAboveThresholdNegative = 1;

        else

            isPulsePeakAboveThresholdNegative = 0;

        end

        %peakData.areaPulse = trapz(peakWindowStart:peakWindowEnd, averageData(peakWindowStart:peakWindowEnd))/10000;    

        [controlPeaks, controlPeakIndices] = findpeaks(negativeResponseSign*averageDataControlWindow);
        [~, maxControlPeakIndex] = max(controlPeaks);
        controlPeakTime = controlPeakIndices(maxControlPeakIndex);
        controlPeakHeight = mean(negativeResponseSign*averageDataControlWindow(max(1,(controlPeakTime-peakAverageWindow/2)):min((controlPeakTime+peakAverageWindow/2-1),peakWindowDuration)));

        heightControlPeakNegative = negativeResponseSign*controlPeakHeight;
        %peakData.areaControl = trapz(controlWindowStart:controlWindowEnd, averageData(controlWindowStart:controlWindowEnd))/10000;

        % Areas
        dataAbove = averageData .* (averageData >= 0);
        crossingTimes = find(diff(dataAbove>0));

        timeFirstPeak = min(timePulsePeakNegative,timePulsePeakPositive);
        crossingTimes = crossingTimes(crossingTimes>=timeFirstPeak);
        crossingTime = crossingTimes(1);
        
        negativeWindowStart = pulseTime;
        negativeWindowEnd = crossingTime;
        
        positiveWindowStart = crossingTime + 1;
        positiveWindowEnd = crossingTime + peakWindowDuration + 1;
        
        areaPulsePositive = trapz(positiveWindowStart:positiveWindowEnd, averageData(positiveWindowStart:positiveWindowEnd))/10000;    
        areaControlPositive = trapz(controlWindowStart:controlWindowEnd, averageData(controlWindowStart:controlWindowEnd))/10000;
        areaPulseNegative = trapz(negativeWindowStart:negativeWindowEnd, averageData(negativeWindowStart:negativeWindowEnd))/10000;    
        areaControlNegative = trapz(controlWindowStart:controlWindowEnd, averageData(controlWindowStart:controlWindowEnd))/10000;
       
        % Positive and negative 
        peakData.heightPulsePeak = [heightPulsePeakPositive, heightPulsePeakNegative];
        peakData.timePulsePeak = [timePulsePeakPositive/10, timePulsePeakNegative/10];
        peakData.isPulsePeakAboveThreshold = [isPulsePeakAboveThresholdPositive, isPulsePeakAboveThresholdNegative];     
        peakData.areaPulse = [areaPulsePositive, areaPulseNegative];       
        peakData.heightControlPeak = [heightControlPeakPositive, heightControlPeakNegative];
        peakData.areaControl = [areaControlPositive, areaControlNegative]; 
        
    end
    
    peakData.trace = averageData;

end