function APdistribution = computeAPdistribution(DataStruct, optoParameters)
    
    APdistribution = [];

    currentPulseStart = optoParameters.delayCurrentPulse;
    currentPulseEnd = currentPulseStart + optoParameters.widthCurrentPulse - currentPulseStart;
    lightPulseStart = optoParameters.delayLightPulse - currentPulseStart;
    lightPulseEnd = optoParameters.delayLightPulse + optoParameters.widthLightPulse - currentPulseStart;
    
    if ~isstruct(DataStruct.phAnalysis.pulseAP{1,1})
        
        APdistribution.nAPtotal = 0;
        APdistribution.nAPpreLight = nan;
        APdistribution.nAPduringLight = nan;
        APdistribution.nAPpostLight = nan;
        APdistribution.rateAPpreLight = nan;
        APdistribution.rateAPduringLight = nan;
        APdistribution.rateAPpostLight = nan;
        APdistribution.delayFirstAPpostLight = nan;
        return;
        
    end

    timesAP = DataStruct.phAnalysis.pulseAP{1,1}.AP_peak_time*10;
    indexFirstAPpostLight = find(timesAP > lightPulseEnd, 1);
    timeFirstAPpostLight = timesAP(indexFirstAPpostLight);
    delayFirstAPpostLight = (timeFirstAPpostLight - lightPulseEnd)/10;
    
    if isempty(timeFirstAPpostLight)
        
        timeFirstAPpostLight = lightPulseEnd;
        delayFirstAPpostLight = nan;
        
    end
    
    nAPpreLight = sum(timesAP >= 0 & timesAP < lightPulseStart);
    nAPduringLight = sum(timesAP >= lightPulseStart & timesAP < lightPulseEnd);
    nAPpostLight = sum(timesAP >= lightPulseEnd & timesAP < currentPulseEnd);
    
    rateAPpreLight = nAPpreLight / (lightPulseStart/10000);
    rateAPduringLight = nAPduringLight / ((lightPulseEnd - lightPulseStart)/10000);
    rateAPpostLight = nAPpostLight / ((currentPulseEnd - timeFirstAPpostLight)/10000);
    
    APdistribution.nAPtotal = size(timesAP,2);
    APdistribution.nAPpreLight = nAPpreLight;
    APdistribution.nAPduringLight = nAPduringLight;
    APdistribution.nAPpostLight = nAPpostLight;
    APdistribution.rateAPpreLight = rateAPpreLight;
    APdistribution.rateAPduringLight = rateAPduringLight;
    APdistribution.rateAPpostLight = rateAPpostLight;
    APdistribution.delayFirstAPpostLight = delayFirstAPpostLight;

end