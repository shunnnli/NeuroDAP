function activeBoxes = defineAnalysisRandomSearch(dataWave0, state)

    % Reorganize pulse parameters and populate substructure

    pulsesFields = ["pulseString_ao0","pulseString_ao1","pulseString_ao2","pulseString_ao3", ...
        "pulseString_do0","pulseString_do1","pulseString_do2","pulseString_do3", ...
        "pulseString_aux0","pulseString_aux1","pulseString_aux2","pulseString_aux3", ...
        "pulseString_aux4","pulseString_aux5","pulseString_aux6","pulseString_aux7"];

    for iPulseField = 1:length(pulsesFields)

        pulseString = state.phys.internal.(pulsesFields(iPulseField));
        iPulseParamStart = [1 strfind(pulseString,';')]+1;
        iPulseParamEntryStart = [strfind(pulseString,'=')]+1;

        for i = 1:(length(iPulseParamStart)-1)

            iPulseParamName = pulseString(iPulseParamStart(i):(iPulseParamEntryStart(i)-2));
            iPulseParamValue = pulseString(iPulseParamEntryStart(i):(iPulseParamStart(i+1)-2));
            state.phys.internal.pulses.(pulsesFields(iPulseField)).(iPulseParamName) = iPulseParamValue;

        end 

    end

    activeChannels = regexp(state.phys.internal.lastLinesUsed, '''(\w+)''', 'tokens');

    optoParameters = [];
    optoParameters.activeChannels = activeChannels;

    if any(cellfun(@(x) contains(x, 'ao1'), activeChannels)) && any(cellfun(@(x) ~contains(x, 'ao2'), activeChannels))

        nPulsesBlue = str2double(state.phys.internal.pulses.pulseString_ao1.numPulses);
        delayFirstPulseBlue = str2double(state.phys.internal.pulses.pulseString_ao1.delay);
        isiBlue = str2double(state.phys.internal.pulses.pulseString_ao1.isi);
        pulseWidthBlue = str2double(state.phys.internal.pulses.pulseString_ao1.pulseWidth)*10;
        amplitudeBlue = str2double(state.phys.internal.pulses.pulseString_ao1.amplitude);
        timeArray = linspace(0,size(dataWave0.data,2)-1,size(dataWave0.data,2));
        pulsesTimeArray = (delayFirstPulseBlue + linspace(0,(nPulsesBlue - 1)*isiBlue,nPulsesBlue))*10;
        pulsesWidth = pulseWidthBlue;
        pulsesTimeDifference = nan;
        whichPulseFirst = nan;

        optoParameters.nPulsesBlue = nPulsesBlue;
        optoParameters.pulseWidthBlue = pulseWidthBlue;
        optoParameters.amplitudeBlue = amplitudeBlue;
        %optoParameters.functionNameBlue = functionNameBlue;
        %optoParameters.boxDefinitionBlue = boxDefinitionBlue;

        optoParameters.nPulsesRed = nan;
        optoParameters.pulseWidthRed = nan;
        optoParameters.amplitudeRed = nan;
        %optoParameters.functionNameRed = nan;
        %optoParameters.boxDefinitionRed = nan;

        optoParameters.whichPulseFirst = whichPulseFirst;
        optoParameters.pulsesTimeDifference = pulsesTimeDifference;
        
        timeBeforePulse = 1000;
        timeBaselineAverage = 500;
        timeAfterPulse = 3000;

    elseif any(cellfun(@(x) contains(x, 'ao2'), activeChannels)) && any(cellfun(@(x) ~contains(x, 'ao1'), activeChannels))

        nPulsesRed = str2double(state.phys.internal.pulses.pulseString_ao2.numPulses);
        delayFirstPulseRed = str2double(state.phys.internal.pulses.pulseString_ao2.delay);
        isiRed = str2double(state.phys.internal.pulses.pulseString_ao2.isi);
        pulseWidthRed = str2double(state.phys.internal.pulses.pulseString_ao2.pulseWidth)*10;
        amplitudeRed = str2double(state.phys.internal.pulses.pulseString_ao2.amplitude);                    
        timeArray = linspace(0,size(dataWave0.data,2)-1,size(dataWave0.data,2));
        pulsesTimeArray = (delayFirstPulseRed + linspace(0,(nPulsesRed - 1)*isiRed,nPulsesRed))*10;
        pulsesWidth = pulseWidthRed;
        pulsesTimeDifference = nan;
        whichPulseFirst = nan;

        optoParameters.nPulsesBlue = nan;
        optoParameters.pulseWidthBlue = nan;
        optoParameters.amplitudeBlue = nan;
        %optoParameters.functionNameBlue = nan;
        %optoParameters.boxDefinitionBlue = nan;

        optoParameters.nPulsesRed = nPulsesRed;
        optoParameters.pulseWidthRed = pulseWidthRed;
        optoParameters.amplitudeRed = amplitudeRed;
        %optoParameters.functionNameRed = functionNameRed;
        %optoParameters.boxDefinitionRed = boxDefinitionRed;

        optoParameters.whichPulseFirst = whichPulseFirst;
        optoParameters.pulsesTimeDifference = pulsesTimeDifference;

    elseif any(cellfun(@(x) contains(x, 'ao1'), activeChannels)) && any(cellfun(@(x) contains(x, 'ao2'), activeChannels))

        nPulsesBlue = str2double(state.phys.internal.pulses.pulseString_ao1.numPulses);
        delayFirstPulseBlue = str2double(state.phys.internal.pulses.pulseString_ao1.delay);
        isiBlue = str2double(state.phys.internal.pulses.pulseString_ao1.isi);
        pulseWidthBlue = str2double(state.phys.internal.pulses.pulseString_ao1.pulseWidth)*10;
        amplitudeBlue = str2double(state.phys.internal.pulses.pulseString_ao1.amplitude);
        timeArray = linspace(0,size(dataWave0.data,2)-1,size(dataWave0.data,2));
        pulsesTimeArrayBlue = (delayFirstPulseBlue + linspace(0,(nPulsesBlue - 1)*isiBlue,nPulsesBlue))*10;

        nPulsesRed = str2double(state.phys.internal.pulses.pulseString_ao2.numPulses);
        delayFirstPulseRed = str2double(state.phys.internal.pulses.pulseString_ao2.delay);
        isiRed = str2double(state.phys.internal.pulses.pulseString_ao2.isi);
        pulseWidthRed = str2double(state.phys.internal.pulses.pulseString_ao2.pulseWidth)*10;
        amplitudeRed = str2double(state.phys.internal.pulses.pulseString_ao2.amplitude);
        pulsesTimeArrayRed = (delayFirstPulseRed + linspace(0,(nPulsesRed - 1)*isiRed,nPulsesRed))*10;

        % Align signal with respect to the second pulse

        if pulsesTimeArrayBlue(1) < pulsesTimeArrayRed(1)

            pulsesTimeArray = pulsesTimeArrayRed;
            pulsesTimeDifference = pulsesTimeArrayRed(1) - pulsesTimeArrayBlue(1);
            whichPulseFirst = 'Blue';

        elseif pulsesTimeArrayBlue(1) >= pulsesTimeArrayRed(1)

            pulsesTimeArray = pulsesTimeArrayBlue;
            pulsesTimeDifference = pulsesTimeArrayBlue(1) - pulsesTimeArrayRed(1);
            whichPulseFirst = 'Red';

        end

        optoParameters.nPulsesBlue = nPulsesBlue;
        optoParameters.pulseWidthBlue = pulseWidthBlue;
        optoParameters.amplitudeBlue = amplitudeBlue;
        %optoParameters.functionNameBlue = functionNameBlue;
        %optoParameters.boxDefinitionBlue = boxDefinitionBlue;

        optoParameters.nPulsesRed = nPulsesRed;
        optoParameters.pulseWidthRed = pulseWidthRed;
        optoParameters.amplitudeRed = amplitudeRed;
        %optoParameters.functionNameRed = functionNameRed;
        %optoParameters.boxDefinitionRed = boxDefinitionRed;

        optoParameters.whichPulseFirst = whichPulseFirst;
        optoParameters.pulsesTimeDifference = pulsesTimeDifference;

    end

    segments = cell(1, length(pulsesTimeArray));
    activeBoxes = length(pulsesTimeArray));
    halfISI = (pulsesTimeArray(2) - pulsesTimeArray(1))/2;

    for iPulse = 1:length(pulsesTimeArray)

        pulseIdx = round(pulsesTimeArray(iPulse));
        startIdx = pulseIdx - halfISI;
        endIdx = pulseIdx + halfISI - 1;
        startIdx = max(1, startIdx);
        endIdx = min(length(dataWave0.data), endIdx);
        baselineAverage = mean(dataWave0.data(startIdx:pulseIdx));
        baselineSD = std(dataWave0.data(startIdx:pulseIdx));
        segments{iPulse} = dataWave0.data(startIdx:endIdx) - baselineAverage;
        segments{iPulse} = preprocessSignalVC(segments{iPulse});

        boxResponse = performAnalysisRandomSearch(segments{iPulse}, baselineSD, halfISI);
        activeBoxes(iPulse) = boxResponse;
        
    end

end
