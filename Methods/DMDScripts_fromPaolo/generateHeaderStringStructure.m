function headerStringStructure = generateHeaderStringStructure(runPath) 
% Split string, create and populate structure according to headerString

    cyclePositionFile = load(runPath);
    [~, baseName] = fileparts(runPath);
    string = cyclePositionFile.(baseName).UserData.headerString;
    iStateStart = strfind(string,'state');
    headerStringStructure = [];

    for i = 1:length(iStateStart)

        subString = string(iStateStart(i):iStateEnd(i,iStateStart,string));
        iStateEntryStart = strfind(subString,'=');
        iStateEntryStart = iStateEntryStart(1);
        iStateName = subString(1:(iStateEntryStart-1));
        subStringEntry = subString((iStateEntryStart+1):length(subString));

        structFields  = strsplit(iStateName,'.');
        headerStringStructure = setfield(headerStringStructure, structFields{1:end}, subStringEntry);

    end

    % Reorganize pulse parameters and populate substructure

    pulsesFields = ["pulseString_ao0","pulseString_ao1","pulseString_ao2","pulseString_ao3", ...
        "pulseString_do0","pulseString_do1","pulseString_do2","pulseString_do3", ...
        "pulseString_aux0","pulseString_aux1","pulseString_aux2","pulseString_aux3", ...
        "pulseString_aux4","pulseString_aux5","pulseString_aux6","pulseString_aux7"];

    for iPulseField = 1:length(pulsesFields)

        pulseString = headerStringStructure.state.phys.internal.(pulsesFields(iPulseField));
        iPulseParamStart = [1 strfind(pulseString,';')]+1;
        iPulseParamEntryStart = [strfind(pulseString,'=')]+1;

        for i = 1:(length(iPulseParamStart)-1)

            iPulseParamName = pulseString(iPulseParamStart(i):(iPulseParamEntryStart(i)-2));
            iPulseParamValue = pulseString(iPulseParamEntryStart(i):(iPulseParamStart(i+1)-2));
            headerStringStructure.state.phys.internal.pulses.(pulsesFields(iPulseField)).(iPulseParamName) = iPulseParamValue;

        end 

    end
    
    activeChannels = regexp(headerStringStructure.state.phys.internal.lastLinesUsed, '''(\w+)''', 'tokens');
                
    if any(cellfun(@(x) contains(x, 'ao1'), activeChannels)) && any(cellfun(@(x) ~contains(x, 'ao2'), activeChannels))
                 
        headerStringStructure.dmdBlue = 'ON';
        headerStringStructure.dmdRed = 'OFF';       
 
    elseif any(cellfun(@(x) contains(x, 'ao2'), activeChannels)) && any(cellfun(@(x) ~contains(x, 'ao1'), activeChannels))

        headerStringStructure.dmdBlue = 'OFF';
        headerStringStructure.dmdRed = 'ON';
        
    elseif any(cellfun(@(x) contains(x, 'ao1'), activeChannels)) && any(cellfun(@(x) contains(x, 'ao2'), activeChannels))

        headerStringStructure.dmdBlue = 'ON';
        headerStringStructure.dmdRed = 'ON';
        
    elseif any(cellfun(@(x) ~contains(x, 'ao1'), activeChannels)) && any(cellfun(@(x) ~contains(x, 'ao2'), activeChannels))
        
        headerStringStructure.dmdBlue = 'OFF';
        headerStringStructure.dmdRed = 'OFF';
        
    end

end

function iStateEnd = iStateEnd(i,iStateStart,string)
    
    if i == length(iStateStart)
        iStateEnd = length(string);
    else
        iStateEnd = iStateStart(i+1)-1;
    end
    
end