disp(['*-*-*-* Running: scriptCollectFeaturesVC *-*-*-*'])

 for iCell = startDir:(startDir+nCells-1)
     
     currentCell = cellList(iCell);
     
     for iSelectedCell = 1:numel(selectedCells)  

        if isequal(currentCell.name, selectedCells(iSelectedCell))
            
            disp(['****** Loaded ' currentCell.name ' ******'])
            cellPath = [expPath filesep cellList(iCell).name];
            processedDataPath = [cellPath filesep 'ProcessedData' filesep 'Data'];
            analysisData = [];
            [spanEpochs, spanCyclePositions] = findSpanEpochsCyclePositions(processedDataPath);
            
            columnNames = {'epoch' 'cyclePosition', 'heightPulsePeak', 'timePulsePeak', 'isPulsePeakAboveThreshold', 'areaPulse', 'heightControlPeak', 'areaControl', 'trace', 'optoParameters'};
            dataTypes = {'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell'};
            analysisTable = table('Size', [0, numel(columnNames)], 'VariableNames', columnNames, 'VariableTypes', dataTypes);
            
            for i = 1:numel(columnNames)
                
                analysisTable.(columnNames{i}) = repmat({[]}, 0, 1);
                
                if strcmp(dataTypes{i}, 'double')
                    
                    analysisTable.(columnNames{i}) = NaN(0, 1);
                
                elseif strcmp(dataTypes{i}, 'logical')
                    
                    analysisTable.(columnNames{i}) = false(0, 1);
                
                end
                
            end
            
            for iCyclePosition = spanCyclePositions
                
                nameCyclePosition = strcat(['cyclePosition', num2str(iCyclePosition)]);
                analysisData.(nameCyclePosition) = [];
                
                for iEpoch = spanEpochs
                    
                    nameEpoch = strcat(['epoch', num2str(iEpoch)]);

                    runPath = [cellPath filesep 'ProcessedData' filesep 'Data' filesep 'Epoch' num2str(iEpoch) '_cyclePosition' num2str(iCyclePosition) '.mat'];

                    if ~isfile(runPath)
                        continue;
                    end

                    disp(['Analyzing: Epoch ' num2str(iEpoch) ' - cyclePosition ' num2str(iCyclePosition)])
                    runData = load(runPath);
                    DataStruct = runData.DataStruct;
                    runInfo = [iEpoch, iCyclePosition];

                    data = DataStruct.data;
                    activeChannels = regexp(DataStruct.state.phys.internal.lastLinesUsed, '''(\w+)''', 'tokens');
                    
                    optoParameters = [];
                    optoParameters.activeChannels = activeChannels;
                    optoParameters.holdingVoltage = DataStruct.holdingVoltage;
 
                    if (any(cellfun(@(x) contains(x, 'ao1'), activeChannels)) && any(cellfun(@(x) contains(x, 'ao2'), activeChannels)))

                        nPulsesBlue = str2double(DataStruct.state.phys.internal.pulses.pulseString_ao1.numPulses);
                        delayFirstPulseBlue = str2double(DataStruct.state.phys.internal.pulses.pulseString_ao1.delay)*10;
                        isiBlue = str2double(DataStruct.state.phys.internal.pulses.pulseString_ao1.isi)*10;
                        pulseWidthBlue = str2double(DataStruct.state.phys.internal.pulses.pulseString_ao1.pulseWidth)*10;
                        amplitudeBlue = str2double(DataStruct.state.phys.internal.pulses.pulseString_ao1.amplitude);
                        timeArray = linspace(0,size(DataStruct.data,2)-1,size(DataStruct.data,2));
                        pulsesTimeArrayBlue = (delayFirstPulseBlue + linspace(0,(nPulsesBlue - 1)*isiBlue,nPulsesBlue));

                        nPulsesRed = str2double(DataStruct.state.phys.internal.pulses.pulseString_ao2.numPulses);
                        delayFirstPulseRed = str2double(DataStruct.state.phys.internal.pulses.pulseString_ao2.delay)*10;
                        isiRed = str2double(DataStruct.state.phys.internal.pulses.pulseString_ao2.isi)*10;
                        pulseWidthRed = str2double(DataStruct.state.phys.internal.pulses.pulseString_ao2.pulseWidth)*10;
                        amplitudeRed = str2double(DataStruct.state.phys.internal.pulses.pulseString_ao2.amplitude);
                        pulsesTimeArrayRed = (delayFirstPulseRed + linspace(0,(nPulsesRed - 1)*isiRed,nPulsesRed));

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
                        
                        functionNameBlue = DataStruct.state.cycle.functionName;
                        functionNameBlue = strrep(functionNameBlue, '''', '');
                        functionNameBlue = strrep(functionNameBlue, ';', '');
                        functionNameBlue = strtrim(functionNameBlue);
                        
                        functionNameRed = DataStruct.state.cycle.functionName;
                        functionNameRed = strrep(functionNameRed, '''', '');
                        functionNameRed = strrep(functionNameRed, ';', '');
                        functionNameRed = strtrim(functionNameRed);

                        optoParameters.nPulsesBlue = nPulsesBlue;
                        optoParameters.pulseWidthBlue = pulseWidthBlue;
                        optoParameters.amplitudeBlue = amplitudeBlue;
                        optoParameters.delayPulseBlue = delayFirstPulseBlue;
                        optoParameters.functionNameBlue = functionNameBlue;

                        optoParameters.nPulsesRed = nPulsesRed;
                        optoParameters.pulseWidthRed = pulseWidthRed;
                        optoParameters.amplitudeRed = amplitudeRed;
                        optoParameters.delayPulseRed = delayFirstPulseRed;
                        optoParameters.functionNameRed = functionNameRed;

                        optoParameters.whichPulseFirst = whichPulseFirst;
                        optoParameters.pulsesTimeDifference = pulsesTimeDifference;
                          
                    elseif (any(cellfun(@(x) contains(x, 'ao1'), activeChannels)) && any(cellfun(@(x) ~contains(x, 'ao2'), activeChannels)))

                        nPulsesBlue = str2double(DataStruct.state.phys.internal.pulses.pulseString_ao1.numPulses);
                        delayFirstPulseBlue = str2double(DataStruct.state.phys.internal.pulses.pulseString_ao1.delay)*10;
                        isiBlue = str2double(DataStruct.state.phys.internal.pulses.pulseString_ao1.isi)*10;
                        pulseWidthBlue = str2double(DataStruct.state.phys.internal.pulses.pulseString_ao1.pulseWidth)*10;
                        amplitudeBlue = str2double(DataStruct.state.phys.internal.pulses.pulseString_ao1.amplitude);
                        timeArray = linspace(0,size(DataStruct.data,2)-1,size(DataStruct.data,2));
                        pulsesTimeArray = (delayFirstPulseBlue + linspace(0,(nPulsesBlue - 1)*isiBlue,nPulsesBlue));
                        pulsesWidth = pulseWidthBlue;
                        pulsesTimeDifference = nan;
                        whichPulseFirst = nan;
                        functionNameBlue = DataStruct.state.cycle.functionName;
                        functionNameBlue = strrep(functionNameBlue, '''', '');
                        functionNameBlue = strrep(functionNameBlue, ';', '');
                        functionNameBlue = strtrim(functionNameBlue);
                        
                        optoParameters.nPulsesBlue = nPulsesBlue;
                        optoParameters.pulseWidthBlue = pulseWidthBlue;
                        optoParameters.amplitudeBlue = amplitudeBlue;
                        optoParameters.delayPulseBlue = delayFirstPulseBlue;
                        optoParameters.functionNameBlue = functionNameBlue;

                        optoParameters.nPulsesRed = nan;
                        optoParameters.pulseWidthRed = nan;
                        optoParameters.amplitudeRed = nan;
                        optoParameters.delayPulseRed = nan;
                        optoParameters.functionNameRed = nan;

                        optoParameters.whichPulseFirst = whichPulseFirst;
                        optoParameters.pulsesTimeDifference = pulsesTimeDifference;                   

                    elseif (any(cellfun(@(x) contains(x, 'ao2'), activeChannels)) && any(cellfun(@(x) ~contains(x, 'ao1'), activeChannels)))

                        nPulsesRed = str2double(DataStruct.state.phys.internal.pulses.pulseString_ao2.numPulses);
                        delayFirstPulseRed = str2double(DataStruct.state.phys.internal.pulses.pulseString_ao2.delay)*10;
                        isiRed = str2double(DataStruct.state.phys.internal.pulses.pulseString_ao2.isi)*10;
                        pulseWidthRed = str2double(DataStruct.state.phys.internal.pulses.pulseString_ao2.pulseWidth)*10;
                        amplitudeRed = str2double(DataStruct.state.phys.internal.pulses.pulseString_ao2.amplitude);                    
                        timeArray = linspace(0,size(DataStruct.data,2)-1,size(DataStruct.data,2));
                        pulsesTimeArray = (delayFirstPulseRed + linspace(0,(nPulsesRed - 1)*isiRed,nPulsesRed));
                        pulsesWidth = pulseWidthRed;
                        pulsesTimeDifference = nan;
                        whichPulseFirst = nan;
                        functionNameRed = DataStruct.state.cycle.functionName;
                        functionNameRed = strrep(functionNameRed, '''', '');
                        functionNameRed = strrep(functionNameRed, ';', '');
                        functionNameRed = strtrim(functionNameRed);

                        optoParameters.nPulsesBlue = nan;
                        optoParameters.pulseWidthBlue = nan;
                        optoParameters.amplitudeBlue = nan;
                        optoParameters.delayPulseBlue = nan;
                        optoParameters.functionNameBlue = nan;

                        optoParameters.nPulsesRed = nPulsesRed;
                        optoParameters.pulseWidthRed = pulseWidthRed;
                        optoParameters.amplitudeRed = amplitudeRed;
                        optoParameters.delayPulseRed = delayFirstPulseRed;
                        optoParameters.functionNameRed = functionNameRed;

                        optoParameters.whichPulseFirst = whichPulseFirst;
                        optoParameters.pulsesTimeDifference = pulsesTimeDifference;

                    end
                    
                    if length(pulsesTimeArray) == 1
                                             
                        timeBeforePulse = 5000;
                        timeBaselineAverage = 500;
                        timeAfterPulse = 4999;%14999;
                        
                    else
                        
                        timeBeforePulse = (pulsesTimeArray(2) - pulsesTimeArray(1))/2 - 500;
                        timeBaselineAverage = 500;
                        timeAfterPulse = (pulsesTimeArray(2) - pulsesTimeArray(1))/2 - 501;
                        
                    end

                    segments = cell(1, length(pulsesTimeArray));

                    for iPulse = 1:length(pulsesTimeArray)

                        pulseIndex = round(pulsesTimeArray(iPulse));
                        startIdx = pulseIndex - timeBeforePulse;
                        endIdx = pulseIndex + timeAfterPulse;
                        endBaselineIdx = startIdx + timeBaselineAverage;
                        startIdx = max(1, startIdx);
                        endIdx = min(length(data), endIdx);
                        baselineAverage = mean([data(5000:10000),data((end-5000):end)]);
                        baselineSD = std([data(5000:10000),data((end-5000):end)]);
                        segments{iPulse} = data(startIdx:endIdx) - baselineAverage;
                        segments{iPulse} = preprocessSignalVC(segments{iPulse});

                    end

                    peakDataStruct = [];
                    peakDataStruct.heightPulsePeak = [];
                    peakDataStruct.timePulsePeak = [];
                    peakDataStruct.isPulsePeakAboveThreshold = [];
                    peakDataStruct.areaPulse = [];
                    peakDataStruct.heightControlPeak = [];
                    peakDataStruct.areaControl = [];   
                    peakDataStruct.trace = [];
                    
                    if DataStruct.holdingVoltage > -30

                        responseSign = 1;

                    else 

                        responseSign = -1;

                    end

                    for iPulse = 1:length(segments)

                        peakData = scriptPerformAnalysisVC(segments{iPulse}, baselineSD, timeBeforePulse, responseSign);

                        if ~isempty(peakData)

                            fields = fieldnames(peakData);

                            for iField = 1:numel(fields)

                                entryPeakData = peakData.(fields{iField});

                                if iPulse == 1

                                    peakDataStruct.(fields{iField}) = [entryPeakData];

                                else

                                    peakDataStruct.(fields{iField}) = [peakDataStruct.(fields{iField}); entryPeakData];

                                end

                            end

                        end

                    end
             
                    runDataStruct = [];  
                    
                    clear fields
                    
                    fieldsPeakDataStruct = fields(peakDataStruct);
                    
                    for iFieldPeakDataStruct = 1:numel(fieldsPeakDataStruct)
                    
                        runDataStruct = setfield(runDataStruct, fieldsPeakDataStruct{iFieldPeakDataStruct}, peakDataStruct.(fieldsPeakDataStruct{iFieldPeakDataStruct}));
                        
                    end
                    
                    runDataStruct.optoParameters = optoParameters;

                    analysisData.(nameCyclePosition).(nameEpoch) = runDataStruct;

                    runDataStructCell = struct2cell(runDataStruct);
                    dataToAdd = [{iEpoch, iCyclePosition}, runDataStructCell'];
                    analysisTable = [analysisTable; dataToAdd];

                end
            
                analyzedDataFolderName = 'Analysis';
                analyzedDataFolderPath = [cellPath filesep analyzedDataFolderName];

                if exist(analyzedDataFolderPath, 'dir') == 0

                    mkdir(analyzedDataFolderPath);
                    addpath(genpath(analyzedDataFolderPath));

                end

            end

            saveResultsName = 'AnalysisStruct.mat';
            saveResultsName = strrep(saveResultsName, ' ', '');
            saveResultsPath = fullfile(analyzedDataFolderPath, saveResultsName);
            save(saveResultsPath, 'analysisData');

            analysisTable = findIdenticalRunsVC(analysisTable);
            saveTableName = 'AnalysisTable.mat';
            saveTableName = strrep(saveTableName, ' ', '');
            saveTablePath = fullfile(analyzedDataFolderPath, saveTableName);
            save(saveTablePath, 'analysisTable');

            disp(['****** Results saved successfully ******']);

            analyzedMouseDataFolderName = 'mouseAnalysis';
            analyzedMouseDataFolderPath = [expPath filesep analyzedMouseDataFolderName];

            if exist(analyzedMouseDataFolderPath, 'dir') == 0

                mkdir(analyzedMouseDataFolderPath);
                addpath(genpath(analyzedMouseDataFolderPath));

            end

            saveMouseTableName = 'MouseAnalysisTable.mat';
            saveMouseTableName = strrep(saveMouseTableName, ' ', '');
            saveMouseTablePath = fullfile(analyzedMouseDataFolderPath, saveMouseTableName);

            if exist(saveMouseTablePath, 'file') == 0
                
                cellNumber = regexp(currentCell.name, '\d+', 'match');
                cellNumber = str2double(cellNumber);
                repeatedCellName = repmat({cellNumber}, size(analysisTable, 1), 1);
                analysisTable = addvars(analysisTable, repeatedCellName, 'Before', 1, 'NewVariableNames', 'cellName');
                analysisTable = removevars(analysisTable, 'runGroup');
                mouseAnalysisTable = analysisTable;
                mouseAnalysisTable = findIdenticalRunsVC(mouseAnalysisTable);
                save(saveMouseTablePath, 'mouseAnalysisTable');

            else 

                existingMouseAnalysisTable = load(fullfile(analyzedMouseDataFolderPath, saveMouseTableName));
                existingMouseAnalysisTable = existingMouseAnalysisTable.mouseAnalysisTable;
                
                if any(strcmp(existingMouseAnalysisTable.Properties.VariableNames, 'runGroup'))

                    existingMouseAnalysisTable = removevars(existingMouseAnalysisTable, 'runGroup');
                
                end

                if any(strcmp(existingMouseAnalysisTable.cellName, currentCell.name))

                    continue;

                end

                cellNumber = regexp(currentCell.name, '\d+', 'match');
                cellNumber = str2double(cellNumber);
                repeatedCellName = repmat({cellNumber}, size(analysisTable, 1), 1);
                analysisTable = addvars(analysisTable, repeatedCellName, 'Before', 1, 'NewVariableNames', 'cellName');
                analysisTable = removevars(analysisTable, 'runGroup');
                
                mouseAnalysisTable = [existingMouseAnalysisTable; analysisTable];
                mouseAnalysisTable = findIdenticalRunsVC(mouseAnalysisTable);
                save(saveMouseTablePath, 'mouseAnalysisTable');

            end
            
        end
        
     end
     
 end

mouseAnalysisTable = load(saveMouseTablePath);
mouseAnalysisTable = mouseAnalysisTable.mouseAnalysisTable;
 
runGroups = unique(mouseAnalysisTable.runGroup(:));
nRunGroups = numel(runGroups);

columnNames = {'runGroup', 'sampleSize', 'activeChannels', 'holdingVoltage', 'nPulsesBlue', 'pulseWidthBlue', 'amplitudeBlue', 'delayPulseBlue', 'functionNameBlue', 'nPulsesRed', 'pulseWidthRed', 'amplitudeRed', 'delayPulseRed', 'functionNameRed', 'whichPulseFirst', 'pulsesTimeDifference'};
dataTypes = {'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell'};
runGroupsParametersTable = table('Size', [0, numel(columnNames)], 'VariableNames', columnNames, 'VariableTypes', dataTypes);

for i = 1:numel(columnNames)

    runGroupsParametersTable.(columnNames{i}) = repmat({[]}, 0, 1);

    if strcmp(dataTypes{i}, 'double')

        runGroupsParametersTable.(columnNames{i}) = NaN(0, 1);

    elseif strcmp(dataTypes{i}, 'logical')

        runGroupsParametersTable.(columnNames{i}) = false(0, 1);

    end

end


for iRunGroup = 1:nRunGroups
    
    rowsRunGroup = find(mouseAnalysisTable.runGroup(:) == iRunGroup);
    sampleSizeRunGroup = size(rowsRunGroup,1);
    runGroupAnalysisTable = mouseAnalysisTable(rowsRunGroup,:);
    
    optoParameters = runGroupAnalysisTable.optoParameters{1};
    optoParametersCell = struct2cell(optoParameters);
    
    activeChannels = mouseAnalysisTable.optoParameters{1,1}.activeChannels{1};
    
    if numel(runGroupAnalysisTable.optoParameters{1,1}.activeChannels) == 2
        
        secondActiveChannel = runGroupAnalysisTable.optoParameters{1,1}.activeChannels{2};
        activeChannels = {[activeChannels{1}, ', ', secondActiveChannel{1}]};
    
    elseif  numel(runGroupAnalysisTable.optoParameters{1,1}.activeChannels) == 3
        
        secondActiveChannel = runGroupAnalysisTable.optoParameters{1,1}.activeChannels{2};
        thirdActiveChannel = runGroupAnalysisTable.optoParameters{1,1}.activeChannels{3};
        activeChannels = {[activeChannels{1}, ', ', secondActiveChannel{1}, ', ', thirdActiveChannel{1}]};
    
    end
    
    dataCellToAdd = [{iRunGroup}, {sampleSizeRunGroup}, activeChannels, optoParametersCell(2:end)'];
    
    runGroupsParametersTable = [runGroupsParametersTable; dataCellToAdd];
    
    saveMouseParametersTableName = 'MouseParametersTable.mat';
    saveMouseParametersTableName = strrep(saveMouseParametersTableName, ' ', '');
    saveMouseParametersTablePath = fullfile(analyzedMouseDataFolderPath, saveMouseParametersTableName);
    save(saveMouseParametersTablePath, 'runGroupsParametersTable');
    
end
