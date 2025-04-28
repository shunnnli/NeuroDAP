disp(['*-*-*-* Running: scriptCollectFeaturesVC *-*-*-*'])

 for iCell = startDir:(startDir+nCells-1)
     
     currentCell = folderContent(iCell);
     
     for iSelectedCell = 1:numel(selectedCells)  

        if isequal(currentCell.name, selectedCells(iSelectedCell))
            
            cellPath = [experimentDirectory filesep folderContent(iCell).name];
            processedDataPath = [cellPath filesep 'ProcessedData' filesep 'Data'];
            analysisData = [];
            [spanEpochs, spanCyclePositions] = findSpanEpochsCyclePositions(processedDataPath);
            
            columnNames = {'epoch' 'cyclePosition', 'heightPulsePeak', 'timePulsePeak', 'isPulsePeakAboveThreshold', 'areaPulse', 'heightControlPeak', 'areaControl', 'trace', 'qualityRS0', 'runProtocol'};
            dataTypes = {'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell'};
            analysisTable = table('Size', [0, numel(columnNames)], 'VariableNames', columnNames, 'VariableTypes', dataTypes);
            
%             for i = 1:numel(columnNames)
%                 
%                 analysisTable.(columnNames{i}) = repmat({[]}, 0, 1);
%                 
%                 if strcmp(dataTypes{i}, 'double')
%                     
%                     analysisTable.(columnNames{i}) = NaN(0, 1);
%                 
%                 elseif strcmp(dataTypes{i}, 'logical')
%                     
%                     analysisTable.(columnNames{i}) = false(0, 1);
%                 
%                 end
%                 
%             end
            
            textprogressbar([currentCell.name ' - Processing acquisitions: ']);            
            runFiles = dir(processedDataPath);
            filesWithCyclePosition = runFiles(~[runFiles.isdir] & contains({runFiles.name}, 'cyclePosition')); 
            nRunFiles = numel(filesWithCyclePosition);
            iRunProgress = 1;
            
            for iCyclePosition = spanCyclePositions
                
                nameCyclePosition = strcat(['cyclePosition', num2str(iCyclePosition)]);
                analysisData.(nameCyclePosition) = [];
                
                for iEpoch = spanEpochs
                    
                    nameEpoch = strcat(['epoch', num2str(iEpoch)]);

                    runPath = [cellPath filesep 'ProcessedData' filesep 'Data' filesep 'Epoch' num2str(iEpoch) '_cyclePosition' num2str(iCyclePosition) '.mat'];

                    if ~isfile(runPath)
                        continue;
                    end
                    
                    textprogressbar(iRunProgress,nRunFiles); iRunProgress = iRunProgress + 1;
                    runData = load(runPath);
                    DataStruct = runData.DataStruct;
                    runInfo = [iEpoch, iCyclePosition];

                    data = DataStruct.data;
                    activeChannels = regexp(DataStruct.state.phys.internal.lastLinesUsed, '''(\w+)''', 'tokens');
                    
                    runProtocol = [];
                    runProtocol.activeChannels = activeChannels;
                    runProtocol.holdingVoltage = DataStruct.holdingVoltage;
                    runProtocol.scopeMagnification = DataStruct.scopeMagnification;
                    runProtocol.generalDrugs = DataStruct.generalDrugs;
                    runProtocol.specificDrugs = DataStruct.specificDrugs;                  
                    runProtocol.cellLocation = DataStruct.cellLocation;
                    runProtocol.cellTemperature = DataStruct.cellTemperature;
                    runProtocol.cellType = DataStruct.cellType;
                    runProtocol.cellFluo = DataStruct.cellFluo;
                    runProtocol.cellCoordinates = DataStruct.cellCoordinates;
                    runProtocol.responseType = DataStruct.responseType;
 
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

                        runProtocol.nPulsesBlue = nPulsesBlue;
                        runProtocol.isiBlue = isiBlue;
                        runProtocol.pulseWidthBlue = pulseWidthBlue;
                        runProtocol.amplitudeBlue = amplitudeBlue;
                        runProtocol.delayPulseBlue = delayFirstPulseBlue;
                        runProtocol.functionNameBlue = functionNameBlue;

                        runProtocol.nPulsesRed = nPulsesRed;
                        runProtocol.isiRed = isiRed;
                        runProtocol.pulseWidthRed = pulseWidthRed;
                        runProtocol.amplitudeRed = amplitudeRed;
                        runProtocol.delayPulseRed = delayFirstPulseRed;
                        runProtocol.functionNameRed = functionNameRed;

                        runProtocol.whichPulseFirst = whichPulseFirst;
                        runProtocol.pulsesTimeDifference = pulsesTimeDifference;
                          
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
                        
                        runProtocol.nPulsesBlue = nPulsesBlue;
                        runProtocol.isiBlue = isiBlue;
                        runProtocol.pulseWidthBlue = pulseWidthBlue;
                        runProtocol.amplitudeBlue = amplitudeBlue;
                        runProtocol.delayPulseBlue = delayFirstPulseBlue;
                        runProtocol.functionNameBlue = functionNameBlue;

                        runProtocol.nPulsesRed = nan;
                        runProtocol.isiRed = nan;
                        runProtocol.pulseWidthRed = nan;
                        runProtocol.amplitudeRed = nan;
                        runProtocol.delayPulseRed = nan;
                        runProtocol.functionNameRed = 'NaN';

                        runProtocol.whichPulseFirst = whichPulseFirst;
                        runProtocol.pulsesTimeDifference = pulsesTimeDifference;                   

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

                        runProtocol.nPulsesBlue = nan;
                        runProtocol.isiBlue = nan;
                        runProtocol.pulseWidthBlue = nan;
                        runProtocol.amplitudeBlue = nan;
                        runProtocol.delayPulseBlue = nan;
                        runProtocol.functionNameBlue = 'NaN';

                        runProtocol.nPulsesRed = nPulsesRed;
                        runProtocol.isiRed = isiRed;
                        runProtocol.pulseWidthRed = pulseWidthRed;
                        runProtocol.amplitudeRed = amplitudeRed;
                        runProtocol.delayPulseRed = delayFirstPulseRed;
                        runProtocol.functionNameRed = functionNameRed;

                        runProtocol.whichPulseFirst = whichPulseFirst;
                        runProtocol.pulsesTimeDifference = pulsesTimeDifference;

                    end
                    
                    if isfield(DataStruct.state.zDMD,'options'); runProtocol.options = DataStruct.state.zDMD.options;
                    else; runProtocol.options = struct; end
                    
                    if length(pulsesTimeArray) == 1
                                             
                        timeBeforePulse = 5000;
                        timeAfterPulse = 14999;
                        
                    else
                        
                        timeBeforePulse = min(5000,(pulsesTimeArray(2) - pulsesTimeArray(1))/2);
                        timeAfterPulse = min(14999,(pulsesTimeArray(2) - pulsesTimeArray(1))/2 - 1);
                        
                    end

                    segments = cell(1, length(pulsesTimeArray));

                    for iPulse = 1:length(pulsesTimeArray)

                        pulseIndex = round(pulsesTimeArray(iPulse));
                        startIdx = pulseIndex - timeBeforePulse;
                        endIdx = pulseIndex + timeAfterPulse;
                        startIdx = max(1, startIdx);
                        endIdx = min(length(data), endIdx);
                                               
                        % No overlap in case standard RC parameters are used
                        timeAfterRCCheck = 5000;
                        startFirstBaseline = min(timeAfterRCCheck,pulsesTimeArray(1) - 1000);
                        endFirstBaseline = pulsesTimeArray(1);
                        startSecondBaseline = min(pulsesTimeArray(end) + timeAfterPulse,length(data) - 1000);
                        endSecondBaseline = length(data);

                        % Average and SD of combined baseline windows
                        baselineAverage = mean([data(startFirstBaseline:endFirstBaseline),data(startSecondBaseline:endSecondBaseline)]);
                        baselineSD = std([data(startFirstBaseline:endFirstBaseline),data(startSecondBaseline:endSecondBaseline)]);

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
                    
                    if DataStruct.holdingVoltage > -20

                        responseSign = 1;

                    elseif DataStruct.holdingVoltage < -60 

                        responseSign = -1;
                        
                    else 
                        
                        responseSign = 2;

                    end

                    for iPulse = 1:length(segments)

                        peakData = scriptPerformAnalysisVCtwoSigns(segments{iPulse}, baselineSD, timeBeforePulse, responseSign);

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
                    
                    runDataStruct.qualityRS0 = DataStruct.qualityRS0;
                    runDataStruct.runProtocol = runProtocol;

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
            
            textprogressbar(' Completed and saved');

            saveResultsName = 'AnalysisStruct.mat';
            saveResultsName = strrep(saveResultsName, ' ', '');
            saveResultsPath = fullfile(analyzedDataFolderPath, saveResultsName);
            save(saveResultsPath, 'analysisData');

            analysisTable = findIdenticalRunsVC(analysisTable);
            saveTableName = 'AnalysisTable.mat';
            saveTableName = strrep(saveTableName, ' ', '');
            saveTablePath = fullfile(analyzedDataFolderPath, saveTableName);
            save(saveTablePath, 'analysisTable');

            analyzedMouseDataFolderName = 'mouseAnalysis';
            analyzedMouseDataFolderPath = [experimentDirectory filesep analyzedMouseDataFolderName];

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
                
                cellNumber = regexp(currentCell.name, '\d+', 'match');
                cellNumber = str2double(cellNumber);

                if any(cell2mat(existingMouseAnalysisTable.cellName) == cellNumber) && overwriteAnalysis == 0

                    continue;

                end
                
                rowsCell = find(cell2mat(existingMouseAnalysisTable.cellName) == cellNumber);
                existingMouseAnalysisTable(rowsCell,:) = [];

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

columnNames = {'runGroup', 'sampleSize', 'activeChannels', 'holdingVoltage', 'scopeMagnification', 'generalDrugs', 'specificDrugs', 'cellLocation', 'cellTemperature', 'cellType', 'cellFluo', 'cellCoordinates', 'responseType', 'nPulsesBlue', 'isiBlue', 'pulseWidthBlue', 'amplitudeBlue', 'delayPulseBlue', 'functionNameBlue', 'nPulsesRed', 'isiRed', 'pulseWidthRed', 'amplitudeRed', 'delayPulseRed', 'functionNameRed', 'whichPulseFirst', 'pulsesTimeDifference', 'options'};
dataTypes = {'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell'};
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
    
    runProtocol = runGroupAnalysisTable.runProtocol{1};
    runProtocolCell = struct2cell(runProtocol);
    
    activeChannels = mouseAnalysisTable.runProtocol{1,1}.activeChannels{1};
    
    if numel(runGroupAnalysisTable.runProtocol{1,1}.activeChannels) == 2
        
        secondActiveChannel = runGroupAnalysisTable.runProtocol{1,1}.activeChannels{2};
        activeChannels = {[activeChannels{1}, ', ', secondActiveChannel{1}]};
    
    elseif  numel(runGroupAnalysisTable.runProtocol{1,1}.activeChannels) == 3
        
        secondActiveChannel = runGroupAnalysisTable.runProtocol{1,1}.activeChannels{2};
        thirdActiveChannel = runGroupAnalysisTable.runProtocol{1,1}.activeChannels{3};
        activeChannels = {[activeChannels{1}, ', ', secondActiveChannel{1}, ', ', thirdActiveChannel{1}]};
    
    end
    
    dataCellToAdd = [{iRunGroup}, {sampleSizeRunGroup}, activeChannels, runProtocolCell(2:end)'];
    
    runGroupsParametersTable = [runGroupsParametersTable; dataCellToAdd];
    
    saveMouseParametersTableName = 'MouseParametersTable.mat';
    saveMouseParametersTableName = strrep(saveMouseParametersTableName, ' ', '');
    saveMouseParametersTablePath = fullfile(analyzedMouseDataFolderPath, saveMouseParametersTableName);
    save(saveMouseParametersTablePath, 'runGroupsParametersTable');
    
end
