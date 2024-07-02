clear all;
close all;
clc;

% User environment (automatic)
user = getenv('USER'); 

% Indicate main directory and add it to path
matlabDirectory = ['/Users/' user '/Desktop/Thesis/Matlab/'];
mainDirectory = ['/Users/' user '/Desktop/Thesis/Matlab/Experiments/']; % including all Matlab toolboxes and necessary paths
addpath(genpath(matlabDirectory));

% Indicate experiment to analyze
experimentName = '';

% Set experiment path
experimentDirectory = [mainDirectory experimentName];

folderContent = dir([experimentDirectory]);
folderContent = folderContent([folderContent.isdir]);
nCells = sum([folderContent.isdir]) - 2;
startDir = 3;

loadCells = struct;
selectedCells = ["cell1","cell2"];

warning('off','MATLAB:unknownObjectNowStruct')

 for iCell = startDir:(startDir+nCells-1)
     
     currentCell = folderContent(iCell);
     
     for iSelectedCell = 1:numel(selectedCells)  

        if isequal(currentCell.name, selectedCells(iSelectedCell))
            
            disp(['****** Loaded ' currentCell.name ' ******'])
            cellPath = [experimentDirectory filesep folderContent(iCell).name];
            processedDataPath = [cellPath filesep 'ProcessedData' filesep 'Data'];
            analysisData = [];
            [spanEpochs, spanCyclePositions] = findSpanEpochsCyclePositions(processedDataPath);
            
            columnNames = {'epoch' 'cyclePosition', 'restMean', 'restMedian', 'restSD', 'restMin', 'restMax', 'pulseAP', 'postAP','reboundAP', 'areaVtPulse', 'areaVtControl', 'nAPtotal', 'nAPpreLight', 'nAPduringLight', 'nAPpostLight', 'rateAPpreLight', 'rateAPduringLight', 'rateAPpostLight', 'delayFirstAPpostLight', 'trace', 'optoParameters'};
            dataTypes = {'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell','cell', 'cell', 'cell', 'cell', 'cell', 'cell'};
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
                    phAnalysis = DataStruct.phAnalysis; 
                    runInfo = [iEpoch, iCyclePosition];

                    activeChannels = regexp(DataStruct.state.phys.internal.lastLinesUsed, '''(\w+)''', 'tokens');
                    
                    optoParameters = [];
                    optoParameters.activeChannels = activeChannels;

                    if any(cellfun(@(x) contains(x, 'ao0'), activeChannels))

                        nCurrentPulses = str2double(DataStruct.state.phys.internal.pulses.pulseString_ao0.numPulses);
                        delayCurrentPulse = str2double(DataStruct.state.phys.internal.pulses.pulseString_ao0.delay)*10;
                        isiCurrentPulse = str2double(DataStruct.state.phys.internal.pulses.pulseString_ao0.isi);
                        currentPulseWidth = str2double(DataStruct.state.phys.internal.pulses.pulseString_ao0.pulseWidth)*10;
                        currentAmplitude = str2double(DataStruct.state.phys.internal.pulses.pulseString_ao0.amplitude);

                        optoParameters.delayCurrentPulse = delayCurrentPulse;
                        optoParameters.widthCurrentPulse = currentPulseWidth;
                        optoParameters.amplitudeCurrentPulse = currentAmplitude;

                    else
                        
                        optoParameters.delayCurrentPulse = nan;
                        optoParameters.widthCurrentPulse = nan;
                        optoParameters.amplitudeCurrentPulse = nan;

                    end

                    if any(cellfun(@(x) contains(x, 'ao1'), activeChannels))

                        nPulsesBlue = str2double(DataStruct.state.phys.internal.pulses.pulseString_ao1.numPulses);
                        delayFirstPulseBlue = str2double(DataStruct.state.phys.internal.pulses.pulseString_ao1.delay)*10;
                        isiBlue = str2double(DataStruct.state.phys.internal.pulses.pulseString_ao1.isi);
                        pulseWidthBlue = str2double(DataStruct.state.phys.internal.pulses.pulseString_ao1.pulseWidth)*10;
                        amplitudeBlue = str2double(DataStruct.state.phys.internal.pulses.pulseString_ao1.amplitude);
                        timeArray = linspace(0,size(DataStruct.data,2)-1,size(DataStruct.data,2));
                        pulsesTimeArray = (delayFirstPulseBlue + linspace(0,(nPulsesBlue - 1)*isiBlue,nPulsesBlue))*10;
                        pulsesWidth = pulseWidthBlue;
                        pulsesTimeDifference = nan;
                        whichPulseFirst = nan;
                        functionNameBlue = DataStruct.state.cycle.functionName;
                        functionNameBlue = strrep(functionNameBlue, '''', '');
                        functionNameBlue = strrep(functionNameBlue, ';', '');

                        optoParameters.nLightPulses = nPulsesBlue;
                        optoParameters.widthLightPulse = pulseWidthBlue;
                        optoParameters.amplitudeLightPulse = amplitudeBlue;
                        optoParameters.delayLightPulse = delayFirstPulseBlue;
                        
                        APdistribution = computeAPdistribution(DataStruct, optoParameters);

                    elseif any(cellfun(@(x) contains(x, 'ao2'), activeChannels))

                        nPulsesRed = str2double(DataStruct.state.phys.internal.pulses.pulseString_ao2.numPulses);
                        delayFirstPulseRed = str2double(DataStruct.state.phys.internal.pulses.pulseString_ao2.delay)*10;
                        isiRed = str2double(DataStruct.state.phys.internal.pulses.pulseString_ao2.isi);
                        pulseWidthRed = str2double(DataStruct.state.phys.internal.pulses.pulseString_ao2.pulseWidth)*10;
                        amplitudeRed = str2double(DataStruct.state.phys.internal.pulses.pulseString_ao2.amplitude);                    
                        timeArray = linspace(0,size(DataStruct.data,2)-1,size(DataStruct.data,2));
                        pulsesTimeArray = (delayFirstPulseRed + linspace(0,(nPulsesRed - 1)*isiRed,nPulsesRed))*10;
                        pulsesWidth = pulseWidthRed;
                        pulsesTimeDifference = nan;
                        whichPulseFirst = nan;
                        functionNameRed = DataStruct.state.cycle.functionName;
                        functionNameRed = strrep(functionNameRed, '''', '');
                        functionNameRed = strrep(functionNameRed, ';', '');

                        optoParameters.nLightPulses = nPulsesRed;
                        optoParameters.widthLightPulse = pulseWidthRed;
                        optoParameters.amplitudeLightPulse = amplitudeRed;
                        optoParameters.delayLightPulse = delayFirstPulseRed;
                        
                        APdistribution = computeAPdistribution(DataStruct, optoParameters);

                    else
                        
                        optoParameters.nLightPulses = 0;
                        optoParameters.widthLightPulse = 5000;
                        optoParameters.amplitudeLightPulse = 0;
                        optoParameters.delayLightPulse = 25500;
                        
                        APdistribution = computeAPdistribution(DataStruct, optoParameters);                        

                    end
                    
                    [areaVtPulse, areaVtControl] = computeAreaVtCC(DataStruct, optoParameters);

                    runDataStruct = [];                
                    runDataStruct = DataStruct.phAnalysis; 
                    fieldsToRemove = {'acqs','checkPulseRpeak','checkPulseRend','checkPulseTau', 'deltaI', 'pulseV', 'pulseAHP','pulseRm', 'nAP'};
                    runDataStruct = rmfield(runDataStruct, fieldsToRemove);
                    runDataStruct.areaVtPulse = areaVtPulse;
                    runDataStruct.areaVtControl = areaVtControl;
                    
                    fieldsAPdistribution = fields(APdistribution);
                    
                    for iFieldAPdistribution = 1: numel(fieldsAPdistribution)
                    
                        runDataStruct = setfield(runDataStruct, fieldsAPdistribution{iFieldAPdistribution}, APdistribution.(fieldsAPdistribution{iFieldAPdistribution}));
                        
                    end

                    runDataStruct.trace = DataStruct.data;
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

            analysisTable = findIdenticalRunsCC(analysisTable);
            saveTableName = 'AnalysisTable.mat';
            saveTableName = strrep(saveTableName, ' ', '');
            saveTablePath = fullfile(analyzedDataFolderPath, saveTableName);
            save(saveTablePath, 'analysisTable');

            disp(['****** Results saved successfully ******']);

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
                mouseAnalysisTable = findIdenticalRunsCC(mouseAnalysisTable);
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
                mouseAnalysisTable = findIdenticalRunsCC(mouseAnalysisTable);
                save(saveMouseTablePath, 'mouseAnalysisTable');

            end
            
        end
        
     end
     
 end

mouseAnalysisTable = load(saveMouseTablePath);
mouseAnalysisTable = mouseAnalysisTable.mouseAnalysisTable;
 
runGroups = unique(mouseAnalysisTable.runGroup(:));
nRunGroups = numel(runGroups);

columnNames = {'runGroup', 'sampleSize', 'activeChannels', 'delayCurrentPulse', 'widthCurrentPulse', 'amplitudeCurrentPulse', 'nLightPulses', 'widthLightPulse', 'amplitudeLightPulse', 'delayLightPulse'};
dataTypes = {'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell'};
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
    
    end
    
    dataCellToAdd = [{iRunGroup}, {sampleSizeRunGroup}, activeChannels, optoParametersCell(2:end)'];
    
    runGroupsParametersTable = [runGroupsParametersTable; dataCellToAdd];
    
    saveMouseParametersTableName = 'MouseParametersTable.mat';
    saveMouseParametersTableName = strrep(saveMouseParametersTableName, ' ', '');
    saveMouseParametersTablePath = fullfile(analyzedMouseDataFolderPath, saveMouseParametersTableName);
    save(saveMouseParametersTablePath, 'runGroupsParametersTable');
    
end

