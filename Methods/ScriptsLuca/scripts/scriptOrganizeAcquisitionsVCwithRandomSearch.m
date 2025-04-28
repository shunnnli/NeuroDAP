disp(['*-*-*-* Running: scriptOrganizeAcquisitionsVCwithRandomSearch *-*-*-*'])

 for iCell = startDir:(startDir+nCells-1)
     
     currentCell = folderContent(iCell);
     
     for iSelectedCell = 1:numel(selectedCells)  

        if isequal(currentCell.name, selectedCells(iSelectedCell))
            
            cellPath = [experimentDirectory filesep folderContent(iCell).name];
            
            if ~exist([cellPath filesep 'RawData'], 'dir'); moveToRawDataFolder(cellPath); end
            if exist([cellPath filesep 'ProcessedData'], 'dir') && overwriteProcessing == 0; continue; end
            
            rawDataFolderPath = [cellPath filesep 'RawData'];
            addpath(genpath(rawDataFolderPath));
            
            fileList = dir(fullfile(rawDataFolderPath, 'AD0_*.mat'));
            existingRuns = {};
            
            for iFile = 1:length(fileList)
                
                [~, fileName, runExt] = fileparts(fileList(iFile).name);
                    
                if ~isempty(regexp(fileName, '^AD0_\d+$', 'once')) && ~contains(fileName, 'AD0_e')
                    
                     existingRuns{end+1} = fileList(iFile).name;
                     
                end
                
            end
             
            processedDataFolderName = 'ProcessedData';
            processedDataFolderPath = [cellPath filesep processedDataFolderName filesep 'Data'];

            if exist(processedDataFolderPath, 'dir') == 0
                
                mkdir(processedDataFolderPath);
                addpath(genpath(processedDataFolderPath));
                
            elseif exist(processedDataFolderPath, 'dir') && overwriteProcessing == 1
                
                oldFileList = dir(fullfile(processedDataFolderPath, '*')); 

                for i = 1:length(oldFileList)

                    if ~oldFileList(i).isdir

                        oldFilePath = fullfile(processedDataFolderPath, oldFileList(i).name);
                        delete(oldFilePath);
                    
                    end
                    
                end
                
            end
            
            textprogressbar([currentCell.name ' - Organizing acquisitions: ']);
            
            for iRun = 1:numel(existingRuns)
                
                textprogressbar(iRun,numel(existingRuns));
                
                runPath = [rawDataFolderPath filesep existingRuns{1,iRun}];              
                headerStringStructure = generateHeaderStringStructure(runPath);
                
                runEpoch = string(headerStringStructure.state.epoch);
                runEpoch = regexp(runEpoch, '\d+', 'match');
                currentCyclePosition = string(headerStringStructure.state.cycle.currentCyclePosition);
                currentCyclePosition = regexp(currentCyclePosition, '\d+', 'match');
                useCyclePosition = string(headerStringStructure.state.cycle.useCyclePos);
                useCyclePosition = regexp(useCyclePosition, '\d+', 'match');
                
                excelPath = [cellPath filesep 'InfoPatching.xlsx'];
                runCyclePosition = findCorrectCyclePosition(excelPath, runEpoch, currentCyclePosition, useCyclePosition);
                
                [generalDrugs, cellLocation, cellTemperature, cellType, cellFluo] = findCellInfoExcel(excelPath);
                [toBeExcluded, holdingVoltage, specificDrugs, qualityRS0, scopeMagnification] = findRunInfoExcel(excelPath, runEpoch, runCyclePosition);
                cellCoordinatesBlue = extractCellCoordinates(headerStringStructure.state.zDMD.tVectorBlue);
                cellCoordinatesRed = extractCellCoordinates(headerStringStructure.state.zDMD.tVectorRed);
                
                fileNumber = str2num(headerStringStructure.state.files.fileCounter);
                fileNumberString = num2str(fileNumber);
                holdingVoltageDataPath = [rawDataFolderPath filesep 'AD1_' fileNumberString '.mat'];
                
                if isfile(holdingVoltageDataPath) 
                    
                    holdingVoltageDataFile = load(holdingVoltageDataPath);
                    holdingVoltageDataFile = holdingVoltageDataFile.(['AD1_' fileNumberString]);
                    holdingVoltageValue = extractHoldingVoltage(holdingVoltageDataFile);
                    holdingVoltage = holdingVoltageValue;
                    
                end
                
                if toBeExcluded; continue; end
                
                activeChannels = regexp(headerStringStructure.state.phys.internal.lastLinesUsed, '''(\w+)''', 'tokens');
                stringArrayActiveChannels = cellfun(@(x) strjoin(x,''), activeChannels, 'UniformOutput', false);
                activeChannels = strjoin(stringArrayActiveChannels,', '); 
                
                if strcmp(responseTypeBlue,'excitatory') && strcmp(activeChannels,'ao0, ao1')
                    
                    if holdingVoltage == -70; responseType = 'excitatory';
                    elseif holdingVoltage >= 0; responseType = 'none';
                    else responseType = 'undefined'; end
                    
                elseif strcmp(responseTypeBlue,'inhibitory') && strcmp(activeChannels,'ao0, ao1')
                    
                    if holdingVoltage == -70; responseType = 'none';
                    elseif holdingVoltage >= 0; responseType = 'inhibitory';
                    else responseType = 'undefined'; end
                    
                elseif strcmp(responseTypeRed,'excitatory') && strcmp(activeChannels,'ao0, ao2')
                    
                    if holdingVoltage == -70; responseType = 'excitatory';
                    elseif holdingVoltage >= 0; responseType = 'none';
                    else responseType = 'undefined'; end
                    
                elseif strcmp(responseTypeRed,'inhibitory') && strcmp(activeChannels,'ao0, ao2')
                    
                    if holdingVoltage == -70; responseType = 'none';
                    elseif holdingVoltage >= 0; responseType = 'inhibitory';
                    else responseType = 'undefined'; end
                    
                elseif strcmp(activeChannels,'ao0, ao1, ao2') && holdingVoltage == -35
                    
                    responseType = 'both';
                    
                else
                    
                    responseType = 'none';
                    
                end

                
                cyclePositionData = load(runPath);
                [~, baseName] = fileparts(runPath);
                
                DataStruct = headerStringStructure;
                DataStruct.epoch = runEpoch;
                DataStruct.cyclePosition = runCyclePosition;
                DataStruct.holdingVoltage = holdingVoltage;
                DataStruct.scopeMagnification = scopeMagnification;
                DataStruct.generalDrugs = generalDrugs;
                DataStruct.specificDrugs = specificDrugs;
                DataStruct.cellLocation = cellLocation;
                DataStruct.cellTemperature = cellTemperature;
                DataStruct.cellType = cellType;
                DataStruct.cellFluo = cellFluo;
                DataStruct.cellCoordinates = [];
                DataStruct.cellCoordinates.blue = cellCoordinatesBlue;
                DataStruct.cellCoordinates.red = cellCoordinatesRed;
                DataStruct.responseType = responseType;
                DataStruct.qualityRS0 = qualityRS0;
                DataStruct.data = cyclePositionData.(baseName).data;
                
                saveFolderPath = [processedDataFolderPath filesep];
                
                if isempty(regexp(DataStruct.state.zDMD.searchDepth, '[a-zA-Z0-9]', 'once')) || ~contains(DataStruct.state.cycle.cycleName,'randomSearch')
                    
                    DataStruct.state.zDMD.searchDepth = NaN;
                    DataStruct.state.zDMD.searchRepetition = NaN;
                    
                    saveRunName = ['Epoch' char(runEpoch) '_cyclePosition' char(runCyclePosition) '.mat'];
                    saveRunName = strrep(saveRunName, ' ', '');
                    saveRunPath = fullfile(saveFolderPath, saveRunName);
                    
                    incrementCyclePosition = 1;
                    
                    while exist(saveRunPath,'file') == 2
                        
                        saveRunName = ['Epoch' char(runEpoch) '_cyclePosition' num2str(str2double(runCyclePosition)+incrementCyclePosition) '.mat'];
                        saveRunName = strrep(saveRunName, ' ', '');
                        saveRunPath = fullfile(saveFolderPath, saveRunName);
                        incrementCyclePosition = incrementCyclePosition + 1;
                    
                    end
                    
                    save(saveRunPath, 'DataStruct');
                    
                else
                    
                    runSearchDepth = string(DataStruct.state.zDMD.searchDepth);
                    runSearchDepth = regexp(runSearchDepth, '\d+', 'match');
                    runSearchRepetition = string(DataStruct.state.zDMD.searchRepetition);
                    runSearchRepetition = regexp(runSearchRepetition, '\d+', 'match');
                    
                    saveRunName = ['Epoch' char(runEpoch) '_searchDepth' char(runSearchDepth) '_searchRepetition' char(runSearchRepetition) '.mat'];
                    saveRunName = strrep(saveRunName, ' ', '');
                    saveRunPath = fullfile(saveFolderPath, saveRunName);
                    save(saveRunPath, 'DataStruct');
                    
                end
                
            end
            
            textprogressbar(' Completed and saved');
            
        end
        
     end
     
 end