disp(['*-*-*-* Running: scriptOrganizeAcquisitionsCC *-*-*-*'])

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
            
            phAnalysisPath = [cellPath filesep 'RawData' filesep '_phAnalysis.mat'];
            phAnalysisFile = load(phAnalysisPath);
            phAnalysisFile = phAnalysisFile.phAnalysis;

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

                responseType = 'CC';
                phAnalysisData = phAnalysisExtractData(phAnalysisFile, runEpoch, existingRuns{1,iRun});
                
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
                DataStruct.phAnalysis = phAnalysisData;
                
                saveFolderPath = [processedDataFolderPath filesep];
                
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

            end
            
            textprogressbar(' Completed and saved');
            
        end
        
     end
     
 end