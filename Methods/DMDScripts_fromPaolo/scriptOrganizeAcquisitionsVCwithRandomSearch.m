disp(['*-*-*-* Running: scriptOrganizeAcquisitionsVCwithRandomSearch *-*-*-*'])

 for iCell = startDir:(startDir+length(cellList)-1)
     
     currentCell = cellList(iCell);
     
     for iSelectedCell = 1:numel(selectedCells)  

        if isequal(currentCell.name, selectedCells(iSelectedCell))
            
            disp(['*** Loaded ' currentCell.name ' ***'])
            cellPath = [expPath filesep cellList(iCell).name];
            
            moveToRawDataFolder(cellPath);
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
                disp(['Folder "', processedDataFolderName, '" created successfully.']);
                
            end
            
            for iRun = 1:numel(existingRuns)
                
                runPath = [rawDataFolderPath filesep existingRuns{1,iRun}];              
                headerStringStructure = generateHeaderStringStructure(runPath);
                disp(['headerString structure for ', existingRuns{1,iRun}, ' created successfully']);
                
                runEpoch = string(headerStringStructure.state.epoch);
                runEpoch = regexp(runEpoch, '\d+', 'match');
                runCyclePosition = string(headerStringStructure.state.cycle.currentCyclePosition);
                runCyclePosition = regexp(runCyclePosition, '\d+', 'match');
                
                excelPath = [cellPath filesep 'InfoPatching.xlsx'];
                toBeExcluded = findExcludedRuns(excelPath, runEpoch, runCyclePosition);
                holdingVoltage = findHoldingVoltage(excelPath, runEpoch, 1);
                
                if toBeExcluded
                    
                    continue;
                    
                end
                
                cyclePositionData = load(runPath);
                [~, baseName] = fileparts(runPath);
                
                DataStruct = headerStringStructure;
                DataStruct.epoch = runEpoch;
                DataStruct.cyclePosition = runCyclePosition;
                DataStruct.holdingVoltage = holdingVoltage;
                DataStruct.data = cyclePositionData.(baseName).data;
                saveFolderPath = [processedDataFolderPath filesep];
                
                if isempty(regexp(DataStruct.state.zDMD.searchDepth, '[a-zA-Z0-9]', 'once'))
                    
                    DataStruct.state.zDMD.searchDepth = NaN;
                    DataStruct.state.zDMD.searchRepetition = NaN;
                    
                    saveRunName = ['Epoch' char(runEpoch) '_cyclePosition' char(runCyclePosition) '.mat'];
                    saveRunName = strrep(saveRunName, ' ', '');
                    saveRunPath = fullfile(saveFolderPath, saveRunName);
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
                
                    disp([saveRunName, ' saved successfully']);

            end
            
        end
        
     end
     
 end