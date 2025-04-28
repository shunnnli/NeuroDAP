% This script should be run for each experiment folder individually

clear all;
clc;

% User environment (automatic)
user = getenv('USER'); 

% Indicate main directory and add it to path
matlabDirectory = ['/Users/' user '/Desktop/Thesis/Matlab/'];
mainDirectory = ['/Users/' user '/Desktop/Thesis/Matlab/Experiments/']; % including all Matlab toolboxes and necessary paths
addpath(genpath(matlabDirectory));

% Indicate experiment to analyze
experimentName = 'vGAT_ZI_Ephys_11';

% Set experiment path
experimentDirectory = [mainDirectory experimentName];

folderContent = dir([experimentDirectory]);
folderContent = folderContent([folderContent.isdir]);
nCells = sum([folderContent.isdir]) - 2;
startDir = 3;

loadCells = struct;
selectedCells = ["cell1","cell2","cell3","cell4","cell5","cell6"];

warning('off','MATLAB:unknownObjectNowStruct')

 for iCell = startDir:(startDir+nCells-1)
     
     currentCell = folderContent(iCell);
     
     for iSelectedCell = 1:numel(selectedCells)  

        if isequal(currentCell.name, selectedCells(iSelectedCell))
            
            disp(['*** Loaded ' currentCell.name ' ***'])
            cellPath = [experimentDirectory filesep folderContent(iCell).name];
            
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
                
                saveRunName = ['Epoch' char(runEpoch) '_cyclePosition' char(runCyclePosition) '.mat'];
                saveRunName = strrep(saveRunName, ' ', '');
                saveRunPath = fullfile(saveFolderPath, saveRunName);
                save(saveRunPath, 'DataStruct');

                disp([saveRunName, ' saved successfully']);

            end
            
        end
        
     end
     
 end