%% Load sessions
clear; close all;
addpath(genpath(osPathSwitch('/Volumes/MICROSCOPE/Shun/Analysis/NeuroDAP/Methods')));

% Select sessions for analysis
% parentPath = osPathSwitch('/Volumes/MICROSCOPE/wengang/Exp_withShun/');
parentPath = osPathSwitch('/Volumes/MICROSCOPE/Shun/Project valence/Patch/');
expPath = uipickfiles('FilterSpec',parentPath,'Prompt','Select experiment folders');
rawDataPath = 'default'; % strcat(parentPath,filesep,'20231221_ally');

% Set comman params
[sessionParams,canceled] = inputSessionParams_singleSlice(expPath,...
                                paradigm=2,redStim=true,...
                                reload=false,calculateQC=false,...
                                timeRange='[-20,100]',nArtifactSamples='10');
taskOptions = {'random','reward pairing','punish pairing'};
task = taskOptions{sessionParams.Paradigm};
timeRange = eval(sessionParams.timeRange);
nArtifactSamples = str2double(sessionParams.nArtifactSamples);

if ~isscalar(expPath); error('Multiple sessions were selected, for multi-session analysis see analyzeSlice pipeline!'); end
expPath = expPath{1};

%% Luca's code

cellList = dir(fullfile(expPathLuca,'cell*'));
nCells = length(cellList);
startDir = 1;

loadCells = struct;
selectedCells = ["cell1"];

queryFeature = 'areaPulse';
savePlots = 1;
responseTimeWindow = num2str(500);

warning('off','MATLAB:unknownObjectNowStruct')

 for iCell = startDir:(startDir+nCells-1)
     
     currentCell = cellList(iCell);
     
     for iSelectedCell = 1:numel(selectedCells)  

        if isequal(currentCell.name, selectedCells(iSelectedCell))
            
            disp(['****** Loaded ' currentCell.name ' ******'])
            cellPath = [expPathLuca filesep cellList(iCell).name];
            processedDataPath = [cellPath filesep 'ProcessedData' filesep 'Data'];
            analyzedDataFolderName = 'Analysis';
            analyzedDataFolderPath = [cellPath filesep analyzedDataFolderName];
            savePlotPath = analyzedDataFolderPath;
            [spanEpochs, spanSearchDepths, spanSearchRepetitions] = findSpanEpochsRandomSearch(processedDataPath);

            for iEpoch = spanEpochs
                
                pathAnalysisTable = [analyzedDataFolderPath filesep 'searchAnalysisTable_Epoch' num2str(iEpoch) '.mat'];
                
                if ~isfile(pathAnalysisTable)
                    continue;
                end
                
                epochAnalysisTable = load(pathAnalysisTable);
                epochAnalysisTable = epochAnalysisTable.analysisTable;
                holdingVoltage = epochAnalysisTable.optoParameters{1,1}.holdingVoltage;
                
                if holdingVoltage > -30

                    responseSign = 1;

                else 

                    responseSign = -1;

                end
                
                if isempty(epochAnalysisTable)
                    continue;
                end
                
                overallResponseMap = nan(684,608);
            
                for iDepth = spanSearchDepths
                    
                    depthAnalysisTable = epochAnalysisTable(cell2mat(epochAnalysisTable.depth) == iDepth,:);
                    depthResponseMap = nan(684,608);
                    
                    if isempty(depthAnalysisTable)
                        continue;
                    end

                    for iSearchSubfieldDepth = 1:height(depthAnalysisTable)
                        
                        xStart = max(0,cell2mat(depthAnalysisTable.xStart(iSearchSubfieldDepth)))+1;
                        xEnd = xStart + cell2mat(depthAnalysisTable.xWidth(iSearchSubfieldDepth))-1;
                        yStart = max(1,cell2mat(depthAnalysisTable.yStart(iSearchSubfieldDepth)))+1;
                        yEnd = yStart + cell2mat(depthAnalysisTable.yHeight(iSearchSubfieldDepth))-1;
                        queryFeatureValue = cell2mat(depthAnalysisTable.(queryFeature)(iSearchSubfieldDepth));
                        queryFeaturePreviousValue = depthResponseMap(round((yStart + yEnd)/2), round((xStart+xEnd)/2));
                        responseTime = cell2mat(depthAnalysisTable.timePulsePeak(iSearchSubfieldDepth));
                        
                        if isnan(queryFeaturePreviousValue) || ((queryFeatureValue > queryFeaturePreviousValue) && (responseSign == 1)) || ((queryFeatureValue < queryFeaturePreviousValue) && (responseSign == -1))
                               
                               depthResponseMap(yStart:yEnd, xStart:xEnd) = queryFeatureValue;
                        
                        end
                        
                    end
                    
                    figure;
                    
                    if responseSign == 1
                        
                        colormap(cool);
                        imagesc(depthResponseMap,[0 100]);
                    
                    elseif responseSign == -1
                        
                        colormap(flipud(cool));
                        depthResponseMapDisplay = depthResponseMap;
                        depthResponseMapDisplay(isnan(depthResponseMapDisplay)) = 0;
                        imagesc(depthResponseMapDisplay,[-100 0]);
                        
                    end
                    
                    colorbar;
                    title(['Response Map - Epoch ' num2str(iEpoch) ' - Depth ' num2str(iDepth)], 'FontSize', 14);
                    
                    saveResponseMapName = ['responseMap_' responseTimeWindow '_' queryFeature '_', 'Epoch' num2str(iEpoch) '_', 'Depth' num2str(iDepth) '.mat'];
                    saveResponseMapName = strrep(saveResponseMapName, ' ', '');
                    saveResponseMapPath = fullfile(analyzedDataFolderPath, saveResponseMapName);
                    save(saveResponseMapPath, 'depthResponseMap');
                    
                    plotName = ['responseMap_' responseTimeWindow '_' queryFeature '_', 'Epoch' num2str(iEpoch) '_', 'Depth' num2str(iDepth)];
                    plotFullPath = fullfile(savePlotPath, plotName);

                    if savePlots == 1
                        
                        saveas(gcf, plotFullPath, 'pdf')
                        saveas(gcf, plotFullPath, 'fig')
                    
                    end
                    
                    overallResponseMap(~isnan(depthResponseMap)) = depthResponseMap(~isnan(depthResponseMap));

                end
                
                figure;
                                    
                if responseSign == 1
                        
                    colormap(cool);
                    imagesc(overallResponseMap,[0 100]);
                    
                elseif responseSign == -1
                        
                    colormap(flipud(cool));
                    overallResponseMapDisplay = overallResponseMap;
                    overallResponseMapDisplay(isnan(overallResponseMapDisplay)) = 0;
                    imagesc(overallResponseMapDisplay,[-100 0]);
                        
                end
                    
                colorbar;
                title(['Response Map - Epoch ' num2str(iEpoch)], 'FontSize', 14);

                saveResponseMapName = ['responseMap_' responseTimeWindow '_' queryFeature '_', 'Epoch' num2str(iEpoch) '.mat'];
                saveResponseMapName = strrep(saveResponseMapName, ' ', '');
                saveResponseMapPath = fullfile(analyzedDataFolderPath, saveResponseMapName);
                save(saveResponseMapPath, 'overallResponseMap');

                plotName = ['responseMap_' responseTimeWindow '_' queryFeature '_', 'Epoch' num2str(iEpoch)];
                plotFullPath = fullfile(savePlotPath, plotName);

                if savePlots == 1

                    saveas(gcf, plotFullPath, 'pdf')
                    saveas(gcf, plotFullPath, 'fig')
                
                end
                
                figure
                hist3([cell2mat(epochAnalysisTable.timePulsePeak),cell2mat(epochAnalysisTable.(queryFeature))],'Edges',{0:10:500 -20:1:20})
                title(['Response Distribution - Epoch ' num2str(iEpoch)], 'FontSize', 14);
                
                xlabel('Time', 'FontSize', 14)
                ylabel('Response magnitude', 'FontSize', 14)
                zlabel('Count', 'FontSize', 14)
                
                plotName = ['responseDistribution_' responseTimeWindow '_' queryFeature '_', 'Epoch' num2str(iEpoch)];
                plotFullPath = fullfile(savePlotPath, plotName);

                if savePlots == 1

                    saveas(gcf, plotFullPath, 'pdf')
                    saveas(gcf, plotFullPath, 'fig')
                
                end
                     
            end
        
        end
        
     end
     
 end