% clear all;
% close all;
% clc;
% 
% % User environment (automatic)
% user = getenv('USER'); 
% 
% % Indicate main directory and add it to path
% matlabDirectory = ['/Users/' user '/Desktop/Thesis/Matlab/'];
% mainDirectory = ['/Users/' user '/Desktop/Thesis/Matlab/Experiments/']; % including all Matlab toolboxes and necessary paths
% addpath(genpath(matlabDirectory));
% 
% % Indicate experiment to analyze
% experimentName = 'vGATFlp_Chrimson_1';
% 
% % Set experiment path
% experimentDirectory = [mainDirectory experimentName];
% 
% folderContent = dir([experimentDirectory]);
% folderContent = folderContent([folderContent.isdir]);
% nCells = sum([folderContent.isdir]) - 2;
% startDir = 3;
% 
% loadCells = struct;
% selectedCells = ["cell4"];
% 
% warning('off','MATLAB:unknownObjectNowStruct')

disp(['*-*-*-* Running: scriptCollectFeaturesRandomSearch *-*-*-*'])

 for iCell = startDir:(startDir+nCells-1)
     
     currentCell = folderContent(iCell);
     
     for iSelectedCell = 1:numel(selectedCells)  

        if isequal(currentCell.name, selectedCells(iSelectedCell))
            
            disp(['****** Loaded ' currentCell.name ' ******'])
            cellPath = [experimentDirectory filesep folderContent(iCell).name];
            rawDataPath = [cellPath filesep 'rawData'];
            processedDataPath = [cellPath filesep 'ProcessedData' filesep 'Data'];
            [spanEpochs, spanSearchDepths, spanSearchRepetitions] = findSpanEpochsRandomSearch(processedDataPath);
            
            columnNames = {'epoch', 'depth', 'repetition', 'xStart', 'xWidth', 'yStart', 'yHeight', 'response', 'baselineSD', 'heightPulsePeak', 'timePulsePeak', 'isPulsePeakAboveThreshold', 'areaPulse', 'heightControlPeak', 'areaControl', 'trace', 'optoParameters'};
            dataTypes = {'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell'};
%            analysisTable = table('Size', [0, numel(columnNames)], 'VariableNames', columnNames, 'VariableTypes', dataTypes);
            analysisTableAllEpochs = table('Size', [0, numel(columnNames)], 'VariableNames', columnNames, 'VariableTypes', dataTypes);
            
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
            
            columnNamesParameters = {'epoch', 'maxDepth', 'activeChannels', 'holdingVoltage', 'nPulsesBlue', 'pulseWidthBlue', 'amplitudeBlue', 'delayPulseBlue', 'functionNameBlue', 'nPulsesRed', 'pulseWidthRed', 'amplitudeRed', 'delayPulseRed', 'functionNameRed', 'whichPulseFirst', 'pulsesTimeDifference'};
            dataTypesParameters = {'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell'};
            searchParametersTable = table('Size', [0, numel(columnNamesParameters)], 'VariableNames', columnNamesParameters, 'VariableTypes', dataTypesParameters);

            for i = 1:numel(columnNamesParameters)

                searchParametersTable.(columnNamesParameters{i}) = repmat({[]}, 0, 1);

                if strcmp(dataTypes{i}, 'double')

                    searchParametersTable.(columnNamesParameters{i}) = NaN(0, 1);

                elseif strcmp(dataTypes{i}, 'logical')

                    searchParametersTable.(columnNamesParameters{i}) = false(0, 1);

                end

            end
            
            for iEpoch = spanEpochs
                
                analysisTable = table('Size', [0, numel(columnNames)], 'VariableNames', columnNames, 'VariableTypes', dataTypes);

                for i = 1:numel(columnNames)

                    analysisTable.(columnNames{i}) = repmat({[]}, 0, 1);

                    if strcmp(dataTypes{i}, 'double')

                        analysisTable.(columnNames{i}) = NaN(0, 1);

                    elseif strcmp(dataTypes{i}, 'logical')

                        analysisTable.(columnNames{i}) = false(0, 1);

                    end

                end
                
                pathFullSearchTable = [rawDataPath filesep 'Epoch' num2str(iEpoch) '_fullSearchTable.mat'];
                
                if ~isfile(pathFullSearchTable)
                    continue;
                end
                
                epochSearchTable = load(pathFullSearchTable);
                epochSearchTable = epochSearchTable.fullSearchTable;
            
                for iDepth = spanSearchDepths
               
                    depthSearchTable = epochSearchTable(epochSearchTable.depth == iDepth,:);
                    iSearchSubfieldDepth = 0;
                    
                    for iRepetition = spanSearchRepetitions

                        runPath = [cellPath filesep 'ProcessedData' filesep 'Data' filesep 'Epoch' num2str(iEpoch) '_searchDepth' num2str(iDepth) '_searchRepetition' num2str(iRepetition) '.mat'];

                        if ~isfile(runPath)
                            continue;
                        end

                        disp(['Analyzing: Epoch ' num2str(iEpoch) ' - searchDepth ' num2str(iDepth) ' - searchRepetition ' num2str(iRepetition)])
                        runData = load(runPath);
                        DataStruct = runData.DataStruct;
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
                    
                        if ~length(pulsesTimeArray) == 1

                            timeBeforePulse = (pulsesTimeArray(2) - pulsesTimeArray(1))/2;
                            timeAfterPulse = (pulsesTimeArray(2) - pulsesTimeArray(1))/2 - 1;

                        else
                            
                            timeBeforePulse = 1000;
                            timeAfterPulse = 999;
                            
                        end

                        segments = cell(1, length(pulsesTimeArray));
                        
                        preBaselineStartIdx = 4000;
                        preBaselineEndIdx = min(9000, pulsesTimeArray(1));
                        
                        postBaselineStartIdx = max(pulsesTimeArray(end)+1000, length(data)-5000);
                        postBaselineEndIdx = length(data);
                        
                        baselineAverage = mean([data(preBaselineStartIdx:preBaselineEndIdx),data(postBaselineStartIdx:postBaselineEndIdx)]);
                        baselineSD = std([data(preBaselineStartIdx:preBaselineEndIdx),data(postBaselineStartIdx:postBaselineEndIdx)]);

                        for iPulse = 1:length(pulsesTimeArray)

                            pulseIndex = round(pulsesTimeArray(iPulse));
                            startIdx = pulseIndex - timeBeforePulse;
                            endIdx = pulseIndex + timeAfterPulse;
                            startIdx = max(1, startIdx);
                            endIdx = min(length(data), endIdx);
                            
                            segments{iPulse} = data(startIdx:endIdx) - baselineAverage;
                            segments{iPulse} = preprocessSignalVC(segments{iPulse});

                        end

%                         peakDataStruct = [];
%                         peakDataStruct.heightPulsePeak = [];
%                         peakDataStruct.timePulsePeak = [];
%                         peakDataStruct.isPulsePeakAboveThreshold = [];
%                         peakDataStruct.areaPulse = [];
%                         peakDataStruct.heightControlPeak = [];
%                         peakDataStruct.areaControl = [];   
%                         peakDataStruct.trace = [];

                        if DataStruct.holdingVoltage > -30

                            responseSign = 1;

                        else 

                            responseSign = -1;

                        end

                        for iPulse = 1:length(segments)

                            iSearchSubfieldDepth = iSearchSubfieldDepth + 1;
                            peakDataStruct = scriptPerformAnalysisRandomSearch(segments{iPulse}, baselineSD, timeBeforePulse, responseSign);

%                             if ~isempty(peakData)
% 
%                                 fields = fieldnames(peakData);
% 
%                                 for iField = 1:numel(fields)
% 
%                                     entryPeakData = peakData.(fields{iField});
%                                     peakDataStruct.(fields{iField}) = [entryPeakData];
% 
%                                 end
% 
%                             end
%                             
%                             clear fields

                            pulseDataStruct = [];  
                            fieldsPeakDataStruct = fields(peakDataStruct);

                            for iFieldPeakDataStruct = 1:numel(fieldsPeakDataStruct)

                                pulseDataStruct = setfield(pulseDataStruct, fieldsPeakDataStruct{iFieldPeakDataStruct}, peakDataStruct.(fieldsPeakDataStruct{iFieldPeakDataStruct}));

                            end

                            pulseDataStruct.optoParameters = optoParameters;
                            pulseDataStructCell = struct2cell(pulseDataStruct);
                            
                            dataToAdd = [{iEpoch, iDepth, iRepetition}, table2cell(depthSearchTable(iSearchSubfieldDepth,2:end)), {baselineSD}, pulseDataStructCell'];
                            analysisTable = [analysisTable; dataToAdd];

                        end
                        
                    end

                    analyzedDataFolderName = 'Analysis';
                    analyzedDataFolderPath = [cellPath filesep analyzedDataFolderName];

                    if exist(analyzedDataFolderPath, 'dir') == 0

                        mkdir(analyzedDataFolderPath);
                        addpath(genpath(analyzedDataFolderPath));

                    end
                    
                    maxDepth = iDepth;
                
                end

                saveTableName = ['searchAnalysisTable_Epoch' num2str(iEpoch) '.mat'];
                saveTableName = strrep(saveTableName, ' ', '');
                saveTablePath = fullfile(analyzedDataFolderPath, saveTableName);
                save(saveTablePath, 'analysisTable');

                disp(['****** Epoch results saved successfully ******']);
                
                analysisTableAllEpochs = [analysisTableAllEpochs; analysisTable];
                
                optoParametersCell = struct2cell(optoParameters);
                activeChannels = optoParameters.activeChannels{1};
    
                if numel(optoParameters.activeChannels) == 2
        
                    secondActiveChannel = optoParameters.activeChannels{2};
                    activeChannels = {[activeChannels{1}, ', ', secondActiveChannel{1}]};
    
                elseif  numel(optoParameters.activeChannels) == 3
        
                    secondActiveChannel = optoParameters.activeChannels{2};
                    thirdActiveChannel = optoParameters.activeChannels{3};
                    activeChannels = {[activeChannels{1}, ', ', secondActiveChannel{1}, ', ', thirdActiveChannel{1}]};
    
                end
    
                dataCellToAdd = [{iEpoch}, {maxDepth}, activeChannels, optoParametersCell(2:end)'];
                searchParametersTable = [searchParametersTable; dataCellToAdd];
                
                clear analysisTable;
                
            end
            
            if ~isempty(analysisTableAllEpochs)

                saveTableName = 'searchAnalysisTable_AllEpochs.mat';
                saveTableName = strrep(saveTableName, ' ', '');
                saveTablePath = fullfile(analyzedDataFolderPath, saveTableName);
                save(saveTablePath, 'analysisTableAllEpochs');

                searchParametersTable = removevars(searchParametersTable, {'nPulsesBlue', 'nPulsesRed'});
                saveTableName = 'searchParametersTable_AllEpochs.mat';
                saveTableName = strrep(saveTableName, ' ', '');
                saveTablePath = fullfile(analyzedDataFolderPath, saveTableName);
                save(saveTablePath, 'searchParametersTable');

                disp(['****** Results saved successfully ******']);
                
            else 
                
                disp(['****** No random search performed ******']);
            
            end
                
        end
        
     end
     
 end
