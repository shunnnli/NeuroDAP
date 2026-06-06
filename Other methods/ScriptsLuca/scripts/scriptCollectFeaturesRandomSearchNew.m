
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
            
            columnNames = {'epoch', 'depth', 'repetition', 'xStart', 'xWidth', 'yStart', 'yHeight', 'searchResponse', 'baselineAverage', 'baselineSD', 'thresholdValue', 'heightPulsePeak', 'timePulsePeak', 'areaPulse', 'isPulsePeakAboveThreshold', 'heightControlPeak', 'areaControl', 'isControlPeakAboveThreshold', 'pulseTrace', 'controlTrace', 'runProtocol'};
            dataTypes = {'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell'};
            analysisTableAllEpochs = table('Size', [0, numel(columnNames)], 'VariableNames', columnNames, 'VariableTypes', dataTypes);
            
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
                    
                    if isempty(depthSearchTable)
                        continue;
                    end
                    
                    for iRepetition = spanSearchRepetitions

                        runPath = [cellPath filesep 'ProcessedData' filesep 'Data' filesep 'Epoch' num2str(iEpoch) '_searchDepth' num2str(iDepth) '_searchRepetition' num2str(iRepetition) '.mat'];

                        if ~isfile(runPath)
                            continue;
                        end

                        disp(['Analyzing: Epoch ' num2str(iEpoch) ' - searchDepth ' num2str(iDepth) ' - searchRepetition ' num2str(iRepetition)])
                        runData = load(runPath);
                        DataStruct = runData.DataStruct; %map
                        data = DataStruct.data;
                        activeChannels = regexp(DataStruct.state.phys.internal.lastLinesUsed, '''(\w+)''', 'tokens');
                    
                        runProtocol = [];
                        runProtocol.activeChannels = activeChannels;
                        runProtocol.holdingVoltage = DataStruct.holdingVoltage;
 
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
                                halfISI = isiRed/2;

                            elseif pulsesTimeArrayBlue(1) >= pulsesTimeArrayRed(1)

                                pulsesTimeArray = pulsesTimeArrayBlue;
                                pulsesTimeDifference = pulsesTimeArrayBlue(1) - pulsesTimeArrayRed(1);
                                whichPulseFirst = 'Red';
                                halfISI = isiBlue/2;

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
                            runProtocol.pulseWidthBlue = pulseWidthBlue;
                            runProtocol.amplitudeBlue = amplitudeBlue;
                            runProtocol.delayPulseBlue = delayFirstPulseBlue;
                            runProtocol.functionNameBlue = functionNameBlue;

                            runProtocol.nPulsesRed = nPulsesRed;
                            runProtocol.pulseWidthRed = pulseWidthRed;
                            runProtocol.amplitudeRed = amplitudeRed;
                            runProtocol.delayPulseRed = delayFirstPulseRed;
                            runProtocol.functionNameRed = functionNameRed;

                            runProtocol.whichPulseFirst = whichPulseFirst;
                            runProtocol.pulsesTimeDifference = pulsesTimeDifference;
                            activeChannel = 'Both';
                          
                        elseif (any(cellfun(@(x) contains(x, 'ao1'), activeChannels)) && any(cellfun(@(x) ~contains(x, 'ao2'), activeChannels)))

                            nPulsesBlue = str2double(DataStruct.state.phys.internal.pulses.pulseString_ao1.numPulses);
                            delayFirstPulseBlue = str2double(DataStruct.state.phys.internal.pulses.pulseString_ao1.delay)*10;
                            isiBlue = str2double(DataStruct.state.phys.internal.pulses.pulseString_ao1.isi)*10;
                            pulseWidthBlue = str2double(DataStruct.state.phys.internal.pulses.pulseString_ao1.pulseWidth)*10;
                            amplitudeBlue = str2double(DataStruct.state.phys.internal.pulses.pulseString_ao1.amplitude);
                            timeArray = linspace(0,size(DataStruct.data,2)-1,size(DataStruct.data,2));
                            pulsesTimeArray = (delayFirstPulseBlue + linspace(0,(nPulsesBlue - 1)*isiBlue,nPulsesBlue));
                            pulsesWidth = pulseWidthBlue;
                            halfISI = isiBlue/2;
                            pulsesTimeDifference = nan;
                            whichPulseFirst = nan;
                            functionNameBlue = DataStruct.state.cycle.functionName;
                            functionNameBlue = strrep(functionNameBlue, '''', '');
                            functionNameBlue = strrep(functionNameBlue, ';', '');
                            functionNameBlue = strtrim(functionNameBlue);
                        
                            runProtocol.nPulsesBlue = nPulsesBlue;
                            runProtocol.pulseWidthBlue = pulseWidthBlue;
                            runProtocol.amplitudeBlue = amplitudeBlue;
                            runProtocol.delayPulseBlue = delayFirstPulseBlue;
                            runProtocol.functionNameBlue = functionNameBlue;

                            runProtocol.nPulsesRed = nan;
                            runProtocol.pulseWidthRed = nan;
                            runProtocol.amplitudeRed = nan;
                            runProtocol.delayPulseRed = nan;
                            runProtocol.functionNameRed = nan;

                            runProtocol.whichPulseFirst = whichPulseFirst;
                            runProtocol.pulsesTimeDifference = pulsesTimeDifference;
                            activeChannel = 'Blue';

                        elseif (any(cellfun(@(x) contains(x, 'ao2'), activeChannels)) && any(cellfun(@(x) ~contains(x, 'ao1'), activeChannels)))

                            nPulsesRed = str2double(DataStruct.state.phys.internal.pulses.pulseString_ao2.numPulses);
                            delayFirstPulseRed = str2double(DataStruct.state.phys.internal.pulses.pulseString_ao2.delay)*10;
                            isiRed = str2double(DataStruct.state.phys.internal.pulses.pulseString_ao2.isi)*10;
                            pulseWidthRed = str2double(DataStruct.state.phys.internal.pulses.pulseString_ao2.pulseWidth)*10;
                            amplitudeRed = str2double(DataStruct.state.phys.internal.pulses.pulseString_ao2.amplitude);                    
                            timeArray = linspace(0,size(DataStruct.data,2)-1,size(DataStruct.data,2));
                            pulsesTimeArray = (delayFirstPulseRed + linspace(0,(nPulsesRed - 1)*isiRed,nPulsesRed));
                            pulsesWidth = pulseWidthRed;
                            halfISI = isiRed/2;
                            pulsesTimeDifference = nan;
                            whichPulseFirst = nan;
                            functionNameRed = DataStruct.state.cycle.functionName;
                            functionNameRed = strrep(functionNameRed, '''', '');
                            functionNameRed = strrep(functionNameRed, ';', '');
                            functionNameRed = strtrim(functionNameRed);

                            runProtocol.nPulsesBlue = nan;
                            runProtocol.pulseWidthBlue = nan;
                            runProtocol.amplitudeBlue = nan;
                            runProtocol.delayPulseBlue = nan;
                            runProtocol.functionNameBlue = nan;

                            runProtocol.nPulsesRed = nPulsesRed;
                            runProtocol.pulseWidthRed = pulseWidthRed;
                            runProtocol.amplitudeRed = amplitudeRed;
                            runProtocol.delayPulseRed = delayFirstPulseRed;
                            runProtocol.functionNameRed = functionNameRed;

                            runProtocol.whichPulseFirst = whichPulseFirst;
                            runProtocol.pulsesTimeDifference = pulsesTimeDifference;
                            activeChannel = 'Red';
                            
                        end
                        
                        if isfield(DataStruct.state.zDMD,'options'); runProtocol.options = DataStruct.state.zDMD.options;
                        else; runProtocol.options = struct; end
                        
                        if DataStruct.holdingVoltage > -30

                            responseSign = 1;

                        else 

                            responseSign = -1;

                        end
                        
                        thresholdFactor = 3;
                        
                        % No overlap in case standard RC parameters are used
                        timeAfterEvent = 5000;
                        startFirstBaseline = min(timeAfterEvent,pulsesTimeArray(1) - 1000);
                        endFirstBaseline = pulsesTimeArray(1);
                        startSecondBaseline = min(pulsesTimeArray(end) + timeAfterEvent,length(data) - 1000);
                        endSecondBaseline = length(data);

                        % Average and SD of combined baseline windows
                        baselineAverage = mean([data(startFirstBaseline:endFirstBaseline),data(startSecondBaseline:endSecondBaseline)]);
                        baselineSD = std([data(startFirstBaseline:endFirstBaseline),data(startSecondBaseline:endSecondBaseline)]);

                        processedData = data - baselineAverage;
                        processedData = preprocessSignalVC(processedData);

                        % Initialize boxResponseData (for plotting [-10,50]ms around stim)
                        % totalLength = 100 + 500 = 600
                        boxResponseData = zeros(length(pulsesTimeArray),601);
                        boxControlData = zeros(length(pulsesTimeArray),502);
                        plotFirstSample = pulsesTimeArray - 100;
                        plotLastSample = pulsesTimeArray + 500;
                        ctrlFirstSample = pulsesTimeArray - 501;

                        for iPulse = 1:length(pulsesTimeArray)

                            pulseIdx = round(pulsesTimeArray(iPulse));
                            startIdx = pulseIdx - halfISI;
                            endIdx = pulseIdx + halfISI - 1;
                            startIdx = max(1, startIdx);
                            endIdx = min(length(data), endIdx);
                            pulseData = processedData(startIdx:endIdx);
                            peakDataStruct = scriptPerformAnalysisRandomSearchNew(pulseData, baselineSD, halfISI, responseSign, thresholdFactor);

                            pulseTrace = processedData(plotFirstSample(iPulse):plotLastSample(iPulse));
                            controlTrace = processedData(ctrlFirstSample(iPulse):pulsesTimeArray(iPulse)-1);
                        
                            pulseDataStruct = [];  
                            fieldsPeakDataStruct = fields(peakDataStruct);

                            for iFieldPeakDataStruct = 1:numel(fieldsPeakDataStruct)

                                pulseDataStruct = setfield(pulseDataStruct, fieldsPeakDataStruct{iFieldPeakDataStruct}, peakDataStruct.(fieldsPeakDataStruct{iFieldPeakDataStruct}));

                            end

                            pulseDataStruct.runProtocol = runProtocol;
                            pulseDataStructCell = struct2cell(pulseDataStruct);
                            
                            iSearchSubfieldDepth = iSearchSubfieldDepth + 1;
                            dataToAdd = [{iEpoch, iDepth, iRepetition}, table2cell(depthSearchTable(iSearchSubfieldDepth,2:end)),{baselineAverage} ,{baselineSD}, pulseDataStructCell(1:end-1)', {pulseTrace}, {controlTrace}, pulseDataStructCell(end)];
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
                    activeChannel = 'Blue';
                    % Visualiztion
                    analysisTableDepth = analysisTable(cell2mat(analysisTable.depth(:)) == iDepth, :);
                    depthResponseData = cell2mat(analysisTableDepth{:,"pulseTrace"});
                    depthControlData = cell2mat(analysisTableDepth{:,"controlTrace"});
                    
                    % Define colors used for plotting
                    %blue = [85, 161, 254]./255; % used for inhibitory
                    blue = [76, 149, 205]./255; % used for blue channel
                    %red = [255, 50, 58]./255; % used for excitatory
                    red = [217, 55, 79]./255; % used for red channel
                    purple = [232 22 224]./255; % used for -35
                    grey = [.95,.95,.95];
                    
                    if responseSign == 1; cmapBlue = createcolormap(grey,blue); cmapRed = createcolormap(grey,red);  
                    elseif responseSign == -1; cmapBlue = createcolormap(blue,grey); cmapRed = createcolormap(red,grey); end 

                    if strcmp(activeChannel,'Blue'); color = blue; 
                    elseif strcmp(activeChannel,'Red'); color = red;
                    else; color = purple; 
                    end

                    % Initialize depth related plotting params
                    plotWindowTime = linspace(-10,50,size(depthResponseData,2));

                    % Build depthCurrentMap
                    % depthCurrentMap = cell{nSpotAtDepth,1}
                    % depthCurrentMap{spot1} = zeros(nRep,nSamples)
                    depthCurrentMap = cell(4^iDepth,1);
                    spotLocationAtDepth = cell2mat(analysisTableDepth{:,["xStart", "yStart", "xWidth", "yHeight"]});
                    for spot = 1:height(spotLocationAtDepth)
                        spotTrace = depthResponseData(spot,:);

                        % Get corresponding index of the spot
                        x_index = floor(spotLocationAtDepth(spot,1)/spotLocationAtDepth(spot,3));
                        y_index = floor(spotLocationAtDepth(spot,2)/spotLocationAtDepth(spot,4));
                        spotIdx = x_index + y_index * 2^iDepth + 1;

                        % Add response to current map
                        depthCurrentMap{spotIdx} = [depthCurrentMap{spotIdx}; spotTrace];
                    end

                    % Build response map
                    depthResponseMap = zeros(684,608);
                    depthAmplitudeMap = zeros(684,608);
                    for spot = 1:height(spotLocationAtDepth)
                        % Get response value
                        yRange = spotLocationAtDepth(spot,1)+1:spotLocationAtDepth(spot,1)+spotLocationAtDepth(spot,3);
                        xRange = spotLocationAtDepth(spot,2)+1:spotLocationAtDepth(spot,2)+spotLocationAtDepth(spot,4);

                        % Get isResponse value
                        originalValue_isResponse = mode(depthResponseMap(xRange,yRange),'all');
                        newValue_isResponse = cell2mat(analysisTableDepth{spot,'isPulsePeakAboveThreshold'});

                        originalValue_amplitude = mode(depthAmplitudeMap(xRange,yRange),'all');
                        newValue_amplitude = cell2mat(analysisTableDepth{spot,'heightPulsePeak'});

                        % Add to isResponseMap
                        if originalValue_isResponse == 0
                            
                            depthResponseMap(xRange,yRange) = newValue_isResponse;
                        
                        else
                            
                            depthResponseMap(xRange,yRange) = originalValue_isResponse || newValue_isResponse;
                        
                        end
                        
                        if originalValue_amplitude == 0
                            
                            depthAmplitudeMap(xRange,yRange) = newValue_amplitude;
                            
                        elseif responseSign == 1 
                            
                            depthAmplitudeMap(xRange,yRange) = max(newValue_amplitude,originalValue_amplitude);
                           
                        elseif responseSign == -1 
                            
                            depthAmplitudeMap(xRange,yRange) = min(newValue_amplitude,originalValue_amplitude);
                                
                        end
                        
                    end

                    % Initialize plotting params
                    initializeFig(0.4,0.7);
                    masterLayout = tiledlayout(2,2);
                    masterLayout.TileSpacing = 'compact';
                    masterLayout.Padding = 'compact';

                    % Plot current map
                    nexttile(masterLayout,1,[1 1]); axis off;
                    sgtitle(['Epoch ', num2str(iEpoch), ' - Depth ', num2str(iDepth), ' - Sign ', num2str(responseSign)]);
                    title('Response traces');
                    depthLayout = tiledlayout(masterLayout,2^iDepth,2^iDepth);
                    depthLayout.Layout.Tile = 1;
                    depthLayout.Layout.TileSpan = [1 1];
                    depthLayout.TileSpacing = 'none'; depthLayout.Padding = 'tight';
                    
                    heightsPulsePeak = cell2mat(analysisTableDepth.heightPulsePeak);
                    spotsPulsePeakAboveThreshold = find(cell2mat(analysisTableDepth.isPulsePeakAboveThreshold) == 1);
                    heightsPulsePeakAboveThreshold = cell2mat(analysisTableDepth{spotsPulsePeakAboveThreshold,"heightPulsePeak"});
                    
                    heightsControlPeak = cell2mat(analysisTableDepth.heightControlPeak);
                    spotsControlPeakAboveThreshold = find(cell2mat(analysisTableDepth.isControlPeakAboveThreshold) == 1);
                    heightsControlPeakAboveThreshold = cell2mat(analysisTableDepth{spotsControlPeakAboveThreshold,"heightControlPeak"});
                    
                    maxSpotResponse = cell2mat(cellfun(@(x) max(x,[],'all'),depthCurrentMap,UniformOutput=false));
                    minSpotResponse = cell2mat(cellfun(@(x) min(x,[],'all'),depthCurrentMap,UniformOutput=false));
                      
                    absMaxCurrent = max([maxSpotResponse -minSpotResponse],[],'all');

                    for spot = 1:4^iDepth
                        nexttile(depthLayout,spot);
                        plotSEM(plotWindowTime,depthCurrentMap{spot},color,plotIndividual=true,...
                                individualColor='same',individualAlpha=0.5);
                        xlim([-10,50]); ylim([-absMaxCurrent,absMaxCurrent]);
                        plotEvent('',5,shadeOnly=true,color=color,FaceAlpha=0.3,percentY=30,zeroValue=0);

                        % Plot axis on the bottom-left tile
                        if spot ~= 4^iDepth - 2^iDepth+1; axis off;
                        else
                            xlabel('ms'); ylabel('pA'); box off;
                            xticks([0,20]); yticks([ceil(-absMaxCurrent)-eps,floor(absMaxCurrent)+eps]);
                        end

                        % Plot thresholds
                        if responseSign == 1
                            if strcmp(activeChannel,'Blue'); yline(baselineSD*thresholdFactor,'-',color=blue);
                            elseif strcmp(activeChannel,'Red'); yline(baselineSD*thresholdFactor,'-',color=red); end
                        elseif responseSign == -1
                            if strcmp(activeChannel,'Blue'); yline(-baselineSD*thresholdFactor,'-',color=blue);
                            elseif strcmp(activeChannel,'Red'); yline(-baselineSD*thresholdFactor,'-',color=red); end
                        else
                            yline(baselineSD*thresholdFactor,'-',color=blue);
                            yline(-baselineSD*thresholdFactor,'-',color=red);
                        end
                    end

                    % Plot response map
                    nexttile(masterLayout, 3);
%                     if strcmp(aiChanStr{2},'ao1')
%                         cellX = 342+round(state.zDMD.tVectorBlue(1));
%                         cellY = 304+round(state.zDMD.tVectorBlue(2));
%                     elseif strcmp(aiChanStr{2},'ao2')
%                         cellX = 342-round(state.zDMD.tVectorRed(2));
%                         cellY = 304+round(state.zDMD.tVectorRed(1));
%                     end
                    
                    cellX = 342;
                    cellY = 304;
                    imagesc(depthResponseMap); axis off; 
                    set(gca, 'clim', [0 1]); hold on;
                  
                    if strcmp(activeChannel,'Blue'); colormap(gca, [grey;blue]);
                    elseif strcmp(activeChannel,'Red'); colormap(gca, [grey;red]); end
                    
                    scatter(cellX,cellY,50,'o','filled','MarkerFaceColor','k','MarkerEdgeColor','none');
                    title('Binary hotspot map');
                    
                    % Plot response map
                    nexttile(masterLayout, 4);
%                     if strcmp(aiChanStr{2},'ao1')
%                         cellX = 342+round(state.zDMD.tVectorBlue(1));
%                         cellY = 304+round(state.zDMD.tVectorBlue(2));
%                     elseif strcmp(aiChanStr{2},'ao2')
%                         cellX = 342-round(state.zDMD.tVectorRed(2));
%                         cellY = 304+round(state.zDMD.tVectorRed(1));
%                     end
                    
                    cellX = 342;
                    cellY = 304;
                    imagesc(depthAmplitudeMap); axis off; hold on;
                    
                    if strcmp(activeChannel,'Blue'); colormap(gca, cmapBlue);
                    elseif strcmp(activeChannel,'Red'); colormap(gca, cmapRed); end
                    
                    colorbar;
                    scatter(cellX,cellY,50,'o','filled','MarkerFaceColor','k','MarkerEdgeColor','none');
                    title('Response amplitude map');
                    
                    % Plot baseline distribution vs max responses
                    nexttile(masterLayout,2);

                    s_ctrl = scatter(zeros(length(heightsControlPeak),1),heightsControlPeak,100,'filled','LineWidth',0.01);
                    s_ctrl.MarkerFaceColor = [.8,.8,.8]; s_ctrl.MarkerEdgeColor = [.8,.8,.8]; hold on;
                    
                    s_ctrlAbove = scatter(zeros(length(heightsControlPeakAboveThreshold),1),heightsControlPeakAboveThreshold,100,'filled','LineWidth',2);
                    s_ctrlAbove.MarkerFaceColor = [.8,.8,.8]; s_ctrlAbove.MarkerEdgeColor = [0,0,0]; hold on;
                    
                    p_ctrl = plot(0.2,mean(heightsControlPeak),'_','Color',[0 0 0]);
                    p_ctrl.MarkerSize = 20; p_ctrl.LineWidth = 1.2; hold on;
                    errorbar(0.2,mean(heightsControlPeak),std(heightsControlPeak),'.','Color',[0 0 0],'LineWidth', 1.2, 'CapSize', 14); hold on;

                    if strcmp(activeChannel,'Blue')
                        
                        s_pulse = scatter(ones(length(heightsPulsePeak),1),heightsPulsePeak,100,'filled','LineWidth',0.01);
                        s_pulse.MarkerFaceColor = blue; s_pulse.MarkerEdgeColor = blue; hold on;
                        s_pulseAbove = scatter(ones(length(heightsPulsePeakAboveThreshold),1),heightsPulsePeakAboveThreshold,100,'filled','LineWidth',2);
                        s_pulseAbove.MarkerFaceColor = blue; s_pulseAbove.MarkerEdgeColor = [0,0,0]; hold on;
                        
                    elseif strcmp(activeChannel,'Red')
                        
                        s_pulse = scatter(ones(length(heightsPulsePeak),1),heightsPulsePeak,100,'filled','LineWidth',0.01);
                        s_pulse.MarkerFaceColor = red; s_pulse.MarkerEdgeColor = red; hold on;
                        s_pulseAbove = scatter(ones(length(heightsPulsePeakAboveThreshold),1),heightsPulsePeakAboveThreshold,100,'filled','LineWidth',2);
                        s_pulseAbove.MarkerFaceColor = red; s_pulse.MarkerEdgeColor = [0,0,0]; hold on;
 
                    end
                    
                    p_pulse = plot(1.2,mean(heightsPulsePeak),'_','Color',[0 0 0]);
                    p_pulse.MarkerSize = 20; p_pulse.LineWidth = 1.2; hold on;
                    errorbar(1.2,mean(heightsPulsePeak),std(heightsPulsePeak),'.','Color',[0 0 0],'LineWidth', 1.2, 'CapSize', 14); hold on;

                    xlim([-0.5,1.5]);
                    xticks([0,1])
                    xticklabels({'Control','Pulse'})
                    plotStats(heightsControlPeak,heightsPulsePeak,[0,1])
                    
                    
%                     h_ctrl = histogram(cell2mat(analysisTableDepth.heightControlPeak));%,Normalization='pdf'
%                     h_ctrl.FaceColor = [.8,.8,.8]; h_ctrl.EdgeColor = [.8,.8,.8];
%                     h_ctrl.FaceAlpha = 0.5; h_ctrl.EdgeAlpha = 0.5; hold on;
%                     
                   
%                     if ~isempty(heightsControlPeakAboveThreshold)
%                         xline(heightsControlPeakAboveThreshold,'--',color=[.8,.8,.8],LineWidth=1.5); box off;
%                     end
                    
%                     % Plot thresholds
%                     if responseSign == 1
%                         xline(baselineSD*thresholdFactor,'-',color=blue,LineWidth=2,...
%                             Label=['Inhi. threshold (',num2str(thresholdFactor),'\sigma)']);
% %                         xline(baselineSD*2,'-',color=blue,Label='2\sigma',LineWidth=2);
% %                         xline(baselineSD*5,'-',color=blue,Label='5\sigma',LineWidth=2);
%                         
%                         if ~isempty(heightsPulsePeakAboveThreshold)
%                             xline(heightsPulsePeakAboveThreshold,'--',color=blue,LineWidth=1.5); box off;
%                         end
%                         
%                         h_pulse = histogram(cell2mat(analysisTableDepth.heightPulsePeak));
%                         h_pulse.FaceColor = blue; h_pulse.EdgeColor = blue;
%                         h_pulse.FaceAlpha = 0.5; h_pulse.EdgeAlpha = 0.5;
%                     elseif responseSign == -1
%                         xline(-baselineSD*thresholdFactor,'-',color=red,LineWidth=2,...
%                             Label=['Exci. threshold (',num2str(thresholdFactor),'\sigma)']);
%                         xline(-baselineSD*2,'-',color=red,Label='2\sigma',LineWidth=2);
%                         xline(-baselineSD*5,'-',color=red,Label='5\sigma',LineWidth=2);
%                         xline(minSpotResponse,'--',color=red); box off;
%                     else
%                         xline(baselineSD*thresholdFactor,'-',color=blue,LineWidth=2,...
%                             Label=['Inhi. threshold (',num2str(thresholdFactor),'\sigma)']);
%                         xline(baselineSD*2,'-',color=blue,Label='2\sigma',LineWidth=2);
%                         xline(baselineSD*5,'-',color=blue,Label='5\sigma',LineWidth=2);
%                         xline(-baselineSD*thresholdFactor,'-',color=red,LineWidth=2,...
%                             Label=['Exci. threshold (',num2str(thresholdFactor),'\sigma)']);
%                         xline(-baselineSD*2,'-',color=red,Label='2\sigma',LineWidth=2);
%                         xline(-baselineSD*5,'-',color=red,Label='5\sigma',LineWidth=2);
%                         xline(maxSpotResponse,'--',color=blue); box off;
%                         xline(minSpotResponse,'--',color=red); box off;
%                     end
                    
                    title('Response amplitude');
                    responseMapName = ['Epoch' num2str(iEpoch) '_Depth' num2str(iDepth) '_responseMap.fig'];
                    responseMapName = strrep(responseMapName, ' ', '');
                    plotFullPath = fullfile(analyzedDataFolderPath, responseMapName);
                    saveas(gcf, plotFullPath)
                
                end

                saveTableName = ['searchAnalysisTable_Epoch' num2str(iEpoch) '.mat'];
                saveTableName = strrep(saveTableName, ' ', '');
                saveTablePath = fullfile(analyzedDataFolderPath, saveTableName);
                save(saveTablePath, 'analysisTable');

                disp(['****** Epoch results saved successfully ******']);
                
                analysisTableAllEpochs = [analysisTableAllEpochs; analysisTable];
                
                runProtocolCell = struct2cell(runProtocol);
                activeChannels = runProtocol.activeChannels{1};
    
                if numel(runProtocol.activeChannels) == 2
        
                    secondActiveChannel = runProtocol.activeChannels{2};
                    activeChannels = {[activeChannels{1}, ', ', secondActiveChannel{1}]};
    
                elseif  numel(runProtocol.activeChannels) == 3
        
                    secondActiveChannel = runProtocol.activeChannels{2};
                    thirdActiveChannel = runProtocol.activeChannels{3};
                    activeChannels = {[activeChannels{1}, ', ', secondActiveChannel{1}, ', ', thirdActiveChannel{1}]};
    
                end
    
                dataCellToAdd = [{iEpoch}, {maxDepth}, activeChannels, runProtocolCell(2:end)'];
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