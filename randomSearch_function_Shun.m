function randomSearch_function(dmdNum, responseSign, selectDepths, maxDepth, maxRepeat, repeatDepthSearch, repeatLastDepthSearch)

    if nargin == 2 || selectDepths == 0
    
        maxDepth = 4; % the maximum achievable depth is 9 (2^9 = 512)
        maxRepeat = 5;
        repeatDepthSearch = 0;
        repeatLastDepthSearch = 1;
        desiredSearchDepths = [1:maxDepth];
        
    elseif selectDepths == 1
        
        maxDepth = 4; % the maximum achievable depth is 9 (2^9 = 512)
        maxRepeat = 10;
        repeatDepthSearch = 1;
        repeatLastDepthSearch = 1;
        desiredSearchDepths = [2,4];
    
    end
    
    global state

    state.epoch = state.epoch + 1;
    updateGuiByGlobal('state.epoch');
    
    nSilentBoxes = 0;
    dmdTemplate = zeros(684, 608);
    
    columnNames = {'depth', 'xStart', 'xWidth', 'yStart', 'yHeight', 'response'};
    dataTypes = {'double', 'double', 'double', 'double', 'double', 'double'};
    fullSearchTable = table('Size', [0, numel(columnNames)], 'VariableNames', columnNames, 'VariableTypes', dataTypes);
    
    columnNames = {'depth', 'xStart', 'xWidth', 'yStart', 'yHeight', 'response', 'nPositiveResponses'};
    dataTypes = {'double', 'double', 'double', 'double', 'double', 'double', 'double'};
    searchData = table('Size', [0, numel(columnNames)], 'VariableNames', columnNames, 'VariableTypes', dataTypes);

    dmdPulseNum = 200;
    aoPulseNum = 4;
    
    saveFolder = state.files.savePath;
    isi = state.cycle.zdmdPeriod1List;
    delay = state.cycle.zdmdDelay1List+round(isi/2);
    
    %% Loop through designated depth
    for depth = desiredSearchDepths
    
        nTestBoxes = determineNumberTestBoxes(searchData, depth);

        % If final depth, do repeats (default repeating for Shun)
        if (depth == maxDepth && repeatLastDepthSearch == 1)  
            repeatDepthSearch = 1;   
        end
    
        for iRepeat = 1:maxRepeat
            disp(['Depth: ' num2str(depth) ' - Repetition: ' num2str(iRepeat) ' - Running ' num2str(nTestBoxes) ' boxes'])
    
           % global iBox dmdPattern

            state.zDMD.searchDepth = depth;
            updateHeaderString('state.zDMD.searchDepth');
            state.zDMD.searchRepetition = iRepeat;
            updateHeaderString('state.zDMD.searchRepetition');
    
            state.zDMD.iBox = 1;
            state.zDMD.dmdPattern = zeros(nTestBoxes, 684, 608);
           
            % Generating the sequence of stimultion for each spot
            [searchData, fullSearchTable] = generateSearchPatterns(dmdTemplate, 0, 0, 608, 684, depth, depth, searchData, fullSearchTable, repeatDepthSearch);


            % Set all the parameters required by scanImage
            % update the pulse pattern
            state.pulses.amplitudeList(dmdPulseNum) = 5;
            state.pulses.numPulsesList(dmdPulseNum) = nTestBoxes;
            state.pulses.delayList(dmdPulseNum) = delay;
            state.pulses.isiList(dmdPulseNum) = isi;
    
            totalDuration = delay + nTestBoxes*isi + 2*delay;
            state.pulses.durationList(dmdPulseNum) = totalDuration;
            state.pulses.durationList(aoPulseNum) = totalDuration;
            
            if dmdNum == 1
              %state.cycle.zdmdImageSeries1List = 'state.zDMD.dmdPattern';
                %state.zDMD.currentSeries1 = 'state.zDMD.dmdPattern';
                %state.zDMD.currentSeries2 = [];
                state.cycle.functionName = 'randomSearchBlue';
                updateHeaderString('state.cycle.functionName');
            elseif dmdNum == 2
                %state.cycle.zdmdImageSeries2List = 'state.zDMD.dmdPattern';
                %state.zDMD.currentSeries2 = 'state.zDMD.dmdPattern';
                %state.zDMD.currentSeries1 = [];
                state.cycle.functionName = 'randomSearchRed';
                updateHeaderString('state.cycle.functionName');
            end
    
            % Writing stim params to pulse maker to make the pattern
            dmdPatternRepeats = 0;
            dmdPatternISI = 0;
            makePulsePattern(dmdPulseNum)
            changePulsePatternNumber(dmdPulseNum)
            makePulsePattern(aoPulseNum)
    
            % Start the acqusition (sweep)
            timerDoOne
    
            % b's code that doesnt matter
            if ~isempty(find(state.timer.packageStatus.*state.timer.activePackages, 1))
                disp('not done!')
            end

            % Extracting data from last acquisition
            wName = physTraceName(0, state.files.lastAcquisition);

            % Decide whether the spot is active
            activeBoxes = defineAnalysisRandomSearch(wName, responseSign, thresholdFactor);
            fullSearchTable.response((end-nTestBoxes+1):end) = activeBoxes;
            searchData.response((end-nTestBoxes+1):end) = activeBoxes;
    
            % Combine different repeats of the same spot
            searchData = compressSearchData(searchData, depth);

            % Repeat the spot if the spot is silent (not for Shun)
            if repeatDepthSearch == 0
                nSilentBoxes = numel(find(searchData.response(:) == 0 & searchData.depth(:) == depth));
                nTestBoxes = nSilentBoxes;
                if nSilentBoxes == 0; break; end
            end
        end
    
        %% Visualization 
        responseMap = zeros(684, 608);
        searchDataDepth = searchData(searchData.depth(:) == depth, :);
        nSearchBoxesDepth = size(searchDataDepth, 1);
    
        for iSearchBox = 1:nSearchBoxesDepth
            if searchDataDepth.response(iSearchBox) == 1
                responseMap((searchDataDepth.yStart(iSearchBox) + 1):(searchDataDepth.yStart(iSearchBox) + searchDataDepth.yHeight(iSearchBox)), (searchDataDepth.xStart(iSearchBox) + 1):(searchDataDepth.xStart(iSearchBox) + searchDataDepth.xWidth(iSearchBox))) = 1;
            end
        end
    

        % Initialize figure
        figure; set(gcf,'Color','w'); box off
        masterLayout = tiledlayout(2,1);
        masterLayout.TileSpacing = 'compact';
        masterLayout.Padding = 'compact';
        hAxes = gca;
        cellMap = zeros(684, 608);
        aiChanStr = state.phys.internal.lastLinesUsed;

        % Plot current trace map
        nexttile(masterLayout,1); axis off;
        title(['Depth ', num2str(depth),': current responses']); 
        depthLayout = tiledlayout(masterLayout,2^depth,2^depth);
        depthLayout.Layout.Tile = 1;
        depthLayout.TileSpacing = 'none'; depthLayout.Padding = 'tight';
        maxCurrent = max(abs(depthCurrentMap),[],'all');
        for t = 1:4^depth
            nexttile(depthLayout,t); color = [236, 73, 67];
            plotSEM(plotWindowTime,depthCurrentMap(t,plotFirstSample:plotLastSample),color,LineWidth=LineWidth);
            xlim([options.timeRange(1), options.timeRange(2)]); ylim([-maxCurrent, maxCurrent]);
            plotEvent('',5,shadeOnly=true,color=[85, 161, 254],FaceAlpha=0.3,...
                      percentY=30,zeroValue=depthCurrentMap(t,eventSample));
            % Plot axis at the bottom-left tile
            if t ~= 4^depth-2^depth+1; axis off
            else
                xlabel('ms'); ylabel('pA'); box off; 
                xticks([0,50]); yticks([ceil(-maxCurrent),floor(maxCurrent)]);
            end
        end
        
        % Plot isResponse map
        if strcmp(aiChanStr{2},'ao1')
            %cellMap = dmdPatternRotationTranslation(cellMap,1);
            cellX = [(342+round(state.zDMD.tVectorBlue(1))-10):(342+round(state.zDMD.tVectorBlue(1))+10)];
            cellY = [(304+round(state.zDMD.tVectorBlue(2))-10):(304+round(state.zDMD.tVectorBlue(2))+10)];
        elseif strcmp(aiChanStr{2},'ao2')
            %cellMap = dmdPatternRotationTranslation(cellMap,2);
            cellX = [(342-round(state.zDMD.tVectorRed(2))-10):(342-round(state.zDMD.tVectorRed(2))+10)];
            cellY = [(304+round(state.zDMD.tVectorRed(1))-10):(304+round(state.zDMD.tVectorRed(1))+10)];
        end

        cellMap(cellY,cellX) = 2;
        mapToDisplay = responseMap+cellMap;
        imagesc(hAxes,mapToDisplay);
        colormap(hAxes,[0 0 0; 1 1 1; 1 0 0])
    
        %responseMapName = ['Epoch' num2str(state.epoch) '_Depth' num2str(depth) '_responseMap.mat'];
        %responseMapPath = fullfile(saveFolder, responseMapName);
        %save(responseMapPath,'responseMap');
    
        if nSearchBoxesDepth == nSilentBoxes
            disp(['No responses found at depth: ' num2str(depth)]);
            break;
        end
    
    end
    
    %% Save
    saveFolder = state.files.savePath;
    %searchDataFileName = ['Epoch' num2str(state.epoch) '_searchData.mat'];
    %searchDataPath = fullfile(saveFolder, searchDataFileName);
    fullSearchTableName = ['Epoch' num2str(state.epoch) '_fullSearchTable.mat'];
    fullSearchTablePath = fullfile(saveFolder, fullSearchTableName);
    %save(searchDataPath,'searchData');
    save(fullSearchTablePath,'fullSearchTable');

    state.zDMD.searchDepth = [];
    updateHeaderString('state.zDMD.searchDepth');
    state.zDMD.searchRepetition = [];
    updateHeaderString('state.zDMD.searchRepetition');

end
