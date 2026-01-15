function hotspotSearch_shun(dmdNum, responseSign, options)

    arguments
        dmdNum double
        responseSign double

        % Search / analysis options
        options.thresholdFactor double = 3
        options.desiredDepths double                        % candidate search depths, default is 1,2,3,4...
        options.maxSearchDepth double = 5                   % hard cap for search depth
        options.maxRepeats double = [3,3,3,5,5]             % max sweeps for each spot before declaring it a non-responsive spot
        options.repeatHotspot logical = [false, false, false, false, false]

        options.analysisWindowLength double = 50            % in ms

        % Hotspot recording options
        options.inhClamp double = 10                      % mV, absolute inhibitory clamp
        options.excClamp double = -70                       % mV, absolute excitatory clamp
        options.nFinalSweeps double = 20                   % sweeps per potential
        options.recoveryTime double = 5                    % sec, recovery time after changing clamping voltage
        options.finalDepthDiff double = 1                  % finalDepth - maxSearchDepth (how much deeper I subdivide hotspots)
        options.nSubSquaresPerHotspot double = 2           % sub-squares per hotspot when going one depth deeper

        % Distribution-matching options
        options.matchTolerance double = 0.25               % relative tolerance for mean / 25% quantile
        options.minMiniSamples double = 50                 % minimum # of minis before we trust baseline
        options.minHotspotSamples double = 10              % minimum # of hotspot responses for matching
        options.matchType string = "median"                % mean or median or kl
    end

    %% Setup params

    % Search params
    thresholdFactor   = options.thresholdFactor;
    maxSearchDepthOpt = options.maxSearchDepth;
    maxRepeats        = options.maxRepeats;
    repeatHotspot     = options.repeatHotspot;
    analysisWindowLength = options.analysisWindowLength*10;  % samples

    % Hotspot recording params
    levelDiff             = options.finalDepthDiff;
    nSubSquaresPerHotspot = options.nSubSquaresPerHotspot;
    nFinalSweeps          = options.nFinalSweeps;

    % Clamp values
    baseClamp = -35;                        % assumed amplifier holding during search

    % Depths to perform search on (upper-bounded by maxSearchDepthOpt)
    if ~isfield(options,'desiredDepths')
        desiredSearchDepths = 1:maxSearchDepthOpt;
    else
        desiredSearchDepths = options.desiredDepths;
        desiredSearchDepths = desiredSearchDepths(desiredSearchDepths <= maxSearchDepthOpt);
    end

    % Check whether maxRepeats and desiredSearchDepth have the same dim 
    if length(maxRepeats) ~= length(desiredSearchDepths)
        if isscalar(maxRepeats)
            maxRepeats = maxRepeats * ones(size(desiredSearchDepths));
        else
            error('Dimension of maxRepeats and desiredDepths should be the same!');
        end
    end

    % Check whether repeatHotspot and desiredSearchDepth have the same dim 
    if length(repeatHotspot) ~= length(desiredSearchDepths)
        error('Dimension of repeatHotspot and desiredDepths should be the same!'); 
    end

    global state

    %% Initialize data structures
    state.epoch = state.epoch + 1;
    updateGuiByGlobal('state.epoch');
    
    nSilentBoxes  = 0;
    dmdTemplate   = zeros(684, 608);

    % NOTE: Keep the canonical 6-variable layout for compatibility with
    % generateSearchPatterns() during the SEARCH phase.
    % We will append optional metadata columns (sweepStage/sweepIndex) only
    % AFTER the search is finished, so we don't break downstream table
    % concatenations inside generateSearchPatterns().
    columnNames = {'depth', 'xStart', 'xWidth', 'yStart', 'yHeight', 'response'};
    dataTypes   = {'double', 'double', 'double', 'double', 'double', 'double'};
    fullSearchTable = table('Size', [0, numel(columnNames)], ...
                            'VariableNames', columnNames, ...
                            'VariableTypes', dataTypes);
    
    columnNames = {'depth', 'xStart', 'xWidth', 'yStart', 'yHeight', 'response', 'nPositiveResponses'};
    dataTypes   = {'double', 'double', 'double', 'double', 'double', 'double', 'double'};
    searchData  = table('Size', [0, numel(columnNames)], ...
                        'VariableNames', columnNames, ...
                        'VariableTypes', dataTypes);

    dmdPulseNum = 200;  % pulse pattern for DMD (no need to change)
    aoPulseNum  = 300;  % pulse pattern for voltage command
    
    % Keep AO offset at 0 mV during the search stage (holding at -35 mV).
    state.pulses.amplitudeList(aoPulseNum) = 0;
    
    % Define save folder and fullSearchTable
    saveFolder = state.files.savePath;
    fullSearchTableName = ['Epoch' num2str(state.epoch) '_fullSearchTable.mat'];
    fullSearchTablePath = fullfile(saveFolder, fullSearchTableName);
    cleanupObj = onCleanup(@() save(fullSearchTablePath,'fullSearchTable'));

    isi   = state.cycle.zdmdPeriod1List;
    delay = state.cycle.zdmdDelay1List + round(isi/2);

    %% Running baseline mini amplitude distribution
    baselineMiniAmps = [];   % will accumulate boxControlData across all depths / repeats

    %% Variables for empirical maxSearchDepth
    foundMaxDepth  = false;
    maxSearchDepth = NaN;
    matchInfoByDepth = struct('depth', {}, 'baselineMean', {}, 'hotspotMean', {}, ...
                              'metric', {}, 'metricName', {});

    %% Perform search depth by depth
    for d = 1:length(desiredSearchDepths)
        
        % Initialize depth related params
        depth      = desiredSearchDepths(d);
        maxRepeat  = maxRepeats(d);
        nTestBoxes = determineNumberTestBoxes(searchData, depth);
        depthResponseData = []; 
        depthControlData  = [];
        repeatDepthSearch = repeatHotspot(d);
    
        % Repeat within depth
        for iRepeat = 1:maxRepeat
            disp(['Depth: ' num2str(depth) ' - Repetition: ' num2str(iRepeat) ...
                  ' - Running ' num2str(nTestBoxes) ' boxes'])
    
            state.zDMD.searchDepth     = depth;
            updateHeaderString('state.zDMD.searchDepth');
            state.zDMD.searchRepetition = iRepeat;
            updateHeaderString('state.zDMD.searchRepetition');
    
            state.zDMD.iBox       = 1;
            state.zDMD.dmdPattern = zeros(nTestBoxes, 684, 608);

            [searchData, fullSearchTable] = generateSearchPatterns( ...
                                                dmdTemplate, 0, 0, 608, 684, depth, depth, ...
                                                searchData, fullSearchTable, repeatDepthSearch);

            % Set all the parameters required by scanImage
            state.pulses.amplitudeList(dmdPulseNum) = 5;
            state.pulses.numPulsesList(dmdPulseNum) = nTestBoxes;
            state.pulses.delayList(dmdPulseNum)     = delay;
            state.pulses.isiList(dmdPulseNum)       = isi;
    
            totalDuration = delay + nTestBoxes*isi + 2*delay;
            state.pulses.durationList(dmdPulseNum) = totalDuration;
            state.pulses.durationList(aoPulseNum)  = totalDuration;

            % AO amplitude (clamp offset) is kept at 0 mV during search.
            if dmdNum == 1
                state.cycle.functionName = 'randomSearchBlue';
            elseif dmdNum == 2
                state.cycle.functionName = 'randomSearchRed';
            else
                state.cycle.functionName = 'randomSearch';
            end
            updateHeaderString('state.cycle.functionName');

            % Headerstring information
            state.zDMD.options = structToString(options);
            updateHeaderString('state.zDMD.options');
    
            dmdPatternRepeats = 0;
            dmdPatternISI     = 0;
    
            makePulsePattern(dmdPulseNum);
            changePulsePatternNumber(dmdPulseNum);
            makePulsePattern(aoPulseNum);
    
            timerDoOne;
  
            if ~isempty(find(state.timer.packageStatus .* state.timer.activePackages, 1))
                disp('not done!')
            end

            % Analyze pulse responses
            wName = physTraceName(0, state.files.lastAcquisition);
            [activeBoxes, boxResponseData, boxControlData, baselineSD, baselineMiniAmpsNew] = ...
                defineAnalysisRandomSearch(wName, responseSign, thresholdFactor, analysisWindowLength);

            % Update tables
            fullSearchTable.response((end-nTestBoxes+1):end) = activeBoxes;
            depthResponseData = [depthResponseData; boxResponseData];
            depthControlData  = [depthControlData; boxControlData];
            searchData.response((end-nTestBoxes+1):end) = activeBoxes;

            % Accumulate baseline mini amplitudes (running distribution)
            baselineMiniAmps = [baselineMiniAmps; baselineMiniAmpsNew];

            % Combine different repeats in searchData
            searchData = compressSearchData(searchData, depth);

            if ~repeatDepthSearch
                nSilentBoxes = numel(find(searchData.response(:) == 0 & searchData.depth(:) == depth));
                nTestBoxes   = nSilentBoxes;
                if nSilentBoxes == 0
                    break; 
                end
            end   
        end  % repeats

        save(fullSearchTablePath,'fullSearchTable');   % <-- SAVE PROGRESS PER DEPTH

        %% Visualization and depth-level response metrics (only if depth <= 5)
        if depth <= 5
            
            fullSearchTableDepth = fullSearchTable(fullSearchTable.depth(:) == depth, :);

            % Define colors used for plotting
            blue   = [85, 161, 254]./255;  % used for inhibitory
            red    = [255, 50, 58]./255;   % used for excitatory
            purple = [232 22 224]./255;    % used for -35
            if dmdNum == 1
                color = blue; 
            elseif dmdNum == 2
                color = red;
            else
                color = purple; 
            end

            % Initialize depth related plotting params
            aiChanStr      = state.phys.internal.lastLinesUsed;
            plotWindowTime = linspace(-10,50,size(depthResponseData,2));
            eventSample    = find(plotWindowTime == 0, 1);
            analysisWindow = eventSample : eventSample + analysisWindowLength;

            % Build depthCurrentMap
            depthCurrentMap     = cell(4^depth,1);
            spotLocationAtDepth = fullSearchTableDepth{:,["xStart", "yStart", "xWidth", "yHeight"]};
            for spot = 1:height(spotLocationAtDepth)
                spotTrace = depthResponseData(spot,:);

                % Get corresponding index of the spot
                x_index = floor(spotLocationAtDepth(spot,1)/spotLocationAtDepth(spot,3));
                y_index = floor(spotLocationAtDepth(spot,2)/spotLocationAtDepth(spot,4));
                spotIdx = x_index + y_index * 2^depth + 1;

                % Add response to current map
                depthCurrentMap{spotIdx} = [depthCurrentMap{spotIdx}; spotTrace];
            end

            % Build response map
            depthResponseMap = zeros(684,608);
            for spot = 1:height(spotLocationAtDepth)
                % Get response value
                yRange = spotLocationAtDepth(spot,1)+1 : spotLocationAtDepth(spot,1)+spotLocationAtDepth(spot,3);
                xRange = spotLocationAtDepth(spot,2)+1 : spotLocationAtDepth(spot,2)+spotLocationAtDepth(spot,4);

                % Get isResponse value
                originalValue_isResponse = mode(depthResponseMap(xRange,yRange),'all');
                newValue_isResponse      = fullSearchTableDepth{spot,'response'};
    
                % Add to isResponseMap
                if originalValue_isResponse == 0
                    depthResponseMap(xRange,yRange) = newValue_isResponse;
                else
                    depthResponseMap(xRange,yRange) = originalValue_isResponse || newValue_isResponse;
                end
            end

            % Initialize plotting params
            initializeFig(0.9,0.8);
            masterLayout              = tiledlayout(4,1);
            masterLayout.TileSpacing  = 'compact';
            masterLayout.Padding      = 'compact';

            % Find max and min current
            optoData = cell2mat(depthCurrentMap); 
            optoData = optoData(:,analysisWindow);
            optoMax  = max(optoData,[],'all'); 
            optoMin  = min(optoData,[],'all');

            % Plot current map
            nexttile(masterLayout,1,[2 1]); axis off;
            title(['Depth ', num2str(depth), ': current responses']);
            depthLayout              = tiledlayout(masterLayout,2^depth,2^depth);
            depthLayout.Layout.Tile  = 1;
            depthLayout.Layout.TileSpan = [2 1];
            depthLayout.TileSpacing  = 'none'; 
            depthLayout.Padding      = 'tight';
            maxSpotResponse = cell2mat(cellfun(@(x) max(x,[],'all'),depthCurrentMap,UniformOutput=false));
            minSpotResponse = cell2mat(cellfun(@(x) min(x,[],'all'),depthCurrentMap,UniformOutput=false));

            for spot = 1:4^depth
                nexttile(depthLayout,spot);
                plotSEM(plotWindowTime, depthCurrentMap{spot}, color, ...
                        plotIndividual=true, individualColor='same', individualAlpha=0.5);
                yaxis_limit = baselineSD * thresholdFactor * 1.5;
                xlim([-10,50]); 
                ylim([-yaxis_limit,yaxis_limit]);
                plotEvent('',5,shadeOnly=true,color=color,FaceAlpha=0.3,percentY=30,zeroValue=0);

                % Plot axis on the bottom-left tile
                if spot ~= 4^depth - 2^depth + 1
                    axis off;
                else
                    xlabel('ms'); ylabel('pA'); box off;
                    xticks([0,50]); 
                    yticks(sort(unique([round(optoMin)-eps,round(optoMax)+eps])));
                end

                % Plot thresholds
                if responseSign == 1
                    xline(baselineSD*thresholdFactor,'-',color=red);
                elseif responseSign == -1
                    xline(-baselineSD*thresholdFactor,'-',color=blue);
                else
                    xline(baselineSD*thresholdFactor,'-',color=red);
                    xline(-baselineSD*thresholdFactor,'-',color=blue);
                end
            end
        
            % Plot response map
            nexttile(masterLayout, 3);
            if strcmp(aiChanStr{2},'ao1')
                cellX = 342+round(state.zDMD.tVectorBlue(1));
                cellY = 304+round(state.zDMD.tVectorBlue(2));
            elseif strcmp(aiChanStr{2},'ao2')
                cellX = 342+round(state.zDMD.tVectorRed(1));
                cellY = 304+round(state.zDMD.tVectorRed(2));
            else
                cellX = 342;
                cellY = 304;
            end
            imagesc(depthResponseMap); axis off; 
            clim([0 1]); colormap([.95,.95,.95;color]); hold on;
            scatter(cellX,cellY,50,'o','filled','MarkerFaceColor','k','MarkerEdgeColor','none');
            title(['Depth ', num2str(depth), ': hotspot map']);
        
            % Plot baseline distribution vs max responses
            nexttile(masterLayout,4);
            h_ctrl = histogram(depthControlData,Normalization='pdf');
            h_ctrl.FaceColor = [.8,.8,.8]; 
            h_ctrl.EdgeColor = [.8,.8,.8];
            h_ctrl.FaceAlpha = 0.5; 
            h_ctrl.EdgeAlpha = 0.5;
            xline(maxSpotResponse,'--',color=red); box off;
            xline(minSpotResponse,'--',color=blue); box off;
            % Plot thresholds
            if responseSign == 1
                xline(baselineSD*thresholdFactor,'-',color=red,LineWidth=2, ...
                    Label=['Inhi. threshold (',num2str(thresholdFactor),'\sigma)']);
                xline(baselineSD*2,'-',color=red,Label='2\sigma',LineWidth=2);
                xline(baselineSD*5,'-',color=red,Label='5\sigma',LineWidth=2);
            elseif responseSign == -1
                xline(-baselineSD*thresholdFactor,'-',color=blue,LineWidth=2, ...
                    Label=['Exci. threshold (',num2str(thresholdFactor),'\sigma)']);
                xline(-baselineSD*2,'-',color=blue,Label='2\sigma',LineWidth=2);
                xline(-baselineSD*5,'-',color=blue,Label='5\sigma',LineWidth=2);
            else
                xline(baselineSD*thresholdFactor,'-',color=red,LineWidth=2, ...
                    Label=['Inhi. threshold (',num2str(thresholdFactor),'\sigma)']);
                xline(baselineSD*2,'-',color=red,Label='2\sigma',LineWidth=2);
                xline(baselineSD*5,'-',color=red,Label='5\sigma',LineWidth=2);
                xline(-baselineSD*thresholdFactor,'-',color=blue,LineWidth=2, ...
                    Label=['Exci. threshold (',num2str(thresholdFactor),'\sigma)']);
                xline(-baselineSD*2,'-',color=blue,Label='2\sigma',LineWidth=2);
                xline(-baselineSD*5,'-',color=blue,Label='5\sigma',LineWidth=2);
            end

            responseMapName = ['Epoch' num2str(state.epoch) '_Depth' num2str(depth) '_responseMap.fig'];
            responseMapPath = fullfile(saveFolder, responseMapName);
            saveas(gcf,responseMapPath);
            set(gcf,'HandleVisibility','off');   % let ScanImage plot afterwards

            %% --- NEW: compute hotspot amplitude distribution & compare to baseline minis
            % Use only responsive spots at this depth
            isHotspot = fullSearchTableDepth.response ~= 0;
            if any(isHotspot)
            
                hotspotAmpsDepth = computeHotspotAmplitudes(depthResponseData, isHotspot, responseSign, analysisWindow);
                [isMatch, matchInfo] = compareHotspotsToMini(hotspotAmpsDepth, baselineMiniAmps,...
                                                              matchTolerance=options.matchTolerance,...
                                                              minMiniSamples=options.minMiniSamples,...
                                                              minHotspotSamples=options.minHotspotSamples,...
                                                              type=options.matchType);

                if ~isempty(matchInfo)
                    matchInfo.depth = depth;
                    matchInfoByDepth(end+1) = matchInfo; 
                end

                if isMatch && ~foundMaxDepth
                    foundMaxDepth  = true;
                    maxSearchDepth = depth;
                    disp(['Empirical maxSearchDepth found at depth ' num2str(depth) '.']);
                end
            end

        end  % depth <= 5 visualization / metrics

        %% Stops search if theres no response
        searchDataDepth   = searchData(searchData.depth(:) == depth, :);
        nSearchBoxesDepth = size(searchDataDepth, 1);
        nSilentBoxes      = numel(find(searchData.response(:) == 0 & searchData.depth(:) == depth));

        if nSearchBoxesDepth == nSilentBoxes
            disp(['No responses found at depth: ' num2str(depth)]);
            maxSearchDepth = depth - 1;
            break;
        end

        % If we have already found the empirical max depth, stop going deeper
        if foundMaxDepth
            break;
        end

    end  % depth loop

    %% Decide final empirical maxSearchDepth if not found explicitly
    if ~foundMaxDepth
        responsiveDepths = unique(searchData.depth(searchData.response ~= 0));
        if ~isempty(responsiveDepths)
            maxSearchDepth = max(responsiveDepths);
        else
            maxSearchDepth = desiredSearchDepths(1);  % fallback
        end
        disp(['No distribution match found. Using maxSearchDepth = ' num2str(maxSearchDepth) '.']);
    end

    % ---------------------------------------------------------------------
    % After SEARCH is complete, append optional metadata columns so that
    % subsequent hotspot validation sweeps (stage 1 & stage 2) can be
    % represented inside fullSearchTable. This allows downstream code
    % (e.g. loadSlicesDMD.m) to use the standard "slice next N rows" logic
    % (lines ~271-285) without special-casing hotspot sweeps.
    %
    % IMPORTANT: Do NOT add these columns before the search loop, because
    % generateSearchPatterns() concatenates tables and may assume the
    % canonical 6-variable layout.
    % ---------------------------------------------------------------------
    if ~ismember('sweepStage', fullSearchTable.Properties.VariableNames)
        fullSearchTable.sweepStage = repmat("search", height(fullSearchTable), 1);
    end
    if ~ismember('sweepIndex', fullSearchTable.Properties.VariableNames)
        fullSearchTable.sweepIndex = nan(height(fullSearchTable), 1);
    end

    %% Initialize hotspotMeta for maxSearchDepth
    hotspotMeta = struct();
    hotspotMeta.epoch          = state.epoch;
    hotspotMeta.timestamp      = char(datetime('today','Format','yyyyMMdd'));
    hotspotMeta.stage          = struct('name', {}, ...
                                        'depth',{}, ...
                                        'targetVhold_mV', {}, ...
                                        'clampCommand_mV', {}, ...
                                        'nSweeps', {}, ...
                                        'acqNums', {}, ...
                                        'nBoxes', {});

    %% ---- Stage 1: record empirical hotspot depth at -70 and +10 mV ----

    % Hotspots at the empirical maxSearchDepth
    maxDepthHotspots = searchData(searchData.depth(:) == maxSearchDepth & ...
                                  searchData.response(:) ~= 0, :);

    if ~isempty(maxDepthHotspots)
        % Save for reference
        maxDepthHotspotName = ['Epoch' num2str(state.epoch) '_maxSearchHotspots_Depth' num2str(maxSearchDepth) '.mat'];
        maxDepthHotspotPath = fullfile(saveFolder, maxDepthHotspotName);
        save(maxDepthHotspotPath, 'maxDepthHotspots', 'matchInfoByDepth');

        % Change clamping to -70mV
        % Clamp command is always expressed relative to baseClamp (-35 mV).
        state.pulses.durationList(aoPulseNum)  = options.recoveryTime;
        state.pulses.amplitudeList(aoPulseNum) = options.excClamp - baseClamp;
        makePulsePattern(aoPulseNum);

        % Record at -70 mV
        disp(['Recording hotspots at depth ' num2str(maxSearchDepth) ' at -70 mV...']);
        [acqNums_exci, fullSearchTable] = runHotspotSweeps(maxDepthHotspots, nFinalSweeps, dmdPulseNum, aoPulseNum, ...
                                                          isi, delay, dmdNum, 'maxDepth_exci', ...
                                                          fullSearchTable, fullSearchTablePath);

        hotspotMeta.stage(end+1) = struct( ...
                                    'name', 'maxDepth_exci', ...
                                    'depth', maxSearchDepth, ...
                                    'targetVhold_mV', options.excClamp, ...
                                    'clampCommand_mV', options.excClamp - baseClamp, ...
                                    'nSweeps', nFinalSweeps, ...
                                    'acqNums', acqNums_exci, ...
                                    'nBoxes', height(maxDepthHotspots));

        % Change clamping to +10mV
        state.pulses.durationList(aoPulseNum)  = options.recoveryTime;
        state.pulses.amplitudeList(aoPulseNum) = options.inhClamp - baseClamp;
        makePulsePattern(aoPulseNum);

        % Record at +10 mV
        disp(['Recording hotspots at depth ' num2str(maxSearchDepth) ' at +10 mV...']);
        [acqNums_inhi, fullSearchTable] = runHotspotSweeps(maxDepthHotspots, nFinalSweeps, dmdPulseNum, aoPulseNum, ...
                                                          isi, delay, dmdNum, 'maxDepth_inhi', ...
                                                          fullSearchTable, fullSearchTablePath);

        hotspotMeta.stage(end+1) = struct( ...
                                    'name', 'maxDepth_inhi', ...
                                    'depth', maxSearchDepth, ...
                                    'targetVhold_mV', options.inhClamp, ...
                                    'clampCommand_mV', options.inhClamp - baseClamp, ...
                                    'nSweeps', nFinalSweeps, ...
                                    'acqNums', acqNums_inhi, ...
                                    'nBoxes', height(maxDepthHotspots));

        % Change back to -35mV
        state.pulses.durationList(aoPulseNum)  = options.recoveryTime;
        state.pulses.amplitudeList(aoPulseNum) = 0;
        makePulsePattern(aoPulseNum);

        % Save hotspotMeta.mat
        hotspotMetaName = ['Epoch' num2str(state.epoch) '_hotspotMeta_maxSearchDepth.mat'];
        hotspotMetaPath = fullfile(saveFolder, hotspotMetaName);
        save(hotspotMetaPath, 'hotspotMeta');

    else
        disp('No responsive hotspots at empirical maxSearchDepth, recording stopped.');
    end

    %% Initialize hotspotMeta for finalDepth
    hotspotMeta = struct();
    hotspotMeta.epoch          = state.epoch;
    hotspotMeta.timestamp      = char(datetime('today','Format','yyyyMMdd'));
    hotspotMeta.stage          = struct('name', {}, ...
                                        'depth',{}, ...
                                        'targetVhold_mV', {}, ...
                                        'clampCommand_mV', {}, ...
                                        'nSweeps', {}, ...
                                        'acqNums', {}, ...
                                        'nBoxes', {});

    %% ---- Stage 2: build next-depth hotspots and record ----

    if levelDiff > 0
        % Build final hotspot map from hotspots at maxSearchDepth by subdividing to finalDepth
        finalHotspots = maxDepthHotspots([],:);  % empty table with the same vars
        nDivPerAxis = 2^levelDiff;               % should be 2
        finalDepth = maxSearchDepth + levelDiff;
    
        for iRow = 1:height(maxDepthHotspots)
            row = maxDepthHotspots(iRow,:);
    
            x0 = row.xStart; 
            y0 = row.yStart;
            w  = row.xWidth; 
            h  = row.yHeight;
    
            subWidth  = w / nDivPerAxis;
            subHeight = h / nDivPerAxis;
    
            % Build all sub-squares for this hotspot at finalDepth
            subs = row;
            subs(1,:) = [];
            for ix = 0:(nDivPerAxis-1)
                for iy = 0:(nDivPerAxis-1)
                    newRow        = row;
                    newRow.depth  = finalDepth;
                    newRow.xStart = x0 + ix * subWidth;
                    newRow.yStart = y0 + iy * subHeight;
                    newRow.xWidth = subWidth;
                    newRow.yHeight = subHeight;
                    subs = [subs; newRow];
                end
            end
    
            % Randomly pick nSubSquaresPerHotspot sub-squares from this hotspot
            nSubs  = height(subs);
            nToPick = min(nSubSquaresPerHotspot, nSubs);
            pickIdx = randperm(nSubs, nToPick);
            finalHotspots = [finalHotspots; subs(pickIdx,:)];
        end
    
        % Remove any duplicated entries (shouldn't normally happen)
        [~, uIdx] = unique(finalHotspots{:, {'depth','xStart','yStart','xWidth','yHeight'}}, 'rows');
        finalHotspots = finalHotspots(uIdx,:);
    
        % Save final hotspot geometry
        if ~isempty(finalHotspots)
            finalHotspotName = ['Epoch' num2str(state.epoch) '_finalHotspots_Depth' num2str(finalDepth) '.mat'];
            finalHotspotPath = fullfile(saveFolder, finalHotspotName);
            save(finalHotspotPath, 'finalHotspots');
        end
    
        % Record from the final hotspot map at -70 mV and +10 mV
        if ~isempty(finalHotspots)
            % Change clamping to -70mV
            state.pulses.durationList(aoPulseNum)  = options.recoveryTime;
            state.pulses.amplitudeList(aoPulseNum) = options.excClamp - baseClamp;
            makePulsePattern(aoPulseNum);
    
            % Record at -70 mV
            disp(['Recording subdepth hotspots at depth ' num2str(finalDepth) ' at -70 mV...']);
            [acqNums_exci, fullSearchTable] = runHotspotSweeps(finalHotspots, nFinalSweeps, dmdPulseNum, aoPulseNum, ...
                                                              isi, delay, dmdNum, 'finalDepth_exci', ...
                                                              fullSearchTable, fullSearchTablePath);

            hotspotMeta.stage(end+1) = struct( ...
                                    'name', 'finalDepth_exci', ...
                                    'depth', finalDepth, ...
                                    'targetVhold_mV', options.excClamp, ...
                                    'clampCommand_mV', options.excClamp - baseClamp, ...
                                    'nSweeps', nFinalSweeps, ...
                                    'acqNums', acqNums_exci, ...
                                    'nBoxes', height(finalHotspots));
    
            % Change clamping to +10mV
            state.pulses.durationList(aoPulseNum)  = options.recoveryTime;
            state.pulses.amplitudeList(aoPulseNum) = options.inhClamp - baseClamp;
            makePulsePattern(aoPulseNum);
    
            % Record at +10 mV
            disp(['Recording subdepth hotspots at depth ' num2str(finalDepth) ' at +10 mV...']);
            [acqNums_inhi, fullSearchTable] = runHotspotSweeps(finalHotspots, nFinalSweeps, dmdPulseNum, aoPulseNum, ...
                                                              isi, delay, dmdNum, 'finalDepth_inhi', ...
                                                              fullSearchTable, fullSearchTablePath);

            hotspotMeta.stage(end+1) = struct( ...
                                    'name', 'finalDepth_inhi', ...
                                    'depth', finalDepth, ...
                                    'targetVhold_mV', options.inhClamp, ...
                                    'clampCommand_mV', options.inhClamp - baseClamp, ...
                                    'nSweeps', nFinalSweeps, ...
                                    'acqNums', acqNums_inhi, ...
                                    'nBoxes', height(finalHotspots));
    
            % Change back to -35mV
            state.pulses.durationList(aoPulseNum)  = options.recoveryTime;
            state.pulses.amplitudeList(aoPulseNum) = 0;
            makePulsePattern(aoPulseNum);

            % Save hotspotMeta.mat
            hotspotMetaName = ['Epoch' num2str(state.epoch) '_hotspotMeta_finalDepth.mat'];
            hotspotMetaPath = fullfile(saveFolder, hotspotMetaName);
            save(hotspotMetaPath, 'hotspotMeta');
        else
            disp('No final hotspots; skipping second hotspot recordings.');
        end
    end

    %% Save data
    save(fullSearchTablePath,'fullSearchTable');
    state.zDMD.searchDepth = [];
    updateHeaderString('state.zDMD.searchDepth');
    state.zDMD.searchRepetition = [];
    updateHeaderString('state.zDMD.searchRepetition');
end



%% --- Helper: sweep a set of DMD boxes with current clamp offset ---
function [acqNums, fullSearchTable] = runHotspotSweeps(finalHotspots, nSweeps, dmdPulseNum, aoPulseNum, ...
                          isi, delay, dmdNum, stageName, fullSearchTable, fullSearchTablePath)

    if isempty(finalHotspots)
        acqNums = [];
        return;
    end

    global state

    % Setup hotspot meta params
    state.zDMD.sweepStage = stageName;
    acqNums = nan(1, nSweeps);

    sweepDepth = finalHotspots.depth(1);  % all boxes share the same depth

    for iSweep = 1:nSweeps
        nBoxes = height(finalHotspots);
        disp(['Hotspot sweep ' num2str(iSweep) '/' num2str(nSweeps) ...
              ', depth ' num2str(sweepDepth) ', ' num2str(nBoxes) ' boxes']);

        state.zDMD.searchDepth     = sweepDepth;
        updateHeaderString('state.zDMD.searchDepth');
        state.zDMD.searchRepetition = iSweep;
        updateHeaderString('state.zDMD.searchRepetition');

        state.zDMD.iBox       = 1;
        state.zDMD.dmdPattern = zeros(nBoxes, 684, 608);

        % -------------------------------------------------------------
        % FULLSEARCHTABLE APPEND (stage 1 / stage 2 hotspot sweeps)
        %
        % Append these sweep boxes to fullSearchTable in the exact order
        % they are presented as DMD pulses. This lets downstream analysis
        % (e.g. loadSlicesDMD) reuse the standard slicing logic.
        %
        % NOTE: hotspot sweeps typically don't have online response labels,
        % so we store response = NaN here.
        % -------------------------------------------------------------
        sweepRows = finalHotspots(:, {'depth','xStart','xWidth','yStart','yHeight'});
        sweepRows.response   = nan(nBoxes, 1);
        sweepRows.sweepStage = repmat(string(stageName), nBoxes, 1);
        sweepRows.sweepIndex = repmat(iSweep, nBoxes, 1);
        fullSearchTable = [fullSearchTable; sweepRows(:, fullSearchTable.Properties.VariableNames)];
        if nargin >= 10 && ~isempty(fullSearchTablePath)
            save(fullSearchTablePath,'fullSearchTable');
        end

        % Build DMD patterns for each hotspot square
        for iBox = 1:nBoxes
            xStart  = finalHotspots.xStart(iBox);
            yStart  = finalHotspots.yStart(iBox);
            xWidth  = finalHotspots.xWidth(iBox);
            yHeight = finalHotspots.yHeight(iBox);
            pattern = zeros(684, 608);

            % Convert DMD coordinates (0-based) + widths to integer pixel indices
            rowStart = floor(yStart) + 1;
            rowEnd   = floor(yStart + yHeight);
            colStart = floor(xStart) + 1;
            colEnd   = floor(xStart + xWidth);
            
            % Clamp to the DMD image bounds: 684 x 608
            rowStart = max(1, min(684, rowStart));
            rowEnd   = max(1, min(684, rowEnd));
            colStart = max(1, min(608, colStart));
            colEnd   = max(1, min(608, colEnd));
            
            if rowEnd >= rowStart && colEnd >= colStart
                rows = rowStart:rowEnd;
                cols = colStart:colEnd;
                pattern(rows, cols) = 1;
            else
                warning('Hotspot box %d produced an empty index range; skipping.', iBox);
            end
            state.zDMD.dmdPattern(iBox,:,:) = pattern;
        end

        % Configure DMD pulse pattern
        state.pulses.amplitudeList(dmdPulseNum) = 5;
        state.pulses.numPulsesList(dmdPulseNum) = nBoxes;
        state.pulses.delayList(dmdPulseNum)     = delay;
        state.pulses.isiList(dmdPulseNum)       = isi;

        totalDuration = delay + nBoxes*isi + 2*delay;
        state.pulses.durationList(dmdPulseNum) = totalDuration;
        state.pulses.durationList(aoPulseNum)  = totalDuration;

        % Keep the same functionName convention as in the search
        if dmdNum == 1
            state.cycle.functionName = 'randomSearchBlue';
        elseif dmdNum == 2
            state.cycle.functionName = 'randomSearchRed';
        else
            state.cycle.functionName = 'randomSearch';
        end
        updateHeaderString('state.cycle.functionName');

        makePulsePattern(dmdPulseNum);
        changePulsePatternNumber(dmdPulseNum);
        makePulsePattern(aoPulseNum);

        timerDoOne;
        acqNums(iSweep) = getLastAcqNumber();

        if ~isempty(find(state.timer.packageStatus .* state.timer.activePackages, 1))
            disp('not done!');
        end
    end
end

%% --- Helper: compute hotspot amplitudes for a given depth ---
function hotspotAmps = computeHotspotAmplitudes(depthResponseData, isHotspot, ...
                                               responseSign, analysisWindow)
    % depthResponseData: nBoxes x nSamples
    % isHotspot: logical vector over nBoxes (true = hotspot)
    % analysisWindow: sample indices (e.g. 50 ms after stim)

    depthResponseData = depthResponseData(isHotspot, :);

    if isempty(depthResponseData)
        hotspotAmps = [];
        return;
    end

    windowTraces = depthResponseData(:, analysisWindow);

    switch sign(responseSign)
        case 1   % inhibitory (positive-going)
            % Use absolute value of minimum in the window
            hotspotAmps = max(windowTraces, [], 2);
        case -1  % excitatory (negative-going)
            hotspotAmps = -min(windowTraces, [], 2);
        otherwise % responseSign == 0, use rectified area over the window (~50 ms)
            hotspotAmps = sum(abs(windowTraces), 2);
    end
end

%% --- Helper: compare hotspot distribution to baseline minis ---
function [isMatch, info] = compareHotspotsToMini(hotspotAmps, baselineMiniAmps, opts)
    % hotspotAmps, baselineMiniAmps: column or row vectors
    % opts.matchTolerance: scalar
    % opts.minMiniSamples: minimum # of baseline samples
    % opts.minHotspotSamples: minimum # of hotspot samples
    % opts.type: 'mean' or 'kl' or 'median'
    %
    % isMatch: true if distributions are "close" under chosen metric
    % info: struct with baselineMean, hotspotMean, metric, metricName

    arguments
        hotspotAmps (:,1) double
        baselineMiniAmps (:,1) double
        opts.matchTolerance double
        opts.minMiniSamples double
        opts.minHotspotSamples double
        opts.type char {mustBeMember(opts.type,{'mean','kl','median'})}
    end

    isMatch = false;
    info    = struct('depth', [], 'baselineMean', [], ...
                     'hotspotMean', [], 'metric', [], ...
                     'metricName', '');

    % Ensure column vectors
    hotspotAmps     = hotspotAmps(:);
    baselineMiniAmps = baselineMiniAmps(:);

    % Not enough samples? -> no decision
    if numel(baselineMiniAmps) < opts.minMiniSamples || ...
       numel(hotspotAmps)     < opts.minHotspotSamples
        return;
    end

    % Common stats we may want to record
    baselineMean = mean(baselineMiniAmps);
    hotspotMean  = mean(hotspotAmps);
    baselineMid  = median(baselineMiniAmps);
    hotspotMid   = median(hotspotAmps);

    info.baselineMean = baselineMean;
    info.hotspotMean  = hotspotMean;

    switch opts.type
        case 'mean'
            % Relative difference of means
            metric = abs(hotspotMean - baselineMean) / max(abs(baselineMean), eps);
            info.metric     = metric;
            info.metricName = 'relMeanDiff';
            isMatch = metric <= opts.matchTolerance;
           
        case 'median'
            % Relative difference of medians
            metric = abs(hotspotMid - baselineMid) / max(abs(baselineMid), eps);
            info.metric     = metric;
            info.metricName = 'relMedianDiff';
            isMatch = metric <= opts.matchTolerance;

        case 'kl'
            % Symmetric KL divergence between histograms
            allData = [baselineMiniAmps; hotspotAmps];
            dataMin = min(allData);
            dataMax = max(allData);

            if dataMax == dataMin
                % Degenerate case: both are essentially delta at same value
                metric = 0;
            else
                nBins = 20;  % you can tune this
                edges = linspace(dataMin, dataMax, nBins+1);

                pCounts = histcounts(hotspotAmps,      edges);
                qCounts = histcounts(baselineMiniAmps, edges);

                % Add small epsilon for numerical stability
                epsProb = 1e-12;
                p = pCounts + epsProb;
                q = qCounts + epsProb;
                p = p / sum(p);
                q = q / sum(q);

                % KL(p||q) and KL(q||p)
                kl_pq = sum(p .* log(p ./ q));
                kl_qp = sum(q .* log(q ./ p));

                % Symmetric KL
                metric = 0.5 * (kl_pq + kl_qp);
            end

            info.metric     = metric;
            info.metricName = 'symKL';
            isMatch = metric <= opts.matchTolerance;
    end
end


%% --- Helper: calculate baseline statistics ---------
function [activeBoxes, boxResponseData, boxControlData, baselineSD, baselineMiniAmps] = ...
          defineAnalysisRandomSearch(wName, responseSign, thresholdFactor, analysisWindowLength)

    global state

    waveData = getWave(wName, 'data');
    hString  = getfield(get(wName, 'UserData'), 'headerString');
    aiChanStr = state.phys.internal.lastLinesUsed;

    pulseString     = valueFromHeaderString(['state.phys.internal.pulseString_' aiChanStr{2}], hString);
    delayFirstPulse = phUtil_parsePulsePatternString(pulseString, 'delay');
    isi             = phUtil_parsePulsePatternString(pulseString, 'isi');
    nPulses         = phUtil_parsePulsePatternString(pulseString, 'numPulses');
    halfISI         = 10*isi/2;

    % times in samples (10 kHz assumed)
    pulsesTimeArray = (delayFirstPulse + linspace(0,(nPulses - 1)*isi,nPulses))*10;

    % activeBoxes will be 0/1 flags per pulse
    activeBoxes = zeros(1, length(pulsesTimeArray));

    % ---- Baseline windows used for mean / SD ----
    timeAfterEvent       = 5000;
    startFirstBaseline   = min(timeAfterEvent, pulsesTimeArray(1) - 1000);
    endFirstBaseline     = pulsesTimeArray(1);
    startSecondBaseline  = min(pulsesTimeArray(end) + timeAfterEvent, length(waveData) - 1000);
    endSecondBaseline    = length(waveData);
    baselineWaveData     = [waveData(startFirstBaseline:endFirstBaseline), waveData(startSecondBaseline:endSecondBaseline)];
    
    % Average and SD of combined baseline windows
    baselineAverage = mean(baselineWaveData);
    baselineSD = std(baselineWaveData);

    % Get baseline mini distribution
    baselineMiniAmps = getMinis(baselineWaveData, baselineSD, responseSign, thresholdFactor);

    % Subtract mean and preprocess
    processedWaveData = waveData - baselineAverage;
    processedWaveData = preprocessSignalVC(processedWaveData);

    % ---- Preallocate response & control matrices ----
    % response: [-10, +50] ms around stim (assuming 10 kHz -> 601 samples)
    boxResponseData = zeros(length(pulsesTimeArray), 601);
    % baseline: ~50 ms before stim (502 samples)
    boxControlData  = zeros(length(pulsesTimeArray), 502);
    % baseline mini amplitudes
    baselineMiniAmps   = zeros(length(pulsesTimeArray), 1);

    plotFirstSample = pulsesTimeArray - 100;       % -10 ms
    plotLastSample  = pulsesTimeArray + 500;       % +50 ms
    ctrlFirstSample = pulsesTimeArray - 502;       % ~-50 ms

    for iPulse = 1:length(pulsesTimeArray)

        pulseIdx = round(pulsesTimeArray(iPulse));
        startIdx = pulseIdx - halfISI;
        endIdx   = pulseIdx + halfISI - 1;
        startIdx = max(1, startIdx);
        endIdx   = min(length(waveData), endIdx);
        pulseWaveData = processedWaveData(startIdx:endIdx);

        % Main hotspot/response classification
        boxResponse            = performAnalysisRandomSearch(pulseWaveData, baselineSD, halfISI, responseSign, thresholdFactor, analysisWindowLength);
        activeBoxes(iPulse)    = boxResponse;

        % Store response snippet [-10, +50] ms
        boxResponseData(iPulse,:) = processedWaveData(plotFirstSample(iPulse):plotLastSample(iPulse));

        % Store baseline snippet (~50 ms before stim)
        boxControlData(iPulse,:)  = processedWaveData(ctrlFirstSample(iPulse):pulsesTimeArray(iPulse)-1);

    end
end


%% --- Helper: calculate minis amplitude from some trace ---
function baselineMinisAmp = getMinis(baselineWaveData, baselineSD, responseSign, thresholdFactor)
    % If responseSign == 0; use rect area 50ms around the peak as amp
    % else, use peak amplitude
    
    thresholdPeak = thresholdFactor * baselineSD;
    
    fs = 10e3;                 % sampling rate (Hz)
    winMs = 50;                % window length in ms
    winSamples = round(fs * winMs / 1000);   % 50 ms -> 500 samples
    halfWin   = floor(winSamples/2);
    
    if responseSign == 1
        % Positive-going response
        [baselineMinisAmp, ~] = findpeaks(baselineWaveData, 'MinPeakHeight', thresholdPeak);
    elseif responseSign == -1
        % Negative-going response
        [baselineMinisAmp, ~] = findpeaks(-baselineWaveData, 'MinPeakHeight', thresholdPeak);
    else
        % Find all peaks
        [peaks, peakLocs] = findpeaks(abs(baselineWaveData), 'MinPeakHeight', thresholdPeak);
        
        nSamples  = numel(baselineWaveData);
        baselineMinisAmp = nan(size(peaks));   % rectified area for each peak
        
        for k = 1:numel(peakLocs)
            center = peakLocs(k);
        
            % Initial symmetric window around the peak
            winStart = center - halfWin;
            winEnd   = center + halfWin - 1;   % so total length ~ winSamples
        
            % Clamp to signal bounds, keeping window length as close as possible
            if winStart < 1
                winStart = 1;
                winEnd   = min(nSamples, winStart + winSamples - 1);
            elseif winEnd > nSamples
                winEnd   = nSamples;
                winStart = max(1, winEnd - winSamples + 1);
            end
        
            seg = baselineWaveData(winStart:winEnd);
            baselineMinisAmp(k) = sum(abs(seg));   % rectified area for this peak
        end
    end
end


function acq = getLastAcqNumber()
    % This is intentionally defensive because ScanImage state fields vary.
    global state
    acq = NaN;

    % Try common field names used in ScanImage variants
    candidates = { ...
        'state.files.lastAcquisition', ...
        'state.files.lastAcqNumber', ...
        'state.files.internal.lastAcquisition', ...
        'state.files.internal.lastAcqNumber', ...
        'state.internal.lastAcquisition', ...
        'state.internal.lastAcqNumber'};

    for k = 1:numel(candidates)
        try
            val = eval(candidates{k});
            if isnumeric(val) && isscalar(val) && ~isnan(val)
                acq = val;
                return;
            end
        catch
            % ignore missing fields
        end
    end
end