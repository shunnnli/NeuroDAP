function outVals = laserShuffleAnalysis(fullFileName, paramVals, selectedCellNames)


    preWindowSize = paramVals.preLaserWindow(2) - paramVals.preLaserWindow(1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                   READ DATA IN FROM .CSV FILE                           %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    sdata = readLaserShuffleData(fullFileName, paramVals.laserColName);
   
    % get cell columns for specified cell names
    cellCol = [];
    cellNames = {};
    for k = 1:length(selectedCellNames)       
        cname = selectedCellNames{k};
        if sdata.cell2column.isKey(cname)
            cellCol(end+1) = sdata.cell2column(cname);
            cellNames{end+1} = cname;
        end
    end
        
    numCells = length(cellCol);
    allCellNums = 1:numCells;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %               MAKE THE FIGURE FOR PLOTTING THE DATA                     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Make the output figure for plotting cell data
    plotFig = figure('color','w');
    rasterAx = subplot(4,1,1,'parent',plotFig);
    psthAx = subplot(4,1,2,'parent',plotFig);
    plotAx = subplot(4,1,3,'parent',plotFig);
    summaryAx = subplot(4,1,4,'parent',plotFig);
    hold(plotAx,'on');
    hold(psthAx,'on');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %               LOOP THROUGH ALL CELLS TO ANALYZE                         %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Cycle through each cell in the data set
    for cellNum = allCellNums
    
        watchon;
        display(sprintf('Processing cell %s...', cellNames{cellNum}));
        
        % Create empty arrays to hold PSTH and randomized/shifted PSTH values for this cell
        psthVals{cellNum} = [];
        psthShiftList{cellNum} = [];

        laserPulseNums = 1:size(sdata.laserTimes, 1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %            FIND AND SHUFFLE BASELINE SPIKE TRAINS                   %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Cycle through all laser pulses at this laser power to pull spikes
        % from the baseline period and make shuffled lists
        for laserPulseNum = laserPulseNums
            % Find the spikes from this cell that fall into the pre-stim
            % window for this laser pulse
            befSpikeNums{laserPulseNum} = find((sdata.allData.data(:, cellCol(cellNum)) > ...
                (sdata.laserTimes(laserPulseNum) + paramVals.preLaserWindow(1))) & ...
                (sdata.allData.data(:, cellCol(cellNum)) < ...
                (sdata.laserTimes(laserPulseNum) + paramVals.preLaserWindow(2))));
            % Find all spikes from this cell that fall into the PSTH window
            % for this laser pulse
            nearPulseNums{laserPulseNum} = find((sdata.allData.data(:, cellCol(cellNum)) > ...
                (sdata.laserTimes(laserPulseNum) + paramVals.psthWindow(1))) & ...
                (sdata.allData.data(:, cellCol(cellNum)) < ...
                (sdata.laserTimes(laserPulseNum) + paramVals.psthWindow(2))));
            % If there are any spikes in the pre-stim window, then shuffle
            % them and add to the shuffled PSTH list
            if any(befSpikeNums{laserPulseNum})
                % Make a list of pre-stim spike times relative to the laser pulse
                subTimes = sdata.allData.data(befSpikeNums{laserPulseNum}, cellCol(cellNum)) - ...
                    (sdata.laserTimes(laserPulseNum) + paramVals.preLaserWindow(1));
                % Replicate the list of pre-stim spikes according to the
                % number of randomized shifts (paramVals.numRandShifts)
                repSubTimes = (diag(subTimes) * eye(size(subTimes,2))) * ...
                    ones(size(subTimes, 1), paramVals.numRandShifts);
                % Generate a randomized list of spike offsets
                shiftVals = diag((paramVals.randShiftMultiplier .* ...
                    rand(1, paramVals.numRandShifts))) * ...
                    ones(paramVals.numRandShifts, size(subTimes, 1));
                % Add the randomized list of spike offsets to the list of
                % pre-stim spikes and use the mod function to "wrap around"
                % spikes that fell off of the edge of the window
                shiftSubTimes{laserPulseNum} = ...
                    mod(repSubTimes + shiftVals', preWindowSize) + ...
                    paramVals.preLaserWindow(1);
                % Add the spikes from this laser pulse to the full list of
                % shifted spike times
                psthShiftList{cellNum} = cat(1, psthShiftList{cellNum}, ...
                    shiftSubTimes{laserPulseNum});
            end
            % If any spikes fell within the PSTH window, add them to the
            % full list of PSTH spikes. Final list is 2-d array. First 
            % column is the plaser pulse number, second column is the spike
            % time relative to that pulse. 
            if any(nearPulseNums{laserPulseNum})
                % Generate list of spike times, filled initally with laser
                % pulse number
                addVals = laserPulseNum .* ones(size(nearPulseNums{laserPulseNum}, 1), 2);
                % Then swap out second column of list to contain spike
                % times relative to the paser pulse
                addVals(:, 2) = sdata.allData.data(nearPulseNums{laserPulseNum}, ...
                    cellCol(cellNum)) - sdata.laserTimes(laserPulseNum);
                % Add to the full list of PSTH vals
                psthVals{cellNum} = cat(1, psthVals{cellNum}, addVals);
            end
        end 

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %            MAKE HISTOGRAM OF SHUFFLED HISTOGRAM BINS                %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Find the PSTH bins corresponding to the pre-stim window and
        % calculate the total number of shuffled PSTH bins for this cell
        histShuffleBins = find(paramVals.psthBins > paramVals.preLaserWindow(1) & ...
            paramVals.psthBins < paramVals.preLaserWindow(2));
        numShuffleHistVals = paramVals.numRandShifts * (size(histShuffleBins,2) - 1);
        shuffleHistVals = zeros(paramVals.numRandShifts, (size(histShuffleBins, 2)));
        numLaserPSTHBins = size(paramVals.laserPSTHBins, 2);

        % Calculate a PSTH for each randomly shifted set of spike times
        if any(psthShiftList{cellNum})
            for randNum = 1:paramVals.numRandShifts
                shuffleHistVals(randNum, :) = ...
                    histc(psthShiftList{cellNum}(:, randNum), paramVals.psthBins(histShuffleBins));
            end
        else
            shuffleHistVals = zeros(paramVals.numRandShifts, ...
                size(paramVals.psthBins(histShuffleBins),2));
        end

        % Drop the last value because it is always zero and reshape into a 1-d array
        shuffleHistVals(:, end) = [];
        shuffleHistValList = reshape(shuffleHistVals, numShuffleHistVals, 1);

        % Calculate a histogram of the PSTH values
        shuffleHistHist = histc(shuffleHistValList, ...
            paramVals.histHistBins) ./ paramVals.numRandShifts;
        [maxVal, maxBin] = max(shuffleHistHist);
        lastBin = find(shuffleHistHist > 0, 1, 'last');


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %    SCORE BINS FROM LASER PULSE RELATIVE TO SHUFFLED BASELINE        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        allDropPulseNums = randperm(size(sdata.laserTimes, 1));
        dropPulsesPerFrac = ceil(size(sdata.laserTimes, 1) / paramVals.numDropFractions);

        % Do the scoring repeatedly by dropping a fraction of the laser pulses
        % and re-calculating the laser triggered histogram. This way a single
        % burst of activity in a low-firing cell does not get counted as a
        % laser response
        for dropFracNum = 1:paramVals.numDropFractions
            firstPulse = 1 + ((dropFracNum - 1) * dropPulsesPerFrac);
            lastPulse = dropFracNum * dropPulsesPerFrac;
            if lastPulse > size(allDropPulseNums, 2)
                lastPulse =  size(allDropPulseNums, 2);
            end
            dropPulses = allDropPulseNums(firstPulse:lastPulse);
            scalePSTH = 1 + (size(dropPulses, 2) / size(sdata.laserTimes, 1));

            % Calculate the PSTH for this cell multiple times by dropping a
            % small fraction of pulses each time
            if any(psthVals{cellNum})
                psthLinesToDrop = [];
                for dropPulseNum = dropPulses
                    psthLinesToDrop = cat(1, psthLinesToDrop, ...
                        find(psthVals{cellNum}(:,1) == dropPulseNum));
                end
                remainPSTHLines = 1:size(psthVals{cellNum},1);
                if any(psthLinesToDrop)
                    remainPSTHLines(psthLinesToDrop) = [];
                end
                dropPSTHBarVals(cellNum, :) = histc(psthVals{cellNum}(remainPSTHLines, 2), ...
                    paramVals.psthBins) .* scalePSTH;
            else
                dropPSTHBarVals(cellNum, :) = zeros(1,size(paramVals.psthBins,2));
            end
            % Run through PSTH bins in laser pulse and calculate percentile score
            % relative to shuffled PSTH bins of baseline firing
            for laserBinNum = 1:numLaserPSTHBins
                % Find the bin of the histogram of PSTH values that corresponds
                % to the PSTH at the counter laserBinNum
                if dropPSTHBarVals(cellNum, paramVals.laserPSTHBins(laserBinNum)) <= ...
                        paramVals.histHistBins(end)
                    histHistLaserBinNum(laserBinNum) = find(paramVals.histHistBins >= ...
                        dropPSTHBarVals(cellNum, paramVals.laserPSTHBins(laserBinNum)), 1, 'first');
                else
                    disp('*****************************************************');
                    disp('* ERROR: Max value in histHistBin array is too low! *');
                    disp('*****************************************************');
                    histHistLaserBinNum(laserBinNum) = paramVals.histHistBins(end);
                end
                % Calculate a percentile score for the likelihood that this
                % PSTH value is significantly HIGHER than the baseline firing
                laserMaxPercentile(laserBinNum) = ...
                    (sum(shuffleHistHist(1:histHistLaserBinNum(laserBinNum))) - ...
                    shuffleHistHist(histHistLaserBinNum(laserBinNum))) / ...
                    sum(shuffleHistHist);
                % Calculate a percentile score for the likelihood that this
                % PSTH value is significantly LOWER than the baseline firing
                laserMinPercentile(laserBinNum) = ...
                    (sum(shuffleHistHist(histHistLaserBinNum(laserBinNum):end)) - ...
                    shuffleHistHist(histHistLaserBinNum(laserBinNum))) / ...
                    sum(shuffleHistHist);
                % Test whether percentile scores exceed the cutoffs
                laserBinMaxDropSig(dropFracNum, laserBinNum) = ...
                    ((laserMaxPercentile(laserBinNum) >= paramVals.percentileCutOff)) && ...
                    ((dropPSTHBarVals(cellNum, paramVals.laserPSTHBins(laserBinNum)) / ...
                    scalePSTH) > paramVals.minSpikesPerSigBin);
                laserBinMinDropSig(dropFracNum, laserBinNum) = ...
                    (laserMinPercentile(laserBinNum) >= paramVals.percentileCutOff);
            end
        end
        if any(psthVals{cellNum})
            psthBarVals(cellNum, :) = histc(psthVals{cellNum}(:, 2), paramVals.psthBins);
        else
            psthBarVals(cellNum, :) = zeros(1,size(paramVals.psthBins,2));
        end
        laserBinMaxSig(cellNum, :) = min(laserBinMaxDropSig, [], 1);
        laserBinMinSig(cellNum, :) = min(laserBinMinDropSig, [], 1);

        % If not enough bins exceeded significance from this cell, set all 
        % significance values to zero for this cell
        if sum(laserBinMaxSig(cellNum, :)) < paramVals.minSigBins
            laserBinMaxSig(cellNum, :) = 0;
        end

        laserBinTimes = paramVals.psthBins(paramVals.laserPSTHBins);
        laserBinTimesPosBins = find(laserBinTimes >= 0);

        if any(laserBinMaxSig(cellNum, laserBinTimesPosBins) > 0)
            firstSigMaxBin(cellNum) = laserBinTimes(laserBinTimesPosBins(...
                find(laserBinMaxSig(cellNum, laserBinTimesPosBins) > 0, 1, 'first')));
        else
            firstSigMaxBin(cellNum) = -1;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                     PLOT DATA TO FIGURE                             %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Clear axes and generate the raster plot of spike times
        cla(rasterAx);
        cla(psthAx);
        cla(plotAx);
        drawnow;
        if any(psthVals{cellNum})
            scatter(psthVals{cellNum}(:,2), psthVals{cellNum}(:,1),'marker','.', ...
                'markerfacecolor', 'k', 'markeredgecolor', 'k', 'parent',rasterAx);
            set(rasterAx,'ydir', 'reverse');
        else
            cla(rasterAx);
        end
        title(rasterAx, ['Cell ', num2str(cellNum), ' of ', ...
            num2str(numCells), ': ', cellNames{cellNum}], 'interpreter', 'none');
        ylabel(rasterAx, 'Cell #');

        % Plot the PSTH and the histogram of shuffled PSTH values
        title(plotAx,['Max: ', num2str(laserMaxPercentile(1)), ...
            '; Min: ', num2str(laserMinPercentile(1))]);
        normPsthBarVals(cellNum, :) = psthBarVals(cellNum, :) ./ ...
            (size(sdata.laserTimes, 1) * paramVals.psthBinSize);
        histHistSum = sum(shuffleHistHist(1:end));
        bar(paramVals.psthBins + (paramVals.psthBinSize / 2), ...
            normPsthBarVals(cellNum, :), ...
            'edgecolor', 'k', 'facecolor', 'k', 'parent',psthAx);
        ylabel(psthAx, 'Rate (Hz)');

        % Plot histogram of shuffled PSTH bin values
        bar(paramVals.histHistBins(1:end) ./ (size(sdata.laserTimes, 1) * paramVals.psthBinSize), ...
            shuffleHistHist(1:end) ./ histHistSum, ...
            'edgecolor', 'k', 'facecolor', 'k', 'parent', plotAx);
        ylabel(plotAx, 'Frac. Bins')

        % Plot the summary data of number of cells exceeding significance at
        % each time bin
        bar(laserBinTimes, [sum(laserBinMaxSig, 1);  sum(laserBinMinSig, 1)]', ...
            'parent', summaryAx, 'BarWidth', 1);
        ylabel('# Sig Bins');

        % Set the x limits of the histogram of PSTH values to reasonable limits
        histHistMaxXLim = max([lastBin ./ (size(sdata.laserTimes, 1) * paramVals.psthBinSize), ...
            normPsthBarVals(cellNum, paramVals.laserPSTHBins(1)), ...
            normPsthBarVals(cellNum, paramVals.laserPSTHBins(2))]);
    %     maxHistHistVal = (10 * floor(max(histHistMaxXLim) / 10)) + 10;
    %     xlimVals = get(plotAx,'xlim');
        xlim(rasterAx,[paramVals.psthWindow(1), paramVals.psthWindow(2)]);
        xlim(psthAx,[paramVals.psthWindow(1), paramVals.psthWindow(2)]);

        numHBins = size(paramVals.highlightBinTimes, 2);
        hColors = jet(numHBins);
        psthYLimVals = get(psthAx,'ylim');
        histHistYLimVals = get(plotAx,'ylim');
        for hBinCt = 1:numHBins
            hBinNum(hBinCt) = find(paramVals.psthBins >= paramVals.highlightBinTimes(hBinCt), 1, 'first');
            plot([paramVals.psthBins(hBinNum(hBinCt)) + (paramVals.psthBinSize / 2), ...
                paramVals.psthBins(hBinNum(hBinCt)) + (paramVals.psthBinSize / 2)], ...
                psthYLimVals, 'color', hColors(hBinCt,:), 'linewidth', 3, 'parent', psthAx);
            histHistBinNum(hBinCt) = ...
                find(paramVals.histHistBins >= psthBarVals(cellNum, hBinNum(hBinCt)), 1, 'first');
            plot([paramVals.histHistBins(histHistBinNum(hBinCt)), ...
                paramVals.histHistBins(histHistBinNum(hBinCt))] ./ ...
                (size(sdata.laserTimes, 1) * paramVals.psthBinSize), ...
                histHistYLimVals, 'color', hColors(hBinCt,:), ...
                'linewidth', 3, 'parent', plotAx);
        end

        maxHistHistVal = ceil(max(cat(2, histHistMaxXLim, ...
            histHistBinNum ./  (size(sdata.laserTimes, 1) * paramVals.psthBinSize))));
        xlim(plotAx,[0, maxHistHistVal]);

        drawnow;

        % Write data to the output variable
        outVals.laserBinMaxSig = laserBinMaxSig;
        outVals.laserBinMinSig = laserBinMinSig;
        outVals.laserMaxSigSum = sum(laserBinMaxSig, 1);
        outVals.laserMinSigSum = sum(laserBinMinSig, 1);
        outVals.numCells       = numCells;
        outVals.cellNames{cellNum} = cellNames{cellNum};
        outVals.psthBarVals = psthBarVals;
        outVals.firstSigMaxBin = firstSigMaxBin;

        % Pause to wait for user
        if paramVals.pauseDur > 0
            watchoff;
            disp('Click the figure and hit any key to continue to the next cell.');
            pause;                    
        end
    end
    watchoff;

end