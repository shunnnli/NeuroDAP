function [sessions, heatmapAccumulator] = analyzeSession_RTPP(sessionpath, options)

    arguments
        sessionpath string
        
        options.box string = "clamp"

        options.heatmapAccumulator
        options.sessions
        options.sessionIdx
    
        options.removeStaticPeriods logical = false
        options.noMovementThreshold double = 20
        options.Fs double = 30
    
        options.timeRange double = [-5, 10]
        options.savePDF logical = true
        options.savePNG logical = false
    end

    %% Initialize variables if needed
    if ~isfield(options,'heatmapAccumulator')
        heatmapAccumulator = containers.Map('KeyType', 'char', 'ValueType', 'any');
    else
        heatmapAccumulator = options.heatmapAccumulator;
    end

    if ~isfield(options,'sessions')
        sessions = repmat(struct( ...
        'name', "", ...
        'path', "", ...
        'sessionType', "", ...
        'condition', "", ...                 % "Baseline" or "Stimulated"
        'stimSideFromFolder', "", ...        % "left", "right", or ""
        'referenceStimSide', "", ...         % what side was used to define stim-side metrics
        'baselineWeightLeft', NaN, ...
        'baselineWeightRight', NaN, ...
        'durationMin', NaN, ...
        'nFramesRaw', NaN, ...
        'nFramesAnalyzed', NaN, ...
        'timePctIfStimRight', NaN, ...
        'timePctIfStimLeft', NaN, ...
        'distStimIfRight', NaN, ...
        'distCtrlIfRight', NaN, ...
        'distStimIfLeft', NaN, ...
        'distCtrlIfLeft', NaN, ...
        'timePctStimSide', NaN, ...
        'distStimSide', NaN, ...
        'distCtrlSide', NaN, ...
        'dataClean', table()), 1, 1);
    else
        sessions = options.sessions;
    end

    if ~isfield(options,'sessionIdx')
        sessionIdx = height(sessions) + 1;
    else
        sessionIdx = options.sessionIdx;
    end

    %% Arena and binning
    switch lower(options.box)
        case 'large'
            yMidpoint = 350;
            xLimits = [300, 650];
            yLimits = [0, 700];
            nXBins = 30;
            nYBins = 60;
        case 'clamp'
            yMidpoint = 234;
            xLimits = [50, 250];
            yLimits = [0, 470];
            nXBins = 20;
            nYBins = 40;
        otherwise
            yMidpoint = 140;
            xLimits = [420, 520];
            yLimits = [15, 280];
            nXBins = 16;
            nYBins = 20;
    end
    
    xEdges = linspace(xLimits(1), xLimits(2), nXBins + 1);
    yEdges = linspace(yLimits(1), yLimits(2), nYBins + 1);

    %% Load paths
    sessionpath = regexprep(sessionpath, [regexptranslate('escape', filesep) '+$'], '');
    [~, sessionName, ext] = fileparts(sessionpath);
    sessionName = sessionName + ext;

    if contains(options.box,'clamp',IgnoreCase=true)
        csvFiles = dir(fullfile(sessionpath, 'camera-processed', 'times-processed.csv'));
    else
        csvFiles = dir(fullfile(sessionpath, 'times-*.csv'));
    end
    assert(~isempty(csvFiles), 'No times-*.csv found in %s', sessionpath);
    
    csvPath = fullfile(csvFiles(1).folder, csvFiles(1).name);
    T = readtable(csvPath);
    
    % Parse session name
    [sessionType, condition, stimSideFromFolder] = parseSessionFolderName(sessionName,options.box);
    
    % Parse camera data
    if contains(options.box,'clamp',IgnoreCase=true)
        rawX = T.X;
        rawY = T.Y;

        T.Date = datetime(string(T.Date), ...
            'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSSSSSSZ', ...
            'TimeZone', 'UTC');
        t0 = T.Date(1);
        T.timeMin = minutes(T.Date - t0);
        T = T(T.timeMin >= 0, :);
    else
        T.Item1 = datetime(T.Item1, ...
            'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSSSSSSZ', ...
            'TimeZone', 'UTC');
        t0 = T.Item1(1);
        T.timeMin = minutes(T.Item1 - t0);
        T = T(T.timeMin >= 0, :);

        rawX = T.Item2_X;
        rawY = T.Item2_Y;
    end
    
    % Store metadata
    sessions(sessionIdx).name = sessionName;
    sessions(sessionIdx).path = string(sessionpath);
    sessions(sessionIdx).sessionType = sessionType;
    sessions(sessionIdx).condition = condition;
    sessions(sessionIdx).stimSideFromFolder = stimSideFromFolder;
    sessions(sessionIdx).nFramesRaw = height(T);
    
    % Session duration for weighting
    sessionDurationMin = max(T.timeMin) - min(T.timeMin);
    if sessionDurationMin <= 0
        sessionDurationMin = height(T) / options.Fs / 60;
    end
    sessions(sessionIdx).durationMin = sessionDurationMin;
    
    % Remove static periods if requested
    
    
    if options.removeStaticPeriods
        staticMask = getStaticPeriod(rawX, rawY, ...
            noMovementThreshold = options.noMovementThreshold, ...
            windowDuration = 30);
        keepMask = ~staticMask;
    else
        keepMask = true(size(rawX));
    end
    
    xClean = rawX(keepMask);
    yClean = rawY(keepMask);
    tClean = T.timeMin(keepMask);
    
    sessions(sessionIdx).nFramesAnalyzed = numel(xClean);
    sessions(sessionIdx).dataClean = table(tClean, xClean, yClean, ...
        'VariableNames', {'timeMin', 'X', 'Y'});
    
    % Compute metrics under both possible definitions of "stimulated side"
    if contains(options.box,'clamp',IgnoreCase=true)
        rightSideMask = (yClean <  yMidpoint);
        leftSideMask  = (yClean >= yMidpoint);
    else
        rightSideMask = (yClean >= yMidpoint);
        leftSideMask  = (yClean <  yMidpoint);
    end
    
    rightIdx = find(rightSideMask);
    leftIdx  = find(leftSideMask);
    
    % If right is the stimulated side
    sessions(sessionIdx).timePctIfStimRight = 100 * nnz(rightSideMask) / max(1, numel(yClean));
    sessions(sessionIdx).distStimIfRight    = getDistancePerMin(xClean, yClean, rightIdx, options.Fs);
    sessions(sessionIdx).distCtrlIfRight    = getDistancePerMin(xClean, yClean, leftIdx, options.Fs);
    
    % If left is the stimulated side
    sessions(sessionIdx).timePctIfStimLeft = 100 * nnz(leftSideMask) / max(1, numel(yClean));
    sessions(sessionIdx).distStimIfLeft    = getDistancePerMin(xClean, yClean, leftIdx, options.Fs);
    sessions(sessionIdx).distCtrlIfLeft    = getDistancePerMin(xClean, yClean, rightIdx, options.Fs);
    
    % Accumulate heatmap by SessionType + Condition
    heatKey = makeHeatmapKey(sessionType, condition);
    H = histcounts2(yClean, xClean, yEdges, xEdges);
    
    if isKey(heatmapAccumulator, heatKey)
        heatmapAccumulator(heatKey) = heatmapAccumulator(heatKey) + H;
    else
        heatmapAccumulator(heatKey) = H;
    end

    %% Extract DA dynamics
end


%% Local functions
function [sessionType, condition, stimSideFromFolder] = parseSessionFolderName(sessionName,box)
    % Expected format:
    %   YYYYMMDD-SessionType-Baseline
    %   YYYYMMDD-SessionType-Left
    %   YYYYMMDD-SessionType-Right

    parts = split(string(sessionName), "-");

    if numel(parts) < 3
        error('Session folder name "%s" does not match expected format YYYYMMDD-SessionType-Condition', sessionName);
    end

    if contains(box,'clamp',IgnoreCase=true)
        sessionType = string(parts(3));
        lastToken = lower(string(parts(end)));
    
        if contains(lastToken,'random',IgnoreCase=true) || contains(lastToken,'baseline',IgnoreCase=true)
            condition = "Baseline";
            stimSideFromFolder = "";
        elseif contains(lastToken,'Left',IgnoreCase=true)
            condition = "Stimulated";
            stimSideFromFolder = "left";
        elseif contains(lastToken,'Right',IgnoreCase=true)
            condition = "Stimulated";
            stimSideFromFolder = "right";
        else
            error('Unrecognized last token "%s" in session folder "%s". Expected baseline/random, Left, or Right.', ...
                lastToken, sessionName);
        end

    else
        sessionType = string(parts(2));
        lastToken = lower(string(parts(end)));
    
        switch lastToken
            case "baseline"
                condition = "Baseline";
                stimSideFromFolder = "";
            case "left"
                condition = "Stimulated";
                stimSideFromFolder = "left";
            case "right"
                condition = "Stimulated";
                stimSideFromFolder = "right";
            otherwise
                error('Unrecognized last token "%s" in session folder "%s". Expected Baseline, Left, or Right.', ...
                    lastToken, sessionName);
        end
    end
end


function key = makeHeatmapKey(sessionType, condition)
    key = sprintf('%s__%s', char(sessionType), char(condition));
end


function dist = getDistancePerMin(X, Y, idx, Fs)
    if isempty(idx)
        total_dist = 0;
    else
        total_dist = getTrajectoryDistance(X, Y, filter = idx);
    end
    dist = total_dist / (length(X)/Fs);
end
