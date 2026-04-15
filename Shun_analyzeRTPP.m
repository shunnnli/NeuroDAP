% Shun_analyzeRTPP_updated_sessionType
% 2026/03/24
%
% Folder naming format expected:
%   YYYYMMDD-SessionType-Baseline
%   YYYYMMDD-SessionType-Left
%   YYYYMMDD-SessionType-Right
%
% Example:
%   20260310-Random-Baseline
%   20260310-Random-Left
%   20260318-Reward-Baseline
%   20260318-Reward-Left
%
% SessionType is parsed from the middle token.
% Condition is parsed from the last token.
%
% Group summary is organized by SessionType rather than stimulated side.

%% Load
clear; close all;
addpath(genpath(osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Analysis/NeuroDAP/Methods')));
[~,~,~,~,~,~,bluePurpleRed] = loadColors;

%% Pick animal folders
animalPaths = uipickfiles( ...
    'FilterSpec', osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Project misc/Recordings'), ...
    'Prompt', 'Select one or MORE animal folders (e.g., SL399, SL400, ...)');

if isempty(animalPaths); return; end
if ischar(animalPaths); animalPaths = {animalPaths}; end
nAnimals = numel(animalPaths);

%% Parameters
box = 'large';                         % 'large' or 'small'
defaultStimSide = "right";             % fallback only if a session type has no stimulated session
removeStaticPeriods = false;
noMovementThreshold = 20;
Fs = 20;                               % used only as fallback for duration

%% Arena and binning
switch lower(box)
    case 'large'
        yMidpoint = 350;
        xLimits = [300, 650];
        yLimits = [0, 700];
        nXBins = 30;
        nYBins = 60;
    otherwise
        yMidpoint = 140;
        xLimits = [420, 520];
        yLimits = [15, 280];
        nXBins = 16;
        nYBins = 20;
end

xEdges = linspace(xLimits(1), xLimits(2), nXBins + 1);
yEdges = linspace(yLimits(1), yLimits(2), nYBins + 1);
xCenters = xEdges(1:end-1) + diff(xEdges)/2;
yCenters = yEdges(1:end-1) + diff(yEdges)/2;

%% Colors
leftColor  = [156, 219, 17]./255;
rightColor = [144, 126, 171]./255;
stimColor  = [89 89 89]./255;
rewardColor = [0 158 115]./255;
punishColor = [135 104 247]./255;
ctrlColor  = [179 179 179]./255;

%% Global accumulators
heatmapAccumulator = containers.Map('KeyType', 'char', 'ValueType', 'any');

allSessionSummaryTable = table();
animalTypeSummaryTable = table();

%% ========================== PER-ANIMAL LOOP ==========================
for animalIdx = 1:nAnimals
    animalFolder = animalPaths{animalIdx};
    [~, animalName] = fileparts(animalFolder);

    rawFolders = dir(animalFolder);
    sessionFolders = rawFolders([rawFolders.isdir]);
    sessionFolders = sessionFolders(~ismember({sessionFolders.name}, {'.','..'}));
    nSessions = numel(sessionFolders);

    if nSessions == 0
        warning('No session subfolders found for animal: %s', animalName);
        continue;
    end

    %% Per-session storage
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
        'dataClean', table()), nSessions, 1);

    %% ---------- Load and preprocess each session ----------
    for sessionIdx = 1:nSessions
        sessionName = string(sessionFolders(sessionIdx).name);
        sessionPath = fullfile(sessionFolders(sessionIdx).folder, sessionFolders(sessionIdx).name);

        csvFiles = dir(fullfile(sessionPath, 'times-*.csv'));
        assert(~isempty(csvFiles), 'No times-*.csv found in %s', sessionPath);

        csvPath = fullfile(csvFiles(1).folder, csvFiles(1).name);
        T = readtable(csvPath);

        % Parse session name
        [sessionType, condition, stimSideFromFolder] = parseSessionFolderName(sessionName);

        % Parse time
        T.Item1 = datetime(T.Item1, ...
            'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSSSSSSZ', ...
            'TimeZone', 'UTC');
        t0 = T.Item1(1);
        T.timeMin = minutes(T.Item1 - t0);
        T = T(T.timeMin >= 0, :);

        % Store metadata
        sessions(sessionIdx).name = sessionName;
        sessions(sessionIdx).path = string(sessionPath);
        sessions(sessionIdx).sessionType = sessionType;
        sessions(sessionIdx).condition = condition;
        sessions(sessionIdx).stimSideFromFolder = stimSideFromFolder;
        sessions(sessionIdx).nFramesRaw = height(T);

        % Session duration for weighting
        sessionDurationMin = max(T.timeMin) - min(T.timeMin);
        if sessionDurationMin <= 0
            sessionDurationMin = height(T) / Fs / 60;
        end
        sessions(sessionIdx).durationMin = sessionDurationMin;

        % Remove static periods if requested
        rawX = T.Item2_X;
        rawY = T.Item2_Y;

        if removeStaticPeriods
            staticMask = getStaticPeriod(rawX, rawY, ...
                noMovementThreshold = noMovementThreshold, ...
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
        rightSideMask = (yClean >= yMidpoint);
        leftSideMask  = (yClean <  yMidpoint);

        rightIdx = find(rightSideMask);
        leftIdx  = find(leftSideMask);

        % If right is the stimulated side
        sessions(sessionIdx).timePctIfStimRight = 100 * nnz(rightSideMask) / max(1, numel(yClean));
        sessions(sessionIdx).distStimIfRight    = safeTrajectoryDistance(xClean, yClean, rightIdx);
        sessions(sessionIdx).distCtrlIfRight    = safeTrajectoryDistance(xClean, yClean, leftIdx);

        % If left is the stimulated side
        sessions(sessionIdx).timePctIfStimLeft = 100 * nnz(leftSideMask) / max(1, numel(yClean));
        sessions(sessionIdx).distStimIfLeft    = safeTrajectoryDistance(xClean, yClean, leftIdx);
        sessions(sessionIdx).distCtrlIfLeft    = safeTrajectoryDistance(xClean, yClean, rightIdx);

        % Accumulate heatmap by SessionType + Condition
        heatKey = makeHeatmapKey(sessionType, condition);
        H = histcounts2(yClean, xClean, yEdges, xEdges);

        if isKey(heatmapAccumulator, heatKey)
            heatmapAccumulator(heatKey) = heatmapAccumulator(heatKey) + H;
        else
            heatmapAccumulator(heatKey) = H;
        end
    end

    %% ---------- Finalize metrics within each session type ----------
    sessionTypes = string({sessions.sessionType});
    sessionConditions = string({sessions.condition});
    uniqueTypesThisAnimal = unique(sessionTypes, 'stable');

    for typeIdx = 1:numel(uniqueTypesThisAnimal)
        currentType = uniqueTypesThisAnimal(typeIdx);

        typeMask = (sessionTypes == currentType);
        stimulatedMaskThisType = typeMask & (sessionConditions == "Stimulated");
        baselineMaskThisType   = typeMask & (sessionConditions == "Baseline");

        stimSidesThisType = string({sessions(stimulatedMaskThisType).stimSideFromFolder});
        stimDurationsThisType = [sessions(stimulatedMaskThisType).durationMin];

        leftStimDuration  = sum(stimDurationsThisType(stimSidesThisType == "left"), 'omitnan');
        rightStimDuration = sum(stimDurationsThisType(stimSidesThisType == "right"), 'omitnan');
        totalStimDuration = leftStimDuration + rightStimDuration;

        if totalStimDuration > 0
            baselineWeightLeft  = leftStimDuration  / totalStimDuration;
            baselineWeightRight = rightStimDuration / totalStimDuration;
        else
            baselineWeightLeft  = double(defaultStimSide == "left");
            baselineWeightRight = double(defaultStimSide == "right");
            warning('Animal %s, session type %s has no stimulated session. Using defaultStimSide = %s.', ...
                animalName, currentType, defaultStimSide);
        end

        % Assign metrics to stimulated sessions in this type
        stimIndices = find(stimulatedMaskThisType);
        for k = 1:numel(stimIndices)
            idx = stimIndices(k);

            sessions(idx).baselineWeightLeft = baselineWeightLeft;
            sessions(idx).baselineWeightRight = baselineWeightRight;

            if sessions(idx).stimSideFromFolder == "right"
                sessions(idx).referenceStimSide = "right";
                sessions(idx).timePctStimSide = sessions(idx).timePctIfStimRight;
                sessions(idx).distStimSide = sessions(idx).distStimIfRight;
                sessions(idx).distCtrlSide = sessions(idx).distCtrlIfRight;
            elseif sessions(idx).stimSideFromFolder == "left"
                sessions(idx).referenceStimSide = "left";
                sessions(idx).timePctStimSide = sessions(idx).timePctIfStimLeft;
                sessions(idx).distStimSide = sessions(idx).distStimIfLeft;
                sessions(idx).distCtrlSide = sessions(idx).distCtrlIfLeft;
            else
                error('Stimulated session does not have valid stimSideFromFolder: %s', sessions(idx).name);
            end
        end

        % Assign metrics to baseline sessions in this type
        baselineIndices = find(baselineMaskThisType);
        for k = 1:numel(baselineIndices)
            idx = baselineIndices(k);

            sessions(idx).baselineWeightLeft = baselineWeightLeft;
            sessions(idx).baselineWeightRight = baselineWeightRight;
            sessions(idx).referenceStimSide = sprintf('weighted_left_%0.3f_right_%0.3f', ...
                baselineWeightLeft, baselineWeightRight);

            sessions(idx).timePctStimSide = ...
                baselineWeightRight * sessions(idx).timePctIfStimRight + ...
                baselineWeightLeft  * sessions(idx).timePctIfStimLeft;

            sessions(idx).distStimSide = ...
                baselineWeightRight * sessions(idx).distStimIfRight + ...
                baselineWeightLeft  * sessions(idx).distStimIfLeft;

            sessions(idx).distCtrlSide = ...
                baselineWeightRight * sessions(idx).distCtrlIfRight + ...
                baselineWeightLeft  * sessions(idx).distCtrlIfLeft;
        end

        %% ---------- Per-animal weighted summary for this session type ----------
        baselineDurations = [sessions(baselineMaskThisType).durationMin]';
        stimulatedDurations = [sessions(stimulatedMaskThisType).durationMin]';

        baselineTimePct = [sessions(baselineMaskThisType).timePctStimSide]';
        baselineDistStim = [sessions(baselineMaskThisType).distStimSide]';
        baselineDistCtrl = [sessions(baselineMaskThisType).distCtrlSide]';

        stimulatedTimePct = [sessions(stimulatedMaskThisType).timePctStimSide]';
        stimulatedDistStim = [sessions(stimulatedMaskThisType).distStimSide]';
        stimulatedDistCtrl = [sessions(stimulatedMaskThisType).distCtrlSide]';

        oneRow = table( ...
            string(animalName), ...
            currentType, ...
            weightedMean(baselineTimePct, baselineDurations), ...
            weightedMean(baselineDistStim, baselineDurations), ...
            weightedMean(baselineDistCtrl, baselineDurations), ...
            weightedMean(stimulatedTimePct, stimulatedDurations), ...
            weightedMean(stimulatedDistStim, stimulatedDurations), ...
            weightedMean(stimulatedDistCtrl, stimulatedDurations), ...
            'VariableNames', { ...
                'Animal', ...
                'SessionType', ...
                'BaselineTimePct', ...
                'BaselineDistStim', ...
                'BaselineDistCtrl', ...
                'StimulatedTimePct', ...
                'StimulatedDistStim', ...
                'StimulatedDistCtrl'});

        animalTypeSummaryTable = [animalTypeSummaryTable; oneRow];
    end

    %% ---------- Per-animal figure ----------
    initializeFig(0.9, 1);
    tl = tiledlayout(2, nSessions + 1, 'TileSpacing', 'compact', 'Padding', 'compact');

    for sessionIdx = 1:nSessions
        X = sessions(sessionIdx).dataClean.X;
        Y = sessions(sessionIdx).dataClean.Y;

        nexttile([2 1]); hold on;

        if sessions(sessionIdx).condition == "Baseline"
            plot(X(Y >= yMidpoint), Y(Y >= yMidpoint), 'Color', rightColor, 'LineWidth', 2);
            plot(X(Y <  yMidpoint), Y(Y <  yMidpoint), 'Color', leftColor,  'LineWidth', 2);
            legend({'Right', 'Left'}, 'Location', 'northeast');
        else
            if sessions(sessionIdx).stimSideFromFolder == "right"
                plot(X(Y >= yMidpoint), Y(Y >= yMidpoint), 'Color', stimColor, 'LineWidth', 2);
                plot(X(Y <  yMidpoint), Y(Y <  yMidpoint), 'Color', ctrlColor, 'LineWidth', 2);
            else
                plot(X(Y >= yMidpoint), Y(Y >= yMidpoint), 'Color', ctrlColor, 'LineWidth', 2);
                plot(X(Y <  yMidpoint), Y(Y <  yMidpoint), 'Color', stimColor, 'LineWidth', 2);
            end
            legend({'Control side', 'Stimulated side'}, 'Location', 'northeast');
        end

        xlim(xLimits); ylim(yLimits);
        xlabel('X Position'); ylabel('Y Position');
        title(sprintf('%s | %s', sessions(sessionIdx).sessionType, sessions(sessionIdx).name), ...
            'Interpreter', 'none');
    end

    % Right-top: time on stimulated side (%)
    nexttile; hold on;
    for sessionIdx = 1:nSessions
        barColor = ctrlColor;
        if sessions(sessionIdx).condition == "Stimulated"
            barColor = stimColor;
        end
        plotScatterBar(sessionIdx, sessions(sessionIdx).timePctStimSide, ...
            Color = barColor, style = 'bar');
    end
    xticks(1:nSessions);
    xticklabels(string({sessions.name}));
    xtickangle(25);
    ylabel('Time in stimulated side (%)');

    % Right-bottom: distance on control vs stimulated side
    nexttile; hold on;
    for sessionIdx = 1:nSessions
        plotScatterBar(2*sessionIdx - 0.25, sessions(sessionIdx).distCtrlSide, ...
            Color = ctrlColor, style = 'bar');
        plotScatterBar(2*sessionIdx + 0.25, sessions(sessionIdx).distStimSide, ...
            Color = stimColor, style = 'bar');
    end
    xticks((1:nSessions) * 2);
    xticklabels(string({sessions.name}));
    xtickangle(25);
    ylabel('Distance traveled (pixels)');

    title(tl, sprintf('%s - Summary', animalName), 'Interpreter', 'none');

    % Save per-animal figure
    animalSummaryPNG = fullfile(animalFolder, sprintf('%s_summary.pdf', animalName));
    try
        exportgraphics(gcf, animalSummaryPNG, 'Resolution', 300);
    catch
        saveas(gcf, animalSummaryPNG);
    end
    close(gcf);

    %% ---------- Save per-animal session summary CSV ----------
    animalColumn = repmat(string(animalName), nSessions, 1);

    sessionSummaryTable = table( ...
        animalColumn, ...
        string({sessions.name})', ...
        string({sessions.sessionType})', ...
        string({sessions.condition})', ...
        string({sessions.stimSideFromFolder})', ...
        string({sessions.referenceStimSide})', ...
        [sessions.baselineWeightLeft]', ...
        [sessions.baselineWeightRight]', ...
        [sessions.durationMin]', ...
        [sessions.nFramesRaw]', ...
        [sessions.nFramesAnalyzed]', ...
        [sessions.timePctStimSide]', ...
        [sessions.distStimSide]', ...
        [sessions.distCtrlSide]', ...
        'VariableNames', { ...
            'Animal', ...
            'SessionName', ...
            'SessionType', ...
            'Condition', ...
            'StimSideFromFolder', ...
            'ReferenceStimSideUsedForMetrics', ...
            'BaselineWeightLeft', ...
            'BaselineWeightRight', ...
            'SessionDurationMin', ...
            'RawFrameCount', ...
            'AnalyzedFrameCount', ...
            'TimePctStimulatedSide', ...
            'DistanceStimulatedSide_pixels', ...
            'DistanceControlSide_pixels'});

    animalSummaryCSV = fullfile(animalFolder, sprintf('%s_sessionSummary.csv', animalName));
    writetable(sessionSummaryTable, animalSummaryCSV);

    allSessionSummaryTable = [allSessionSummaryTable; sessionSummaryTable];
end

%% ========================== GROUP SUMMARY FIGURE ==========================
if isempty(allSessionSummaryTable)
    warning('No valid session data found. Group summary not generated.');
    return;
end

sessionTypeList = unique(allSessionSummaryTable.SessionType, 'stable');
nSessionTypes = numel(sessionTypeList);

[groupRoot, ~] = fileparts(animalPaths{1});
groupOutputFolder = fullfile(groupRoot, 'GroupSummary');
if ~exist(groupOutputFolder, 'dir')
    mkdir(groupOutputFolder);
end

initializeFig(1, 1);
tl = tiledlayout(2, nSessionTypes + 1, 'TileSpacing', 'compact', 'Padding', 'compact');

%% ---------- Top row: Baseline heatmaps by session type ----------
for typeIdx = 1:nSessionTypes
    currentType = sessionTypeList(typeIdx);
    H = getHeatmapOrZero(heatmapAccumulator, currentType, "Baseline", nYBins, nXBins);
    H = normalizeHeatmap(H);

    nexttile; hold on;
    imagesc(xCenters, yCenters, H);
    axis xy;
    xlim(xLimits); ylim(yLimits);
    xlabel('X'); ylabel('Y');
    title(sprintf('%s - Baseline', currentType), 'Interpreter', 'none');
    colormap(bluePurpleRed);
    colorbar;
end

%% ---------- Top-right: pooled baseline vs stim/reward/punish ----------
nexttile; hold on;

animalList = unique(string(allSessionSummaryTable.Animal), 'stable');
nAnimalsGroup = numel(animalList);
sessionTypeNames = string(allSessionSummaryTable.SessionType);

hasRandom = any(strcmpi(sessionTypeNames, "Random"));
hasReward = any(strcmpi(sessionTypeNames, "Reward"));
hasPunish = any(strcmpi(sessionTypeNames, "Punish"));

% Per-animal values
pooledBaselinePerAnimal = nan(nAnimalsGroup, 1);
stimPerAnimal   = nan(nAnimalsGroup, 1);   % Random stim
rewardPerAnimal = nan(nAnimalsGroup, 1);
punishPerAnimal = nan(nAnimalsGroup, 1);

for animalIdx = 1:nAnimalsGroup
    currentAnimal = animalList(animalIdx);

    animalMask = string(allSessionSummaryTable.Animal) == currentAnimal;
    baselineMask = animalMask & (string(allSessionSummaryTable.Condition) == "Baseline");

    % pooled baseline across ALL baseline sessions, weighted by duration
    pooledBaselinePerAnimal(animalIdx) = weightedMean( ...
        allSessionSummaryTable.TimePctStimulatedSide(baselineMask), ...
        allSessionSummaryTable.SessionDurationMin(baselineMask));

    % Random -> Stim
    if hasRandom
        randomStimMask = animalMask & ...
                         (string(allSessionSummaryTable.SessionType) == "Random") & ...
                         (string(allSessionSummaryTable.Condition) == "Stimulated");

        stimPerAnimal(animalIdx) = weightedMean( ...
            allSessionSummaryTable.TimePctStimulatedSide(randomStimMask), ...
            allSessionSummaryTable.SessionDurationMin(randomStimMask));
    end

    % Reward
    if hasReward
        rewardStimMask = animalMask & ...
                         (string(allSessionSummaryTable.SessionType) == "Reward") & ...
                         (string(allSessionSummaryTable.Condition) == "Stimulated");

        rewardPerAnimal(animalIdx) = weightedMean( ...
            allSessionSummaryTable.TimePctStimulatedSide(rewardStimMask), ...
            allSessionSummaryTable.SessionDurationMin(rewardStimMask));
    end

    % Punish
    if hasPunish
        punishStimMask = animalMask & ...
                         (string(allSessionSummaryTable.SessionType) == "Punish") & ...
                         (string(allSessionSummaryTable.Condition) == "Stimulated");

        punishPerAnimal(animalIdx) = weightedMean( ...
            allSessionSummaryTable.TimePctStimulatedSide(punishStimMask), ...
            allSessionSummaryTable.SessionDurationMin(punishStimMask));
    end
end

% Plot order
xBaseline = 1;
xStim     = 2;
xReward   = 3;
xPunish   = 4;

% Bars
plotScatterBar(xBaseline, pooledBaselinePerAnimal(~isnan(pooledBaselinePerAnimal)), ...
    Color = ctrlColor, style = 'bar');


if hasRandom
    plotScatterBar(xStim, stimPerAnimal(~isnan(stimPerAnimal)), ...
        Color = stimColor, style = 'bar');
end

if hasReward
    plotScatterBar(xReward, rewardPerAnimal(~isnan(rewardPerAnimal)), ...
        Color = rewardColor, style = 'bar');
end

if hasPunish
    plotScatterBar(xPunish, punishPerAnimal(~isnan(punishPerAnimal)), ...
        Color = punishColor, style = 'bar');
end

xticks([xBaseline, xStim, xReward, xPunish]);
xticklabels({'Baseline', 'Stim', 'Reward', 'Punish'});
xtickangle(20);
ylabel('Time in stimulated side (%)');
title('Across animals');

% Baseline vs Baseline stimulated
if hasRandom
    baselineForStim = pooledBaselinePerAnimal(~isnan(pooledBaselinePerAnimal));
    stimVals = stimPerAnimal(~isnan(stimPerAnimal));

    if ~isempty(baselineForStim) && ~isempty(stimVals)
        plotStats(baselineForStim, stimVals, [xBaseline, xStim], testType = 'kstest');
    end
end

% Baseline vs Reward
if hasRandom && hasReward 
    rewardVals = rewardPerAnimal(~isnan(rewardPerAnimal));
    baselineForStim = pooledBaselinePerAnimal(~isnan(pooledBaselinePerAnimal));

    if ~isempty(baselineForStim) && ~isempty(rewardVals)
        plotStats(baselineForStim, rewardVals, [xBaseline, xReward], testType = 'kstest');
    end
end

% Baseline vs Punish
if hasRandom && hasPunish 
    punishVals = punishPerAnimal(~isnan(punishPerAnimal));
    baselineForStim = pooledBaselinePerAnimal(~isnan(pooledBaselinePerAnimal));

    if ~isempty(baselineForStim) && ~isempty(punishVals)
        plotStats(baselineForStim, punishVals, [xBaseline, xPunish], testType = 'kstest');
    end
end

% Reward vs Punish
if hasReward && hasPunish
    rewardVals = rewardPerAnimal(~isnan(rewardPerAnimal));
    punishVals = punishPerAnimal(~isnan(punishPerAnimal));

    if ~isempty(rewardVals) && ~isempty(punishVals)
        plotStats(rewardVals, punishVals, [xReward, xPunish], testType = 'kstest');
    end
end

%% ---------- Bottom row: Stimulated heatmaps by session type ----------
for typeIdx = 1:nSessionTypes
    currentType = sessionTypeList(typeIdx);
    H = getHeatmapOrZero(heatmapAccumulator, currentType, "Stimulated", nYBins, nXBins);
    H = normalizeHeatmap(H);

    nexttile; hold on;
    imagesc(xCenters, yCenters, H);
    axis xy;
    xlim(xLimits); ylim(yLimits);
    xlabel('X'); ylabel('Y');
    title(sprintf('%s - Stimulated', currentType), 'Interpreter', 'none');
    colormap(bluePurpleRed);
    colorbar;
end

%% ---------- Bottom-right: Distance plot by session type ----------
% For each session type:
%   baseline pair  = control vs stimulated side
%   stimulated pair = control vs stimulated side
nexttile; hold on;
distGapOuter = 0.22;
distGapInner = 0.08;

for typeIdx = 1:nSessionTypes
    currentType = sessionTypeList(typeIdx);
    rowMask = (animalTypeSummaryTable.SessionType == currentType);

    baselineCtrl = animalTypeSummaryTable.BaselineDistCtrl(rowMask);
    baselineCtrl = baselineCtrl(~isnan(baselineCtrl));

    baselineStim = animalTypeSummaryTable.BaselineDistStim(rowMask);
    baselineStim = baselineStim(~isnan(baselineStim));

    stimulatedCtrl = animalTypeSummaryTable.StimulatedDistCtrl(rowMask);
    stimulatedCtrl = stimulatedCtrl(~isnan(stimulatedCtrl));

    stimulatedStim = animalTypeSummaryTable.StimulatedDistStim(rowMask);
    stimulatedStim = stimulatedStim(~isnan(stimulatedStim));

    xBaseCtrl = typeIdx - distGapOuter - distGapInner;
    xBaseStim = typeIdx - distGapOuter + distGapInner;
    xStimCtrl = typeIdx + distGapOuter - distGapInner;
    xStimStim = typeIdx + distGapOuter + distGapInner;

    if ~isempty(baselineCtrl)
        plotScatterBar(xBaseCtrl, baselineCtrl, Color = ctrlColor, style = 'bar');
    end
    if ~isempty(baselineStim)
        plotScatterBar(xBaseStim, baselineStim, Color = stimColor, style = 'bar');
    end
    if ~isempty(stimulatedCtrl)
        plotScatterBar(xStimCtrl, stimulatedCtrl, Color = ctrlColor, style = 'bar');
    end
    if ~isempty(stimulatedStim)
        plotScatterBar(xStimStim, stimulatedStim, Color = stimColor, style = 'bar');
    end

    if ~isempty(baselineCtrl) && ~isempty(baselineStim)
        plotStats(baselineCtrl, baselineStim, [xBaseCtrl, xBaseStim], testType = 'kstest');
    end
    if ~isempty(stimulatedCtrl) && ~isempty(stimulatedStim)
        plotStats(stimulatedCtrl, stimulatedStim, [xStimCtrl, xStimStim], testType = 'kstest');
    end
end

xticks(1:nSessionTypes);
xticklabels(sessionTypeList);
xtickangle(20);
ylabel('Distance traveled (pixels)');
title('Distance by session type');
legend({'Control side', 'Stimulated side'}, 'Location', 'best');

title(tl, sprintf('RTPP Group Summary (n = %d animals)', nAnimals));

groupSummaryPNG = fullfile(groupOutputFolder, 'AcrossAnimals_summary_bySessionType.pdf');
try
    exportgraphics(gcf, groupSummaryPNG, 'Resolution', 300);
catch
    saveas(gcf, groupSummaryPNG);
end

disp(['Saved group summary: ' groupSummaryPNG]);

%% Optional: save the per-animal weighted session-type summary table
groupSummaryCSV = fullfile(groupOutputFolder, 'AcrossAnimals_summary_bySessionType.csv');
writetable(animalTypeSummaryTable, groupSummaryCSV);
disp(['Saved group summary CSV: ' groupSummaryCSV]);

%% Optional: pick corresponding recording sessions before RTPP

groupOutputFolder = '/Volumes/Neurobio/MICROSCOPE/Shun/Project valence/Recordings/202511-RTPP-EPStim/GroupSummary';

sessionList_reward = uipickfiles( ...
    'FilterSpec', osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Project valence/Recordings/202511-RTPP-EPStim/'), ...
    'Prompt', 'Select reward session folders');

sessionList_punish = uipickfiles( ...
    'FilterSpec', osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Project valence/Recordings/202511-RTPP-EPStim/'), ...
    'Prompt', 'Select punish session folders');

sessionList_punish = sessionList_punish([2 3 4 1 5 6 7 8]);
sessionList_reward = sessionList_reward';
sessionList_punish = sessionList_punish';

save(fullfile(groupOutputFolder,'sessionList.mat'),'sessionList_reward','sessionList_punish');

%% Correlation between behavioral and location preference
% Step 1: find corresponding recording sessions for each RTPP sessions
% Step 2: read behavior_sessionName.mat and load trialTable
% Step 3: calculate average of relevant behavior variables per trial type
% Step 4: find a way to store these things across animals

behVariables = {'ReactionTime', 'nLicks', 'nAnticipatoryLicks', 'nBaselineLicks', 'nOutcomeLicks', 'stageAmp'};
leftColor  = [156, 219, 17]./255;
rightColor = [144, 126, 171]./255;
stimColor  = [89 89 89]./255;
rewardColor = [0 158 115]./255;
punishColor = [135 104 247]./255;
ctrlColor  = [179 179 179]./255;

% Read RTPP summary csv
groupOutputFolder = '/Volumes/Neurobio/MICROSCOPE/Shun/Project valence/Recordings/202511-RTPP-EPStim/GroupSummary';
csvpath = fullfile(groupOutputFolder, 'AcrossAnimals_summary_bySessionType.csv');
rtppSummary = readtable(csvpath);

% Get reward session table
load(fullfile(groupOutputFolder,'sessionList.mat'));
rewardSummaryTable = getBehaviorSummary(sessionList_reward, statsType='avg', variables=behVariables);
punishSummaryTable = getBehaviorSummary(sessionList_punish, statsType='avg', variables=behVariables);

save(fullfile(groupOutputFolder,'sessionSummaryTable.mat'),'rewardSummaryTable','punishSummaryTable');

%% Final: correlate averages with position of mouse in box
for i = 1:length(behVariables)
    behVar = behVariables{i};
    behCol = "StimOnly_" + behVar; % change to Pair_ if that is what you want

    initializeFig(0.5,0.5);
    tiledlayout(1,2);

    % Panel 1: behavior variable vs stimulated percentage
    nexttile; hold on;

    x_reward = rewardSummaryTable.(behCol);
    y_reward = rtppSummary.StimulatedTimePct(strcmpi(rtppSummary.SessionType, 'Reward'));

    x_punish = punishSummaryTable.(behCol);
    y_punish = rtppSummary.StimulatedTimePct(strcmpi(rtppSummary.SessionType, 'Punish'));

    scatter(x_reward, y_reward, 150, rewardColor, 'filled');
    scatter(x_punish, y_punish, 150, punishColor, 'filled');
    plotFit(x_reward, y_reward, color=rewardColor, LineWidth=5);
    plotFit(x_punish, y_punish, color=punishColor, LineWidth=5);

    % Correlation stats
    [r_reward, p_reward] = corr(x_reward, y_reward, 'rows', 'complete');
    [r_punish, p_punish] = corr(x_punish, y_punish, 'rows', 'complete');

    xlabel(strrep(behCol, '_', '\_'));
    ylabel('Stimulated side (%)');
    title(sprintf('%s vs Stim%%', behVar), 'Interpreter', 'none');
    legend({sprintf('reward (r=%.2f, p=%.3f)', r_reward, p_reward), ...
            sprintf('punish (r=%.2f, p=%.3f)', r_punish, p_punish)}, ...
            'Location', 'best');


    % Panel 2: behavior variable vs distance
    nexttile; hold on;

    x_reward = rewardSummaryTable.(behCol);
    y_reward = rtppSummary.StimulatedDistStim(strcmpi(rtppSummary.SessionType, 'Reward'));

    x_punish = punishSummaryTable.(behCol);
    y_punish = rtppSummary.StimulatedDistStim(strcmpi(rtppSummary.SessionType, 'Punish'));

    scatter(x_reward, y_reward, 150, rewardColor, 'filled');
    scatter(x_punish, y_punish, 150, punishColor, 'filled');
    plotFit(x_reward, y_reward, color=rewardColor, LineWidth=5);
    plotFit(x_punish, y_punish, color=punishColor, LineWidth=5);

    % Correlation stats
    [r_reward, p_reward] = corr(x_reward, y_reward, 'rows', 'complete');
    [r_punish, p_punish] = corr(x_punish, y_punish, 'rows', 'complete');

    xlabel(strrep(behCol, '_', '\_'));
    ylabel('Distance');
    title(sprintf('%s vs Distance', behVar), 'Interpreter', 'none');
    legend({sprintf('reward (r=%.2f, p=%.3f)', r_reward, p_reward), ...
            sprintf('punish (r=%.2f, p=%.3f)', r_punish, p_punish)}, ...
            'Location', 'best');

    saveFigures(gcf,['Correlation-',behVar],groupOutputFolder);
end





%% ========================== LOCAL FUNCTIONS ==========================

function [sessionType, condition, stimSideFromFolder] = parseSessionFolderName(sessionName)
    % Expected format:
    %   YYYYMMDD-SessionType-Baseline
    %   YYYYMMDD-SessionType-Left
    %   YYYYMMDD-SessionType-Right

    parts = split(string(sessionName), "-");

    if numel(parts) < 3
        error('Session folder name "%s" does not match expected format YYYYMMDD-SessionType-Condition', sessionName);
    end

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

function key = makeHeatmapKey(sessionType, condition)
    key = sprintf('%s__%s', char(sessionType), char(condition));
end

function H = getHeatmapOrZero(heatmapAccumulator, sessionType, condition, nYBins, nXBins)
    key = makeHeatmapKey(sessionType, condition);
    if isKey(heatmapAccumulator, key)
        H = heatmapAccumulator(key);
    else
        H = zeros(nYBins, nXBins);
    end
end

function H = normalizeHeatmap(H)
    totalCount = sum(H(:));
    if totalCount <= 0
        return;
    end
    H = H ./ totalCount;
end

function value = weightedMean(values, weights)
    values = values(:);
    weights = weights(:);

    validMask = ~isnan(values) & ~isnan(weights) & (weights > 0);

    if ~any(validMask)
        value = NaN;
        return;
    end

    value = sum(values(validMask) .* weights(validMask)) / sum(weights(validMask));
end

function dist = safeTrajectoryDistance(X, Y, idx)
    if isempty(idx)
        dist = 0;
    else
        dist = getTrajectoryDistance(X, Y, filter = idx);
    end
end



%% Method fo extracting behavior sesssion summary

function results = getBehaviorSummary(sessionList, options)

arguments
    sessionList
    options.variables = {'ReactionTime', 'nLicks', 'nAnticipatoryLicks', 'nBaselineLicks', 'nOutcomeLicks', 'stageAmp'}
    options.statsType string = 'avg' % can also be median
end

sessionList = string(sessionList);
vars = options.variables;
nVars = numel(vars);
nSessions = numel(sessionList);

animalIDs = strings(nSessions, 1);
tasks = strings(nSessions, 1);
allStimOnly = nan(nSessions, nVars);
allPair = nan(nSessions, nVars);

for s = 1:nSessions
    sessionpath = sessionList(s);
    pathsplit = strsplit(sessionpath, '/');
    sessionName = pathsplit{end};
    sessionsplit = strsplit(sessionName,'-');
    animalIDs(s) = string(sessionsplit{2});
    
    % Find task of the session
    if contains(sessionpath, 'reward', IgnoreCase=true)
        tasks(s) = "reward";
    elseif contains(sessionpath, 'punish', IgnoreCase=true)
        tasks(s) = "punish";
    else
        tasks(s) = "random";
    end
    
    % Load session files
    load(fullfile(sessionpath,['analysis_',sessionName,'.mat']), 'analysis');
    summary_stimOnly = nan(1, nVars);
    summary_pair = nan(1, nVars);
    
    % Get stim only or pair trials row index
    stimOnlyRow_dLight = find(strcmpi({analysis.event}, 'Stim only') & strcmpi({analysis.name}, 'dLight'), 1);
    stimOnlyRow_lick = find(strcmpi({analysis.event}, 'Stim only') & strcmpi({analysis.name}, 'Lick'), 1);
    pairRow_dLight = find(strcmpi({analysis.event}, 'Pair') & strcmpi({analysis.name}, 'dLight'), 1);
    pairRow_lick = find(strcmpi({analysis.event}, 'Pair') & strcmpi({analysis.name}, 'Lick'), 1);

    % Get stim only trial lick data
    trialTable_stimOnly = analysis(stimOnlyRow_lick).trialInfo.trialTable;
    trialTable_pair = analysis(pairRow_lick).trialInfo.trialTable;

    % Get stageAmp for stim only trial dLight data
    stageMax = analysis(stimOnlyRow_dLight).stageMax.data;
    stageMin = analysis(stimOnlyRow_dLight).stageMin.data;
    stageAmp_stimOnly = getAmplitude(stageMax, stageMin);

    % Get stageAmp for pair trial dLight data
    stageMax = analysis(pairRow_dLight).stageMax.data;
    stageMin = analysis(pairRow_dLight).stageMin.data;
    stageAmp_pair = getAmplitude(stageMax, stageMin);

    % Calculate mean/median for each variable
    for v = 1:nVars
        varName = vars{v};
        if ~contains(varName,'stage')
            behCol_stimOnly = trialTable_stimOnly.(varName);
            behCol_pair = trialTable_pair.(varName);
        
            if strcmpi(options.statsType, 'mean') || strcmpi(options.statsType, 'avg')
                summary_stimOnly(v) = mean(behCol_stimOnly);
                summary_pair(v) = mean(behCol_pair);
            elseif strcmpi(options.statsType, 'median')
                summary_stimOnly(v) = median(behCol_stimOnly);
                summary_pair(v) = median(behCol_pair);
            else
                error('options.statsType must be "avg", "mean", or "median".')
            end
        else
            if strcmpi(options.statsType, 'mean') || strcmpi(options.statsType, 'avg')
                summary_stimOnly(v) = mean(stageAmp_stimOnly(:,2));
                summary_pair(v) = mean(stageAmp_pair(:,2));
            elseif strcmpi(options.statsType, 'median')
                summary_stimOnly(v) = median(stageAmp_stimOnly(:,2));
                summary_pair(v) = median(stageAmp_pair(:,2));
            else
                error('options.statsType must be "avg", "mean", or "median".')
            end
    
        end
    end

    % Store session data
    allStimOnly(s, :) = summary_stimOnly;
    allPair(s, :) = summary_pair;

    % Remove loaded analysis session
    clearvars analysis
    disp(['Finished: session ', num2str(s), ': ',sessionName]);
end

% Build results table dynamically
results = table(animalIDs, tasks, 'VariableNames', {'Animal', 'Task'});

stimNames = matlab.lang.makeValidName("StimOnly_" + string(vars));
pairNames = matlab.lang.makeValidName("Pair_" + string(vars));

for v = 1:nVars
    results.(stimNames(v)) = allStimOnly(:, v);
    results.(pairNames(v)) = allPair(:, v);
end

end