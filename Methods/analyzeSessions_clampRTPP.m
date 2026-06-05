function [allSessionSummaryTable, animalTypeSummaryTable] = analyzeSessions_clampRTPP(animalPaths, options)

arguments
    animalPaths = string.empty

    options.outputName string = "clampRTPP"
    options.outputFolder string = ""

    options.box string = "large"                  % 'large' or 'small'
    options.removeStaticPeriods logical = false
    options.noMovementThreshold double = 20
    options.Fs double = 20                        % fallback frame rate for times*.csv

    options.daSignalName string = "dLight"
    options.daTimeRange double = [-5, 10]
    options.plotDA logical = true

    options.plotAnimal logical = true
    options.plotGroup logical = true
    options.savePDF logical = true
    options.savePNG logical = false
    options.closeFigures logical = true
end

%% Notes
% analyzeSessions_clampRTPP
% Shun Li / Codex, 2026/06/02
%
% RTPP-style analysis for clamp sessions.
% Unlike Shun_analyzeRTPP.m, the stimulated/control label is read from
% times*.csv ClampON instead of the session folder side token.

%% Pick animal folders if not provided

if isempty(animalPaths)
    animalPaths = uipickfiles( ...
        'FilterSpec', osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Project misc/Recordings'), ...
        'Prompt', 'Select one or MORE animal folders');
end

if isempty(animalPaths)
    allSessionSummaryTable = table();
    animalTypeSummaryTable = table();
    return
end

animalPaths = string(animalPaths);
animalPaths = animalPaths(:);
nAnimals = numel(animalPaths);

%% Arena and binning

switch lower(options.box)
    case "large"
        xLimits = [300, 650];
        yLimits = [0, 700];
        nXBins = 30;
        nYBins = 60;
    otherwise
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

[~,~,~,~,~,~,bluePurpleRed] = loadColors;
clampColor = [.232 .76 .58];
nonClampColor = [217, 237, 223]./255;
ctrlColor = [179 179 179]./255;

%% Global accumulators

heatmapAccumulator = containers.Map('KeyType', 'char', 'ValueType', 'any');
allSessionSummaryTable = table();
animalTypeSummaryTable = table();

%% ========================== PER-ANIMAL LOOP ==========================

for animalIdx = 1:nAnimals
    inputPath = char(animalPaths(animalIdx));

    if hasTimesCsv(inputPath)
        sessionPathInput = inputPath;
        [animalFolder, sessionFolderName] = fileparts(sessionPathInput);
        [~, animalName] = fileparts(animalFolder);
        sessionFolders = struct( ...
            'name', sessionFolderName, ...
            'folder', animalFolder, ...
            'date', '', ...
            'bytes', 0, ...
            'isdir', true, ...
            'datenum', NaN);
    else
        animalFolder = inputPath;
        [~, animalName] = fileparts(animalFolder);

        rawFolders = dir(animalFolder);
        sessionFolders = rawFolders([rawFolders.isdir]);
        sessionFolders = sessionFolders(~ismember({sessionFolders.name}, {'.','..','GroupSummary'}));
    end

    nRawSessions = numel(sessionFolders);

    if nRawSessions == 0
        warning('No session subfolders found for animal: %s', animalName);
        continue
    end

    sessions = repmat(emptySessionStruct(), 0, 1);

    %% ---------- Load and preprocess each session ----------
    for sessionIdx = 1:nRawSessions
        sessionName = string(sessionFolders(sessionIdx).name);
        sessionPath = fullfile(sessionFolders(sessionIdx).folder, sessionFolders(sessionIdx).name);

        try
            [T, csvPath] = readClampTimesCsv(sessionPath, options.Fs);
        catch ME
            warning('Skipping %s: %s', sessionPath, ME.message);
            continue
        end

        [sessionType, folderCondition] = parseClampSessionFolderName(sessionName);

        sessionDurationMin = max(T.timeMin) - min(T.timeMin);
        if sessionDurationMin <= 0 || isnan(sessionDurationMin)
            sessionDurationMin = height(T) / options.Fs / 60;
        end

        rawX = T.X;
        rawY = T.Y;
        rawClampON = T.ClampON;

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
        clampClean = rawClampON(keepMask);
        nonClampClean = ~clampClean;

        timePctClamp = 100 * nnz(clampClean) / max(1, numel(clampClean));
        timePctNonClamp = 100 * nnz(nonClampClean) / max(1, numel(nonClampClean));
        distClamp = safeTrajectoryDistance(xClean, yClean, clampClean);
        distNonClamp = safeTrajectoryDistance(xClean, yClean, nonClampClean);

        clampDurationMin = sessionDurationMin * nnz(clampClean) / max(1, numel(clampClean));
        nonClampDurationMin = sessionDurationMin * nnz(nonClampClean) / max(1, numel(nonClampClean));

        HNonClamp = histcounts2(yClean(nonClampClean), xClean(nonClampClean), yEdges, xEdges);
        HClamp = histcounts2(yClean(clampClean), xClean(clampClean), yEdges, xEdges);
        heatmapAccumulator = addHeatmap(heatmapAccumulator, sessionType, "NonClamp", HNonClamp);
        heatmapAccumulator = addHeatmap(heatmapAccumulator, sessionType, "Clamp", HClamp);

        daStats = analyzeDopamineForClampSession( ...
            sessionPath, sessionName, T, rawClampON, ...
            options, clampColor, nonClampColor);

        session = emptySessionStruct();
        session.name = sessionName;
        session.path = string(sessionPath);
        session.csvPath = string(csvPath);
        session.sessionType = sessionType;
        session.folderCondition = folderCondition;
        session.durationMin = sessionDurationMin;
        session.clampDurationMin = clampDurationMin;
        session.nonClampDurationMin = nonClampDurationMin;
        session.nFramesRaw = height(T);
        session.nFramesAnalyzed = numel(xClean);
        session.nClampFrames = nnz(clampClean);
        session.nNonClampFrames = nnz(nonClampClean);
        session.timePctClamp = timePctClamp;
        session.timePctNonClamp = timePctNonClamp;
        session.distClamp = distClamp;
        session.distNonClamp = distNonClamp;
        session.dataClean = table(tClean, xClean, yClean, clampClean, ...
            'VariableNames', {'timeMin', 'X', 'Y', 'ClampON'});
        session.da = daStats;

        sessions(end+1, 1) = session; %#ok<AGROW>
    end

    if isempty(sessions)
        warning('No valid clamp RTPP sessions found for animal: %s', animalName);
        continue
    end

    %% ---------- Per-animal weighted summary by session type ----------
    sessionTypes = string({sessions.sessionType});
    uniqueTypesThisAnimal = unique(sessionTypes, 'stable');

    for typeIdx = 1:numel(uniqueTypesThisAnimal)
        currentType = uniqueTypesThisAnimal(typeIdx);
        typeMask = sessionTypes == currentType;
        typeSessions = sessions(typeMask);
        durations = [typeSessions.durationMin]';

        daClampSamples = arrayfun(@(s) s.da.nClampSamples, typeSessions)';
        daNonClampSamples = arrayfun(@(s) s.da.nNonClampSamples, typeSessions)';

        oneRow = table( ...
            string(animalName), ...
            currentType, ...
            weightedMean([typeSessions.timePctClamp]', durations), ...
            weightedMean([typeSessions.timePctNonClamp]', durations), ...
            weightedMean([typeSessions.distClamp]', durations), ...
            weightedMean([typeSessions.distNonClamp]', durations), ...
            weightedMean(arrayfun(@(s) s.da.meanClamp, typeSessions)', daClampSamples), ...
            weightedMean(arrayfun(@(s) s.da.varClamp, typeSessions)', daClampSamples), ...
            weightedMean(arrayfun(@(s) s.da.meanNonClamp, typeSessions)', daNonClampSamples), ...
            weightedMean(arrayfun(@(s) s.da.varNonClamp, typeSessions)', daNonClampSamples), ...
            sum(daClampSamples, 'omitnan'), ...
            sum(daNonClampSamples, 'omitnan'), ...
            'VariableNames', { ...
                'Animal', ...
                'SessionType', ...
                'TimePctClamp', ...
                'TimePctNonClamp', ...
                'DistanceClamp_pixels', ...
                'DistanceNonClamp_pixels', ...
                'DAMeanClamp', ...
                'DAVarianceClamp', ...
                'DAMeanNonClamp', ...
                'DAVarianceNonClamp', ...
                'NDAClampSamples', ...
                'NDANonClampSamples'});

        animalTypeSummaryTable = [animalTypeSummaryTable; oneRow]; %#ok<AGROW>
    end

    %% ---------- Per-animal figure ----------
    if options.plotAnimal
        plotAnimalSummary(sessions, animalName, xLimits, yLimits, ...
            clampColor, nonClampColor, ctrlColor, options);
    end

    %% ---------- Save per-animal session summary CSV ----------
    sessionSummaryTable = buildSessionSummaryTable(sessions, animalName);
    animalSummaryCSV = fullfile(animalFolder, sprintf('%s_%s_sessionSummary.csv', animalName, options.outputName));
    writetable(sessionSummaryTable, animalSummaryCSV);

    allSessionSummaryTable = [allSessionSummaryTable; sessionSummaryTable]; %#ok<AGROW>
end

%% ========================== GROUP SUMMARY ==========================

if isempty(allSessionSummaryTable)
    warning('No valid session data found. Group summary not generated.');
    return
end

if options.outputFolder == ""
    [groupRoot, ~] = fileparts(char(animalPaths(1)));
    groupOutputFolder = fullfile(groupRoot, 'GroupSummary');
else
    groupOutputFolder = char(options.outputFolder);
end

if ~exist(groupOutputFolder, 'dir')
    mkdir(groupOutputFolder);
end

allSessionCSV = fullfile(groupOutputFolder, sprintf('AcrossAnimals_%s_sessionSummary.csv', options.outputName));
writetable(allSessionSummaryTable, allSessionCSV);

groupSummaryCSV = fullfile(groupOutputFolder, sprintf('AcrossAnimals_%s_summary_bySessionType.csv', options.outputName));
writetable(animalTypeSummaryTable, groupSummaryCSV);

if options.plotGroup
    nAnimalsGroup = numel(unique(allSessionSummaryTable.Animal, 'stable'));
    plotGroupSummary(allSessionSummaryTable, animalTypeSummaryTable, heatmapAccumulator, ...
        xCenters, yCenters, xLimits, yLimits, nXBins, nYBins, ...
        bluePurpleRed, clampColor, nonClampColor, nAnimalsGroup, ...
        groupOutputFolder, options);
end

disp(['Saved clamp RTPP session summary: ' allSessionCSV]);
disp(['Saved clamp RTPP group summary CSV: ' groupSummaryCSV]);

end

%% ========================== LOCAL FUNCTIONS ==========================

function tf = hasTimesCsv(folderPath)
tf = isfolder(folderPath) && ~isempty(dir(fullfile(folderPath, 'times*.csv')));
end

function session = emptySessionStruct()
session = struct( ...
    'name', "", ...
    'path', "", ...
    'csvPath', "", ...
    'sessionType', "", ...
    'folderCondition', "", ...
    'durationMin', NaN, ...
    'clampDurationMin', NaN, ...
    'nonClampDurationMin', NaN, ...
    'nFramesRaw', NaN, ...
    'nFramesAnalyzed', NaN, ...
    'nClampFrames', NaN, ...
    'nNonClampFrames', NaN, ...
    'timePctClamp', NaN, ...
    'timePctNonClamp', NaN, ...
    'distClamp', NaN, ...
    'distNonClamp', NaN, ...
    'dataClean', table(), ...
    'da', emptyDAStats());
end

function daStats = emptyDAStats()
daStats = struct( ...
    'hasTimeseries', false, ...
    'hasDASignal', false, ...
    'timeseriesFile', "", ...
    'signalName', "", ...
    'signalSystem', "", ...
    'signalFs', NaN, ...
    'meanClamp', NaN, ...
    'varClamp', NaN, ...
    'meanNonClamp', NaN, ...
    'varNonClamp', NaN, ...
    'nClampSamples', NaN, ...
    'nNonClampSamples', NaN, ...
    'nClampToNonClampTransitions', NaN, ...
    'nNonClampToClampTransitions', NaN);
end

function [T, csvPath] = readClampTimesCsv(sessionPath, fallbackFs)
csvFiles = dir(fullfile(sessionPath, 'times*.csv'));
if isempty(csvFiles)
    error('No times*.csv found');
end

[~, order] = sort({csvFiles.name});
csvFiles = csvFiles(order);
csvPath = fullfile(csvFiles(1).folder, csvFiles(1).name);
raw = readtable(csvPath);

if width(raw) < 4
    error('Expected times*.csv to contain at least Date, X, Y, ClampON columns');
end

varNames = string(raw.Properties.VariableNames);
dateIdx = find(contains(lower(varNames), ["date", "time", "item1"]), 1);
if isempty(dateIdx); dateIdx = 1; end

clampIdx = find(contains(lower(varNames), "clamp"), 1);
if isempty(clampIdx); clampIdx = 4; end

xIdx = find(strcmpi(varNames, "X") | endsWith(varNames, "_X", 'IgnoreCase', true), 1);
yIdx = find(strcmpi(varNames, "Y") | endsWith(varNames, "_Y", 'IgnoreCase', true), 1);

if isempty(xIdx) || isempty(yIdx)
    isNumericVar = varfun(@(x) isnumeric(x) || islogical(x), raw, 'OutputFormat', 'uniform');
    numericIdx = find(isNumericVar);
    numericIdx = numericIdx(numericIdx ~= clampIdx & numericIdx ~= dateIdx);
    if numel(numericIdx) < 2
        error('Could not identify X/Y columns in %s', csvPath);
    end
    xIdx = numericIdx(1);
    yIdx = numericIdx(2);
end

timeSec = parseTimeSeconds(raw{:, dateIdx}, height(raw), fallbackFs);
X = double(raw{:, xIdx});
Y = double(raw{:, yIdx});
ClampON = parseLogicalColumn(raw{:, clampIdx});

validMask = ~isnan(timeSec) & isfinite(X) & isfinite(Y) & ~isnan(double(ClampON));
T = table(timeSec(validMask), timeSec(validMask)./60, X(validMask), Y(validMask), ClampON(validMask), ...
    'VariableNames', {'timeSec', 'timeMin', 'X', 'Y', 'ClampON'});
end

function timeSec = parseTimeSeconds(timeColumn, nRows, fallbackFs)
timeSec = nan(nRows, 1);

if isdatetime(timeColumn)
    t = timeColumn;
else
    t = NaT(nRows, 1);
    timeStrings = string(timeColumn);
    formats = [ ...
        "yyyy-MM-dd'T'HH:mm:ss.SSSSSSSXXX", ...
        "yyyy-MM-dd'T'HH:mm:ss.SSSSSSXXX", ...
        "yyyy-MM-dd'T'HH:mm:ss.SSSXXX", ...
        "yyyy-MM-dd'T'HH:mm:ss.SSSSSSS'Z'", ...
        "yyyy-MM-dd'T'HH:mm:ss.SSSSSS'Z'", ...
        "yyyy-MM-dd'T'HH:mm:ss.SSS'Z'", ...
        "yyyy-MM-dd'T'HH:mm:ss"];

    for fmt = formats
        try
            tCandidate = datetime(timeStrings, 'InputFormat', fmt, 'TimeZone', 'UTC');
            if all(~isnat(tCandidate))
                t = tCandidate;
                break
            end
        catch
        end
    end

    if all(isnat(t))
        try
            t = datetime(timeStrings, 'TimeZone', 'local');
        catch
        end
    end
end

if all(~isnat(t))
    if isempty(t.TimeZone)
        t.TimeZone = 'UTC';
    end
    timeSec = seconds(t - t(1));
else
    timeSec = (0:nRows-1)' ./ fallbackFs;
end
end

function logicalColumn = parseLogicalColumn(rawColumn)
if islogical(rawColumn)
    logicalColumn = rawColumn;
elseif isnumeric(rawColumn)
    logicalColumn = rawColumn > 0.5;
else
    values = lower(strtrim(string(rawColumn)));
    logicalColumn = ismember(values, ["true", "t", "1", "yes", "y", "on", "clamp"]);
end

logicalColumn = logicalColumn(:);
end

function [sessionType, folderCondition] = parseClampSessionFolderName(sessionName)
parts = split(string(sessionName), ["-", "_"]);
parts = parts(parts ~= "");

if numel(parts) >= 3 && isDateLike(parts(1))
    sessionType = parts(2);
    folderCondition = parts(end);
elseif numel(parts) >= 2 && isDateLike(parts(1))
    sessionType = parts(2);
    folderCondition = "";
elseif numel(parts) >= 2
    sessionType = parts(1);
    folderCondition = parts(end);
else
    sessionType = string(sessionName);
    folderCondition = "";
end
end

function tf = isDateLike(token)
token = char(token);
tf = numel(token) == 8 && all(isstrprop(token, 'digit'));
end

function accumulator = addHeatmap(accumulator, sessionType, stateName, H)
key = makeHeatmapKey(sessionType, stateName);
if isKey(accumulator, key)
    accumulator(key) = accumulator(key) + H;
else
    accumulator(key) = H;
end
end

function key = makeHeatmapKey(sessionType, stateName)
key = sprintf('%s__%s', char(sessionType), char(stateName));
end

function H = getHeatmapOrZero(heatmapAccumulator, sessionType, stateName, nYBins, nXBins)
key = makeHeatmapKey(sessionType, stateName);
if isKey(heatmapAccumulator, key)
    H = heatmapAccumulator(key);
else
    H = zeros(nYBins, nXBins);
end
end

function H = normalizeHeatmap(H)
totalCount = sum(H(:));
if totalCount <= 0
    return
end
H = H ./ totalCount;
end

function value = weightedMean(values, weights)
values = values(:);
weights = weights(:);

validMask = ~isnan(values) & ~isnan(weights) & (weights > 0);

if ~any(validMask)
    value = NaN;
    return
end

value = sum(values(validMask) .* weights(validMask)) / sum(weights(validMask));
end

function dist = safeTrajectoryDistance(X, Y, mask)
if isempty(mask) || nnz(mask) < 2
    dist = 0;
else
    dist = getTrajectoryDistance(X, Y, filter = mask);
end
end

function daStats = analyzeDopamineForClampSession(sessionPath, sessionName, T, clampON, options, clampColor, nonClampColor)
daStats = emptyDAStats();

tsFile = findMatchingMatFile(sessionPath, "timeseries", sessionName);
if tsFile == ""
    return
end

daStats.hasTimeseries = true;
daStats.timeseriesFile = tsFile;

loaded = load(tsFile, 'timeSeries');
if ~isfield(loaded, 'timeSeries') || isempty(loaded.timeSeries)
    warning('timeseries file found but no timeSeries variable: %s', char(tsFile));
    return
end

[signal, signalFs, signalSystem, signalName] = getDASignal(loaded.timeSeries, options.daSignalName);
if isempty(signal)
    warning('No DA signal found in %s', char(tsFile));
    return
end

daStats.hasDASignal = true;
daStats.signalName = signalName;
daStats.signalSystem = signalSystem;
daStats.signalFs = signalFs;

[cameraTimeSec, daTimeSec] = getDATiming(sessionPath, sessionName, T, signal, signalFs, signalSystem);
validCamera = isfinite(cameraTimeSec) & cameraTimeSec >= 0 & cameraTimeSec <= daTimeSec(end);

if nnz(validCamera) < 2
    warning('Could not align ClampON to DA timing for %s', char(sessionName));
    return
end

[cameraTimeForInterp, uniqueIdx] = unique(cameraTimeSec(validCamera), 'stable');
clampForInterp = double(clampON(validCamera));
clampForInterp = clampForInterp(uniqueIdx);
[cameraTimeForInterp, sortIdx] = sort(cameraTimeForInterp);
clampForInterp = clampForInterp(sortIdx);

if numel(cameraTimeForInterp) < 2
    warning('Could not align ClampON to DA timing for %s', char(sessionName));
    return
end

clampAtDA = interp1(cameraTimeForInterp, clampForInterp, daTimeSec, 'nearest', NaN);
validDA = ~isnan(clampAtDA) & ~isnan(signal(:));
clampMaskDA = clampAtDA > 0.5;

DAClamp = signal(validDA & clampMaskDA);
DANonClamp = signal(validDA & ~clampMaskDA);

daStats.meanClamp = mean(DAClamp, 'omitnan');
daStats.varClamp = var(DAClamp, 0, 'omitnan');
daStats.meanNonClamp = mean(DANonClamp, 'omitnan');
daStats.varNonClamp = var(DANonClamp, 0, 'omitnan');
daStats.nClampSamples = numel(DAClamp);
daStats.nNonClampSamples = numel(DANonClamp);

clampToNonClampIdx = find(diff(clampON(:)) < 0) + 1;
nonClampToClampIdx = find(diff(clampON(:)) > 0) + 1;
clampToNonClampEvents = cameraTimeSec(clampToNonClampIdx);
nonClampToClampEvents = cameraTimeSec(nonClampToClampIdx);
clampToNonClampEvents = clampToNonClampEvents(isfinite(clampToNonClampEvents) & clampToNonClampEvents >= 0);
nonClampToClampEvents = nonClampToClampEvents(isfinite(nonClampToClampEvents) & nonClampToClampEvents >= 0);

[clampToNonClampTraces, timestamp] = extractAlignedTracesByTime( ...
    signal, signalFs, clampToNonClampEvents, options.daTimeRange);
[nonClampToClampTraces, ~] = extractAlignedTracesByTime( ...
    signal, signalFs, nonClampToClampEvents, options.daTimeRange);

daStats.nClampToNonClampTransitions = size(clampToNonClampTraces, 1);
daStats.nNonClampToClampTransitions = size(nonClampToClampTraces, 1);

if options.plotDA
    plotDATransitions(sessionPath, sessionName, signalName, timestamp, ...
        clampToNonClampTraces, nonClampToClampTraces, clampColor, nonClampColor, options);
end
end

function matFile = findMatchingMatFile(sessionPath, prefix, sessionName)
files = dir(fullfile(sessionPath, sprintf('%s_*.mat', prefix)));
if isempty(files)
    matFile = "";
    return
end

names = string({files.name});
matchIdx = find(contains(names, string(sessionName)), 1);
if isempty(matchIdx)
    matchIdx = 1;
end
matFile = string(fullfile(files(matchIdx).folder, files(matchIdx).name));
end

function [signal, signalFs, signalSystem, signalName] = getDASignal(timeSeries, requestedName)
signal = [];
signalFs = NaN;
signalSystem = "";
signalName = "";

names = string({timeSeries.name});
systems = strings(size(names));
if isfield(timeSeries, 'system')
    systems = string({timeSeries.system});
end

idx = find(strcmpi(names, requestedName), 1);

if isempty(idx)
    idx = find(contains(names, "dLight", 'IgnoreCase', true), 1);
end

if isempty(idx)
    idx = find(contains(names, ["DA", "dopamine"], 'IgnoreCase', true), 1);
end

if isempty(idx)
    photometryIdx = find(contains(systems, ["NI", "LJ", "labjack"], 'IgnoreCase', true));
    nonClampIdx = photometryIdx(~contains(names(photometryIdx), ["clamp", "laser"], 'IgnoreCase', true));
    if ~isempty(nonClampIdx)
        idx = nonClampIdx(1);
    end
end

if isempty(idx)
    return
end

signal = double(timeSeries(idx).data(:));
if isfield(timeSeries, 'finalFs') && ~isempty(timeSeries(idx).finalFs)
    signalFs = timeSeries(idx).finalFs;
else
    signalFs = 50;
    warning('timeSeries.%s has no finalFs; using signalFs = 50 Hz.', char(names(idx)));
end
if isfield(timeSeries, 'system') && ~isempty(timeSeries(idx).system)
    signalSystem = string(timeSeries(idx).system);
else
    signalSystem = "lj";
end
signalName = string(timeSeries(idx).name);
end

function [cameraTimeSec, daTimeSec] = getDATiming(sessionPath, sessionName, T, signal, signalFs, signalSystem)
cameraTimeSec = T.timeSec;
daTimeSec = (0:numel(signal)-1)' ./ signalFs;

syncFile = findMatchingMatFile(sessionPath, "sync", sessionName);
if syncFile == ""
    return
end

loaded = load(syncFile, 'params');
if ~isfield(loaded, 'params') || ~isfield(loaded.params, 'sync')
    return
end

params = loaded.params;
if ~isfield(params.sync, 'timeCamera') || isempty(params.sync.timeCamera)
    return
end

targetTime = getSyncTimeForSystem(params, signalSystem);
if isempty(targetTime)
    targetStart = 0;
else
    targetStart = targetTime(1);
end

timeCamera = params.sync.timeCamera(:);
if numel(timeCamera) == height(T)
    cameraTimeSec = timeCamera - targetStart;
elseif numel(timeCamera) > 1
    cameraTimeSec = interp1(linspace(1, height(T), numel(timeCamera)), ...
        timeCamera, (1:height(T))', 'linear', 'extrap') - targetStart;
end
end

function targetTime = getSyncTimeForSystem(params, systemName)
targetTime = [];
if ~isfield(params, 'sync')
    return
end

if strcmpi(systemName, "lj") || strcmpi(systemName, "labjack")
    if isfield(params.sync, 'timePhotometry')
        targetTime = params.sync.timePhotometry(:);
    end
elseif strcmpi(systemName, "ni")
    if isfield(params.sync, 'timeNI')
        targetTime = params.sync.timeNI(:);
    end
elseif contains(systemName, "cam", 'IgnoreCase', true)
    if isfield(params.sync, 'timeCamera')
        targetTime = params.sync.timeCamera(:);
    end
end
end

function [traces, timestamp] = extractAlignedTracesByTime(signal, signalFs, eventTimesSec, timeRange)
timestamp = timeRange(1):(1/signalFs):timeRange(2);
traces = nan(numel(eventTimesSec), numel(timestamp));

for eventIdx = 1:numel(eventTimesSec)
    firstBin = round((eventTimesSec(eventIdx) + timeRange(1)) * signalFs) + 1;
    lastBin = firstBin + numel(timestamp) - 1;

    if firstBin < 1 || lastBin > numel(signal)
        continue
    end

    traces(eventIdx, :) = signal(firstBin:lastBin);
end

traces = traces(~any(isnan(traces), 2), :);
end

function plotDATransitions(sessionPath, sessionName, signalName, timestamp, ...
    clampToNonClampTraces, nonClampToClampTraces, clampColor, nonClampColor, options)

initializeFig(0.67, 0.45);
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile; hold on
if isempty(clampToNonClampTraces)
    text(0.5, 0.5, 'No valid transitions', 'Units', 'normalized', 'HorizontalAlignment', 'center');
else
    plotSEM(timestamp, clampToNonClampTraces, nonClampColor, label = "Clamp -> non-clamp");
end
xline(0, '--', 'Color', [.3 .3 .3], 'HandleVisibility', 'off');
xlabel('Time from transition (s)');
ylabel(sprintf('%s DA', char(signalName)), 'Interpreter', 'none');
title('Clamp -> non-clamp');

nexttile; hold on
if isempty(nonClampToClampTraces)
    text(0.5, 0.5, 'No valid transitions', 'Units', 'normalized', 'HorizontalAlignment', 'center');
else
    plotSEM(timestamp, nonClampToClampTraces, clampColor, label = "Non-clamp -> clamp");
end
xline(0, '--', 'Color', [.3 .3 .3], 'HandleVisibility', 'off');
xlabel('Time from transition (s)');
ylabel(sprintf('%s DA', char(signalName)), 'Interpreter', 'none');
title('Non-clamp -> clamp');

sgtitle(sprintf('%s | %s', char(sessionName), char(signalName)), 'Interpreter', 'none');
saveFigures(gcf, sprintf('DA_clampTransitions_%s', char(signalName)), sessionPath, ...
    savePDF = options.savePDF, savePNG = options.savePNG);

if options.closeFigures
    close(gcf);
end
end

function plotAnimalSummary(sessions, animalName, xLimits, yLimits, clampColor, nonClampColor, ctrlColor, options)
nSessions = numel(sessions);
animalFolder = char(fileparts(sessions(1).path));

initializeFig(0.9, 1);
tl = tiledlayout(2, nSessions + 1, 'TileSpacing', 'compact', 'Padding', 'compact');

for sessionIdx = 1:nSessions
    X = sessions(sessionIdx).dataClean.X;
    Y = sessions(sessionIdx).dataClean.Y;
    clampON = sessions(sessionIdx).dataClean.ClampON;

    nexttile([2 1]); hold on
    plot(X(~clampON), Y(~clampON), 'Color', nonClampColor, 'LineWidth', 2);
    plot(X(clampON), Y(clampON), 'Color', clampColor, 'LineWidth', 2);
    legend({'Non-clamp', 'Clamp'}, 'Location', 'northeast');

    xlim(xLimits); ylim(yLimits);
    xlabel('X Position'); ylabel('Y Position');
    title(sprintf('%s | %s', sessions(sessionIdx).sessionType, sessions(sessionIdx).name), ...
        'Interpreter', 'none');
end

nexttile; hold on
for sessionIdx = 1:nSessions
    plotScatterBar(sessionIdx, sessions(sessionIdx).timePctClamp, ...
        color = clampColor, style = 'bar');
end
xticks(1:nSessions);
xticklabels(string({sessions.name}));
xtickangle(25);
ylabel('Time in clamp region (%)');
title('Clamp occupancy');

nexttile; hold on
for sessionIdx = 1:nSessions
    plotScatterBar(2*sessionIdx - 0.25, sessions(sessionIdx).distNonClamp, ...
        color = ctrlColor, style = 'bar');
    plotScatterBar(2*sessionIdx + 0.25, sessions(sessionIdx).distClamp, ...
        color = clampColor, style = 'bar');
end
xticks((1:nSessions) * 2);
xticklabels(string({sessions.name}));
xtickangle(25);
ylabel('Distance traveled (pixels)');
title('Distance by clamp state');
legend({'Non-clamp', 'Clamp'}, 'Location', 'best');

title(tl, sprintf('%s - clamp RTPP summary', animalName), 'Interpreter', 'none');
saveFigures(gcf, sprintf('%s_%s_summary', animalName, options.outputName), animalFolder, ...
    savePDF = options.savePDF, savePNG = options.savePNG);

if options.closeFigures
    close(gcf);
end
end

function sessionSummaryTable = buildSessionSummaryTable(sessions, animalName)
nSessions = numel(sessions);
animalColumn = repmat(string(animalName), nSessions, 1);

sessionSummaryTable = table( ...
    animalColumn, ...
    string({sessions.name})', ...
    string({sessions.sessionType})', ...
    string({sessions.folderCondition})', ...
    string({sessions.csvPath})', ...
    [sessions.durationMin]', ...
    [sessions.clampDurationMin]', ...
    [sessions.nonClampDurationMin]', ...
    [sessions.nFramesRaw]', ...
    [sessions.nFramesAnalyzed]', ...
    [sessions.nClampFrames]', ...
    [sessions.nNonClampFrames]', ...
    [sessions.timePctClamp]', ...
    [sessions.timePctNonClamp]', ...
    [sessions.distClamp]', ...
    [sessions.distNonClamp]', ...
    arrayfun(@(s) string(s.da.timeseriesFile), sessions)', ...
    arrayfun(@(s) string(s.da.signalName), sessions)', ...
    arrayfun(@(s) string(s.da.signalSystem), sessions)', ...
    arrayfun(@(s) s.da.signalFs, sessions)', ...
    arrayfun(@(s) s.da.meanClamp, sessions)', ...
    arrayfun(@(s) s.da.varClamp, sessions)', ...
    arrayfun(@(s) s.da.meanNonClamp, sessions)', ...
    arrayfun(@(s) s.da.varNonClamp, sessions)', ...
    arrayfun(@(s) s.da.nClampSamples, sessions)', ...
    arrayfun(@(s) s.da.nNonClampSamples, sessions)', ...
    arrayfun(@(s) s.da.nClampToNonClampTransitions, sessions)', ...
    arrayfun(@(s) s.da.nNonClampToClampTransitions, sessions)', ...
    'VariableNames', { ...
        'Animal', ...
        'SessionName', ...
        'SessionType', ...
        'FolderCondition', ...
        'TimesCsvPath', ...
        'SessionDurationMin', ...
        'ClampDurationMin', ...
        'NonClampDurationMin', ...
        'RawFrameCount', ...
        'AnalyzedFrameCount', ...
        'ClampFrameCount', ...
        'NonClampFrameCount', ...
        'TimePctClamp', ...
        'TimePctNonClamp', ...
        'DistanceClamp_pixels', ...
        'DistanceNonClamp_pixels', ...
        'TimeseriesFile', ...
        'DASignalName', ...
        'DASignalSystem', ...
        'DASignalFs', ...
        'DAMeanClamp', ...
        'DAVarianceClamp', ...
        'DAMeanNonClamp', ...
        'DAVarianceNonClamp', ...
        'NDAClampSamples', ...
        'NDANonClampSamples', ...
        'NDAClampToNonClampTransitions', ...
        'NDANonClampToClampTransitions'});
end

function plotGroupSummary(allSessionSummaryTable, animalTypeSummaryTable, heatmapAccumulator, ...
    xCenters, yCenters, xLimits, yLimits, nXBins, nYBins, ...
    bluePurpleRed, clampColor, nonClampColor, nAnimals, groupOutputFolder, options)

sessionTypeList = unique(allSessionSummaryTable.SessionType, 'stable');
nSessionTypes = numel(sessionTypeList);

initializeFig(1, 1);
tl = tiledlayout(2, nSessionTypes + 1, 'TileSpacing', 'compact', 'Padding', 'compact');

%% ---------- Top row: non-clamp heatmaps ----------
for typeIdx = 1:nSessionTypes
    currentType = sessionTypeList(typeIdx);
    H = getHeatmapOrZero(heatmapAccumulator, currentType, "NonClamp", nYBins, nXBins);
    H = normalizeHeatmap(H);

    nexttile; hold on
    imagesc(xCenters, yCenters, H);
    axis xy;
    xlim(xLimits); ylim(yLimits);
    xlabel('X'); ylabel('Y');
    title(sprintf('%s - non-clamp', currentType), 'Interpreter', 'none');
    colormap(bluePurpleRed);
    colorbar;
end

%% ---------- Top-right: clamp occupancy ----------
nexttile; hold on
for typeIdx = 1:nSessionTypes
    currentType = sessionTypeList(typeIdx);
    rowMask = animalTypeSummaryTable.SessionType == currentType;
    values = animalTypeSummaryTable.TimePctClamp(rowMask);
    values = values(~isnan(values));
    if isempty(values); continue; end
    plotScatterBar(typeIdx, values, color = clampColor, style = 'bar');
end
xticks(1:nSessionTypes);
xticklabels(sessionTypeList);
xtickangle(20);
ylabel('Time in clamp region (%)');
title('Across animals');

%% ---------- Bottom row: clamp heatmaps ----------
for typeIdx = 1:nSessionTypes
    currentType = sessionTypeList(typeIdx);
    H = getHeatmapOrZero(heatmapAccumulator, currentType, "Clamp", nYBins, nXBins);
    H = normalizeHeatmap(H);

    nexttile; hold on
    imagesc(xCenters, yCenters, H);
    axis xy;
    xlim(xLimits); ylim(yLimits);
    xlabel('X'); ylabel('Y');
    title(sprintf('%s - clamp', currentType), 'Interpreter', 'none');
    colormap(bluePurpleRed);
    colorbar;
end

%% ---------- Bottom-right: distance plot by session type ----------
nexttile; hold on
distGap = 0.14;

for typeIdx = 1:nSessionTypes
    currentType = sessionTypeList(typeIdx);
    rowMask = animalTypeSummaryTable.SessionType == currentType;

    nonClampDist = animalTypeSummaryTable.DistanceNonClamp_pixels(rowMask);
    nonClampDist = nonClampDist(~isnan(nonClampDist));

    clampDist = animalTypeSummaryTable.DistanceClamp_pixels(rowMask);
    clampDist = clampDist(~isnan(clampDist));

    xNonClamp = typeIdx - distGap;
    xClamp = typeIdx + distGap;

    if ~isempty(nonClampDist)
        plotScatterBar(xNonClamp, nonClampDist, color = nonClampColor, style = 'bar');
    end
    if ~isempty(clampDist)
        plotScatterBar(xClamp, clampDist, color = clampColor, style = 'bar');
    end

    if ~isempty(nonClampDist) && ~isempty(clampDist)
        plotStats(nonClampDist, clampDist, [xNonClamp, xClamp], testType = 'kstest');
    end
end

xticks(1:nSessionTypes);
xticklabels(sessionTypeList);
xtickangle(20);
ylabel('Distance traveled (pixels)');
title('Distance by clamp state');
legend({'Non-clamp', 'Clamp'}, 'Location', 'best');

title(tl, sprintf('Clamp RTPP Group Summary (n = %d animals)', nAnimals));
saveFigures(gcf, sprintf('AcrossAnimals_%s_summary_bySessionType', options.outputName), groupOutputFolder, ...
    savePDF = options.savePDF, savePNG = options.savePNG);

if options.closeFigures
    close(gcf);
end

end
