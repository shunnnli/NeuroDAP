% Shun_analyzeRTPP

% 2025/01/31

% Plots trajectory of RTPP centroid and performs other analysis

%% Load the CSV file

clear; close all;
addpath(genpath(osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Analysis/NeuroDAP/Methods')));
[twoColors,~,~,~,~,~,bluePurpleRed] = loadColors;

%% Pick one or more animal folders
animalPaths = uipickfiles('FilterSpec', osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Project misc/Recordings'), ...
                          'Prompt', 'Select one or MORE animal folders (e.g., SL399, SL400, ...)');
if isempty(animalPaths); return; end
if ischar(animalPaths); animalPaths = {animalPaths}; end
nAnimals = numel(animalPaths);

%% Params (edit as needed)
box           = 'large';                 % 'large' or 'small'
defaultStim   = 'right';                 % used for Baseline (so % = time on right)
removeStatic  = false;
noMoveThresh  = 20;                      % for getStaticPeriod
Fs            = 20;                      % frame rate (not used below, kept for consistency)

% arena and binning
switch lower(box)
    case 'large'
        Y_midpoint = 350;
        xlimit = [300,650]; ylimit = [0,700];
        nbx = 30; nby = 60;              % heatmap bins
    otherwise
        Y_midpoint = 140;
        xlimit = [420,520]; ylimit = [15,280];
        nbx = 16; nby = 20;
end
xedges = linspace(xlimit(1), xlimit(2), nbx+1);
yedges = linspace(ylimit(1), ylimit(2), nby+1);

% colors
leftColor = [156, 219, 17]./255;
rightColor = [144, 126, 171]./255;
stimColor  = [7, 162, 222]./255;
ctrlColor  = [.3 .3 .3];

%% Group accumulators
grpHeat.Baseline = zeros(nby, nbx);   % histcounts2 returns [ny, nx]
grpHeat.Left     = zeros(nby, nbx);
grpHeat.Right    = zeros(nby, nbx);

% we’ll collect per-animal means (so each animal contributes one point per session-type)
grpStimPct.Baseline = []; grpStimPct.Left = []; grpStimPct.Right = [];
grpDistStim.Baseline = []; grpDistStim.Left = []; grpDistStim.Right = [];
grpDistCtrl.Baseline = []; grpDistCtrl.Left = []; grpDistCtrl.Right = [];

%% ---------- PER-ANIMAL LOOP ----------
for a = 1:nAnimals
    animalDir = animalPaths{a};
    [parentDir, animalName] = fileparts(animalDir);

    % session folders in this animal folder
    raw = dir(animalDir);
    sessionList = raw([raw.isdir]);
    sessionList = sessionList(~ismember({sessionList.name},{'.','..'}));
    nSessions = numel(sessionList);

    if nSessions==0
        warning('No session subfolders for animal: %s', animalName);
        continue;
    end

    % holders
    sessions = struct([]);
    stim_pct = nan(nSessions,1);
    side_dist = nan(nSessions,2); % [stim ctrl]
    sessTypes = strings(nSessions,1); % 'Baseline'|'Left'|'Right'

    % --------- load & preprocess each session ----------
    for s = 1:nSessions
        cur_session = dir(fullfile(sessionList(s).folder, sessionList(s).name));
        % find CSV named like times-*.csv
        hit = contains({cur_session.name}, 'times-') & endsWith({cur_session.name}, '.csv');
        assert(any(hit), 'No times-*.csv in %s', fullfile(sessionList(s).folder, sessionList(s).name));
        tableName = cur_session(hit).name;
        tablePath = cur_session(hit).folder;

        T = readtable(fullfile(tablePath, tableName));

        % session meta
        dirsplit = split(tablePath, filesep);
        sessName = dirsplit{end};
        sessions(s).name = sessName;
        sessions(s).path = tablePath;

        lowname = lower(sessName);
        if contains(lowname, 'left')
            sessTypes(s) = "Left";  sessions(s).stimSide = 'left';
        elseif contains(lowname, 'right')
            sessTypes(s) = "Right"; sessions(s).stimSide = 'right';
        elseif contains(lowname, 'base')
            sessTypes(s) = "Baseline"; sessions(s).stimSide = defaultStim; % use default
        else
            % fallback (treat as baseline-like)
            sessTypes(s) = "Baseline"; sessions(s).stimSide = defaultStim;
        end
        sessions(s).box = box;

        % add time column (align to first frame)
        T.Item1 = datetime(T.Item1, 'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSSSSSSZ','TimeZone','UTC');
        t0 = T.Item1(1);
        minuteSinceStart = minutes(T.Item1 - t0);
        T.time = minuteSinceStart; % no cutoffs across animals; edit if needed

        % keep only t>=0 (safety) and drop static if needed
        cur_data = T(T.time >= 0, :);
        X_raw = cur_data.Item2_X; Y_raw = cur_data.Item2_Y;

        if removeStatic
            staticMask = getStaticPeriod(X_raw, Y_raw, noMovementThreshold=noMoveThresh, windowDuration=30);
            keepMask = ~staticMask;
        else
            keepMask = true(size(X_raw));   % keep everything
        end
        
        X = X_raw(keepMask);
        Y = Y_raw(keepMask);


        % --- metrics (time & distance per chamber) ---
        stimSide = sessions(s).stimSide;
        if strcmpi(stimSide,'right')
            stimIdx = find(Y >= Y_midpoint);
            ctrlIdx = find(Y <  Y_midpoint);
        else % 'left'
            stimIdx = find(Y <  Y_midpoint);
            ctrlIdx = find(Y >= Y_midpoint);
        end

        side_dist(s,1) = getTrajectoryDistance(X, Y, filter=stimIdx);
        side_dist(s,2) = getTrajectoryDistance(X, Y, filter=ctrlIdx);
        stim_pct(s)    = numel(stimIdx) / max(1,numel(Y)) * 100;

        % --- accumulate session heatmap into the right bucket ---
        H = histcounts2(Y, X, yedges, xedges);  % note order (y,x)
        switch sessTypes(s)
            case "Baseline", grpHeat.Baseline = grpHeat.Baseline + H;
            case "Left",     grpHeat.Left     = grpHeat.Left     + H;
            case "Right",    grpHeat.Right    = grpHeat.Right    + H;
        end

        % stash the cleaned data for plotting the per-animal figure
        sessions(s).dataClean = table(cur_data.time(keepMask), X, Y, ...
                                      'VariableNames', {'time','X','Y'});
    end

    %% --------- PER-ANIMAL FIGURE (and save) ----------
    initializeFig(0.7,1); tl = tiledlayout(2, nSessions+1, 'TileSpacing','compact', 'Padding','compact');

    for s = 1:nSessions
        X = sessions(s).dataClean.X; Y = sessions(s).dataClean.Y;
        nexttile([2 1]); hold on;

        % color by side
        if sessTypes(s)=="Baseline"
            plot(X(Y>=Y_midpoint), Y(Y>=Y_midpoint), 'Color', rightColor, 'LineWidth', 2);
            plot(X(Y< Y_midpoint), Y(Y< Y_midpoint), 'Color', leftColor,  'LineWidth', 2);
            legend({'Right','Left'}, 'Location','northeast');
        else
            if strcmpi(sessions(s).stimSide,'right')
                plot(X(Y>=Y_midpoint), Y(Y>=Y_midpoint), 'Color', stimColor, 'LineWidth', 2);
                plot(X(Y< Y_midpoint), Y(Y< Y_midpoint), 'Color', ctrlColor, 'LineWidth', 2);
            else
                plot(X(Y>=Y_midpoint), Y(Y>=Y_midpoint), 'Color', ctrlColor, 'LineWidth', 2);
                plot(X(Y< Y_midpoint), Y(Y< Y_midpoint), 'Color', stimColor, 'LineWidth', 2);
            end
            legend({'Stim OFF','Stim ON'}, 'Location','northeast');
        end
        xlim(xlimit); ylim(ylimit);
        xlabel('X Position'); ylabel('Y Position');
        title(sprintf('%s', sessions(s).name), 'Interpreter','none');
    end

    % Right-top: time spent in stimulated side (%)
    nexttile; hold on;
    for s = 1:nSessions
        thisColor = ctrlColor;               % default for baseline bar color
        if sessTypes(s) ~= "Baseline", thisColor = stimColor; end
        plotScatterBar(s, stim_pct(s), Color=thisColor, style='bar');
    end
    xticks(1:nSessions); xticklabels({sessions.name}); xtickangle(20);
    ylabel('Time spent in stimulated side (%)');

    % Right-bottom: distance traveled stim vs ctrl
    nexttile; hold on;
    for s = 1:nSessions
        plotScatterBar(2*s-0.5, side_dist(s,1), Color=stimColor, style='bar');
        plotScatterBar(2*s+0.5, side_dist(s,2), Color=ctrlColor, style='bar');
    end
    xticks((1:nSessions)*2); xticklabels({sessions.name}); xtickangle(20);
    ylabel('Distance traveled (pixels)');

    title(tl, sprintf('%s — Summary', animalName), 'Interpreter','none');

    % save in the animal folder
    outPNG = fullfile(animalDir, sprintf('%s_summary.png', animalName));
    try
        exportgraphics(gcf, outPNG, 'Resolution', 300);
    catch
        saveas(gcf, outPNG);
    end
    close(gcf);

    %% --------- push this animal into GROUP arrays (per-type means) ----------
    % Per animal mean across sessions of the same type
    for typ = ["Baseline","Left","Right"]
        mStim = mean(stim_pct(sessTypes==typ), 'omitnan');
        mDstS = mean(side_dist(sessTypes==typ,1), 'omitnan');
        mDstC = mean(side_dist(sessTypes==typ,2), 'omitnan');

        if ~isnan(mStim), eval(sprintf('grpStimPct.%s(end+1) = mStim;', typ)); end
        if ~isnan(mDstS),  eval(sprintf('grpDistStim.%s(end+1) = mDstS;', typ)); end
        if ~isnan(mDstC),  eval(sprintf('grpDistCtrl.%s(end+1) = mDstC;', typ)); end
    end
end

%% ---------- GROUP SUMMARY FIGURE ----------
% normalize heatmaps to probability (each type independently)
normHeat = @(H) H ./ max(1, sum(H(:)));
HB = normHeat(grpHeat.Baseline);
HL = normHeat(grpHeat.Left);
HR = normHeat(grpHeat.Right);

% where to save
[groupRoot, ~] = fileparts(animalPaths{1});
groupOutDir = fullfile(groupRoot, 'GroupSummary');
if ~exist(groupOutDir,'dir'), mkdir(groupOutDir); end

initializeFig(0.7, 1);
tl = tiledlayout(2,nSessions+1, 'TileSpacing','compact','Padding','compact');

% three summary heatmaps (span both rows)
nexttile([2 1]); imagesc(xedges, yedges, HB); axis xy; xlim(xlimit); ylim(ylimit);
title('Baseline occupancy'); xlabel('X'); ylabel('Y'); colormap(sky); colorbar;

nexttile([2 1]); imagesc(xedges, yedges, HL); axis xy; xlim(xlimit); ylim(ylimit);
title('Left (stim on left) occupancy'); xlabel('X'); ylabel('Y'); colormap(sky); colorbar;

nexttile([2 1]); imagesc(xedges, yedges, HR); axis xy; xlim(xlimit); ylim(ylimit);
title('Right (stim on right) occupancy'); xlabel('X'); ylabel('Y'); colormap(sky); colorbar;

% right-top: time spent in stimulated side (%) across animals
types = ["Baseline","Left","Right"];
typeX  = 1:numel(types);
nexttile; hold on;
for i = 1:numel(types)
    vals = grpStimPct.(char(types(i)));
    if isempty(vals), continue; end
    baseColor = ctrlColor; 
    if types(i)~="Baseline", baseColor = stimColor; end
    plotScatterBar(2*typeX(i), vals, Color=baseColor, style='bar');

    % Plot significance
    for j = (i+1):numel(types)
        plotStats(vals, grpStimPct.(char(types(j))), [2*typeX(i), 2*typeX(j)], testType='kstest');
    end
end
xticks((1:length(types))*2); xticklabels(types);
ylabel('Time spent in stimulated side (%)');
title('Across animals');

% right-bottom: distance traveled — stim vs ctrl for each type (bars = mean; dots = animals)
nexttile; hold on;
gap = 0.3; % half-gap between stim and ctrl bars within each type
for i = 1:numel(types)
    vS = grpDistStim.(char(types(i)));
    vC = grpDistCtrl.(char(types(i)));
    if isempty(vS) || isempty(vC), continue; end
    plotScatterBar(2*typeX(i)-gap, vC, Color=ctrlColor, style='bar');
    plotScatterBar(2*typeX(i)+gap, vS, Color=stimColor, style='bar');
    plotStats(vC, vS, [2*typeX(i)+gap, 2*typeX(i)+gap], testType='kstest');
end
xticks((1:length(types))*2); xticklabels(types);
ylabel('Distance traveled (pixels)');
title('Across animals');

title(tl, sprintf('RTPP Group Summary  (n=%d animals)', nAnimals));

% save group figure
groupPNG = fullfile(groupOutDir, 'AcrossAnimals_summary.png');
try
    exportgraphics(gcf, groupPNG, 'Resolution', 300);
catch
    saveas(gcf, groupPNG);
end
disp(['Saved group summary: ' groupPNG]);
