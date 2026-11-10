% Shun_analyzeExperiments_template
% 2023/12/04

% Template for multi-session analysis of an experiment

% To start, copy this matlab file and replace template with specific
% experiments. In theory, there will be a specific analyzeExperiments file
% for each individual experiments, as the specific needs for analysis
% varies between different experiments.

%% Analysis pipeline
% The pipeline in general is the following:

% 1. Select whether to load a previously animals struct (described below)
% or select individual session to combine.

% 2. After selecting ALL SESSIONS from an experiments, the pipeline will
% automatically concatenate analysis.mat for each recording sessions.
% Rename properties as needed in order to facilitate further analysis.

% 3. Run getAnimalStruct.m function to recreate animals struct from
% summary. This combines all sessions from the same animals together while
% cutoffs between individual sessions are also recorded.

% 4. Save animals struct if needed. Note: saving summary struct will take
% extremely long (>5hrs) so while saving animals struct is much shorter 
% (~2min). animals struct should contain information that satisfies MOST 
% plotting requirements so saving summary struct is not needed.

% 5. Data analysis and plotting. This part is designed to vary across
% experiments. Thus, following codes are just for demonstration of 
% essential functions.

%% Essential functions

% getAnimalStruct(summary)
% combine sessions of the same animal, from the same task,
% of the same event, recorded from the same signal (eg NAc, LHb, cam,
% Lick) together. As described above, animals struct will be the MOST
% IMPORTANT struct that stores information about the experiments.

% combineTraces(animals,options)
% combine traces and their relevant statstics of selected animals, 
% selected tasks, selected trialRange, selected totalTrialRange, 
% selected signals, and selected events together. 
% This is the MOST IMPORTANT and USED function in this script. 
% Important features are listed as follows:
    % 1. the function returns a structure with fields. data fields stores
    % the data (photometry, cam, lick rate traces) of selected sessions.
    % 2. Field stats stores stageAvg/Max/Min of each traces at selected stage
    % time (often determined when creating analysis.mat but can modify later).
    % 3. Field options contains following important variables:
        % options.empty: true if no session is found that fits the input criteria. 
        % Should skip during plotting or further analysis 
        % options.animalStartIdx: Records index (in field data) of the
        % first trace for each animals. Used in plotGroupTraces
        % options.sessionStartIdx: Records index (in field data) of the
        % first trace for each session.
    % 4. totalTrialRange and trialRange
        % totalTrialRange selects the ACTUAL trial number within each
        % session while trialRange selects the samples across selected
        % sessions. For example: I have 3 session where I inhibit CaMKII
        % activity for the first 60 trials of each session. Within these
        % first 60 trials, 30% of them are stim-only trials. If I want to
        % only plot the 50-100th stim-only trials with CaMKII inhibition 
        % across all sessions, I will set totalTrialRange=[1,60] and
        % trialRange=[50,100]. Detailed description and automatic handling
        % of edge cases is documented within the method.
    % 5. The function can take both animals and summary struct as inputs.

% plotGroupTraces(combined.data,combined.timestamp,options)
% While plotGroupTraces is also used in analyzeSessions.m; here, we can
% plot traces across all animal easily (see code below). Key options are as
% follows:
    % 1. groupSize and nGroups
        % You need to provide either groupSize or nGroups for the function to
        % run. If you provide both, plotGroupTraces will plot to the maximum
        % number of groups based on groupSize. Thus, for a input with 50
        % trials and groupSize = 10, the function will automatically plot 5
        % lines even when nGroups=10
    % 2. options.animalStartIdx
        % Use this to reorganize input data so that its plotted based on
        % animals. eg when I want to plot Trial 1-10, 11-20 for EACH animal
        % across all sessions
    % 3. options.remaining
        % There inevitably will be some traces that does not fully form a group
        % (eg 5 traces remaining for a groupSize of 50 traces). These traces,
        % if plotted separately, can induce lines with great variations and
        % error bars. To address this, one can either set remaining='include'
        % to include these traces to prev group; set to 'exclude' to not plot
        % these traces, or 'separate' if you really want to plot these traces
        % separately

%% Intro of sample data set

% The sample data set is recorded by Shun Li in 2023. It contains 4
% animals, with 1 animals with off-target expression ('SL137'). dLight
% signals in NAc, pupil/Eye area, and lick are simultaneously recorded for
% all sessions.

% There are 5 major phases:
    % 1. Random: water, airpuff, EP stim, and tone (75dB) are delivered
    % randomly
    % 2. Reward1/2: where EP stim and tone are paired with water
    % 3. Punish1/2: where EP stim and tone are paired with airpuff
    % 4. Timeline: Random (2 sessions) -> Reward1 (3 sessions) -> Punish1
    % (3 sessions) -> Reward2 (3 sessions) -> 1 week rest -> Punish2 (3 sessions but 3 animals)

%% Setup

clear; close all;
loadNeuroDAP;
[~,~,~,~,~,~,bluePurpleRed] = loadColors;
clampColor = [.232 .76 .58];
unclampColor = [165, 209, 178]./255;

% Define result directory
resultspath = osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Project clamping/Results');

% Building summary struct from selected sessions
answer = questdlg('Group sessions or load combined data?','Select load sources',...
                  'Group single sessions','Load combined data','Load sample data','Load combined data');

if strcmpi(answer,'Group single sessions')
    sessionList = uipickfiles('FilterSpec',osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Project clamping/Recordings'))';
    groupSessions = true;
    % Update resultspath
    dirsplit = strsplit(sessionList{1},filesep); projectName = dirsplit{end-1}; 
    resultspath = strcat(resultspath,filesep,projectName);
    % Create resultspath if necessary
    if isempty(dir(resultspath)); mkdir(resultspath); end

elseif strcmpi(answer,'Load combined data')
    fileList = uipickfiles('FilterSpec',osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Project clamping/Results'))';
    groupSessions = false;
    % Update resultspath
    dirsplit = strsplit(fileList{1},filesep); projectName = dirsplit{end-1}; 
    resultspath = strcat(resultspath,filesep,projectName);
    % Load selected files
    for file = 1:length(fileList)
        dirsplit = strsplit(fileList{file},filesep);
        disp(['Ongoing: loading ',dirsplit{end}]);
        load(fileList{file});
        disp(['Finished: loaded ',dirsplit{end}]);
    end

elseif strcmpi(answer,'Load sample data')
    groupSessions = false;
    % Update resultspath
    dirsplit = strsplit(fileList{1},filesep); projectName = dirsplit{end-1}; 
    resultspath = strcat(osPathSwith('/Volumes/Neurobio/MICROSCOPE/Shun/Analysis/NeuroDAP/Tutorials/Sample data/Results'),filesep,projectName);
    % Load selected files
    for file = 1:length(fileList)
        dirsplit = strsplit(fileList{file},filesep);
        disp(['Ongoing: loading ',dirsplit{end}]);
        load(fileList{file});
        disp(['Finished: loaded ',dirsplit{end}]);
    end
end





%% 1. Generate clampQuality for the currently loaded animals.mat

clear currentClampQuality; % do not retain an older dataset if generation fails

% Clamp quality: what is the ratio of clamp transients vs unclamp ctrl
% Option1: RMS of clamp divided by unclamp RMS
% Option2: peak amplitude of response between clamp vs unclamp

% Water trials might be the dominating factor, which is reasonable!

% Plot the quality at the x axis, animal at y axis
qualitySignal = 'NAc-clamp'; % recorded signal on the clamped side
qualityEvents = {'Water','Airpuff','Tone'};
qualityBaseline = [-0.5,0]; % baseline includes start, excludes event onset
qualityWindow = [0,1.5]; % quality response window for Water, Airpuff, and Tone
qualityConditions = {'clamp','unclamp'};

% Compare performing trials in Random tasks, including RandomClampMix.
qualityRows = contains({animals.task},'Random',IgnoreCase=true) & ...
              strcmpi({animals.name},qualitySignal);
qualityAnimals = unique(string({animals(qualityRows).animal}));
qualityAnimals = qualityAnimals(:);
qualityValues = nan(numel(qualityAnimals),2,2); % animal x condition x metric
qualityTrialCount = zeros(numel(qualityAnimals),2);
qualitySummaryRows = cell(0,15);
clampQualityTrials = table('Size',[0,7], ...
    'VariableTypes',{'string','string','string','double','double','double','double'}, ...
    'VariableNames',{'Animal','Event','Condition','AnimalRow','TrialRow','RMS','Peak'});

for a = 1:numel(qualityAnimals)
    pooledMetrics = {zeros(0,2),zeros(0,2)};
    for e = 1:numel(qualityEvents)
        eventMetrics = {zeros(0,2),zeros(0,2)};
        for c = 1:2
            rows = find(qualityRows & ...
                strcmpi({animals.animal},qualityAnimals(a)) & ...
                strcmpi({animals.event},[qualityEvents{e},' (',qualityConditions{c},')']));
            for r = rows
                row = animals(r);
                t = row.timestamp(:)';
                baselineIdx = t >= qualityBaseline(1) & t < qualityBaseline(2);
                responseIdx = t >= qualityWindow(1) & t <= qualityWindow(2);
                if isempty(t) || t(1) > qualityBaseline(1) || ...
                        t(end) < qualityWindow(2) || ...
                        ~any(baselineIdx) || ~any(responseIdx)
                    warning('ClampQuality:Window', ...
                        'Skipping row %d: baseline/response window unavailable.',r);
                    continue;
                end

                % Keep the same performing-trial filter as the trace plots.
                keepTrials = row.trialInfo.trialTable.performing == 1;
                traces = row.data(keepTrials,:);
                % Require complete windows so RMS/peak use comparable data.
                validTrials = all(isfinite(traces(:,baselineIdx)),2) & ...
                              all(isfinite(traces(:,responseIdx)),2);
                trialRows = find(keepTrials);
                trialRows = trialRows(validTrials);
                traces = traces(validTrials,:);
                if isempty(traces); continue; end
                response = traces(:,responseIdx) - mean(traces(:,baselineIdx),2);
                trialRMS = sqrt(mean(response.^2,2));
                trialPeak = max(abs(response),[],2); % includes negative transients
                eventMetrics{c} = [eventMetrics{c}; trialRMS,trialPeak];
                n = numel(trialRMS);
                clampQualityTrials = [clampQualityTrials; table( ...
                    repmat(qualityAnimals(a),n,1),repmat(string(qualityEvents{e}),n,1), ...
                    repmat(string(qualityConditions{c}),n,1),repmat(r,n,1), ...
                    trialRows(:),trialRMS,trialPeak, ...
                    'VariableNames',clampQualityTrials.Properties.VariableNames)];
            end
        end
        % Save each event/condition even when its counterpart is missing.
        pairedEvent = ~isempty(eventMetrics{1}) && ~isempty(eventMetrics{2});
        for c = 1:2
            metrics = eventMetrics{c};
            n = size(metrics,1);
            metricMean = [NaN,NaN]; metricSD = [NaN,NaN]; metricSEM = [NaN,NaN];
            if n > 0; metricMean = mean(metrics,1); end
            if n > 1
                metricSD = std(metrics,0,1);
                metricSEM = metricSD / sqrt(n);
            end
            metricSum = sum(metrics,1);
            metricSumSquares = sum(metrics.^2,1);
            qualitySummaryRows(end+1,:) = {qualityAnimals(a),string(qualityEvents{e}), ...
                string(qualityConditions{c}),n,metricMean(1),metricSD(1),metricSEM(1), ...
                metricSum(1),metricSumSquares(1),metricMean(2),metricSD(2), ...
                metricSEM(2),metricSum(2),metricSumSquares(2),pairedEvent};
        end
        % Pool only event types with usable trials in BOTH conditions.
        if ~pairedEvent; continue; end
        for c = 1:2
            pooledMetrics{c} = [pooledMetrics{c}; eventMetrics{c}];
        end
    end
    for c = 1:2
        qualityTrialCount(a,c) = size(pooledMetrics{c},1);
        if isempty(pooledMetrics{c}); continue; end
        qualityValues(a,c,:) = reshape(mean(pooledMetrics{c},1),1,1,2);
    end
end

% Ratio of mean trial metrics, pooling events/trials (no event normalization).
% Thus large responses, e.g. Water, can dominate. Trial counts may differ
% between conditions; these are pooled comparisons, not paired trials.
qualityRatios = nan(numel(qualityAnimals),2);
for m = 1:2
    validControl = isfinite(qualityValues(:,2,m)) & qualityValues(:,2,m) > 0;
    qualityRatios(validControl,m) = qualityValues(validControl,1,m) ./ ...
                                   qualityValues(validControl,2,m);
end
clampQuality = table(qualityAnimals,qualityValues(:,1,1),qualityValues(:,2,1), ...
    qualityRatios(:,1),qualityValues(:,1,2),qualityValues(:,2,2), ...
    qualityRatios(:,2),qualityTrialCount(:,1),qualityTrialCount(:,2), ...
    'VariableNames',{'Animal','ClampRMS','UnclampRMS','RMSRatio', ...
    'ClampPeak','UnclampPeak','PeakRatio','NClamp','NUnclamp'});
disp(clampQuality);

% Long-format event summaries: one row per animal/event/condition.
% N, Sum and SumSquares allow correct pooling across exports; do not average
% saved ratios or SEMs. SD/SEM are undefined (NaN) for fewer than two trials.
clampQualityEventSummary = cell2table(qualitySummaryRows, ...
    'VariableNames',{'Animal','Event','Condition','N','RMSMean','RMSSD','RMSSEM', ...
    'RMSSum','RMSSumSquares','PeakMean','PeakSD','PeakSEM','PeakSum', ...
    'PeakSumSquares','IncludedInPooledComparison'});
clampQualityEventSummary.Animal = string(clampQualityEventSummary.Animal);
clampQualityEventSummary.Event = string(clampQualityEventSummary.Event);
clampQualityEventSummary.Condition = string(clampQualityEventSummary.Condition);

% Common count weights prevent different event proportions in clamp/control.
% Keep the original pooled table for comparison and save the weighted result.
[clampQualityWeighted,clampQualityWeights] = weightClampQuality( ...
    clampQualityEventSummary,qualityAnimals,string(qualityEvents));
disp('Trial-count-weighted clamp quality:'); disp(clampQualityWeighted);

clampQualitySettings = struct('version',1,'signal',qualitySignal, ...
    'events',{qualityEvents},'baselineWindow',qualityBaseline, ...
    'responseWindow',qualityWindow,'performingOnly',true, ...
    'requireCompleteWindows',true,'peakDefinition','max(abs(baseline-corrected trial))', ...
    'pooling','Ratio of mean trial metrics; only events present in both conditions', ...
    'createdAt',datetime('now'),'sourceFiles',strings(0,1));
clampQualitySettings.weighting = 'Common event weights = min(NClamp,NUnclamp); ratio of weighted means';
if exist('fileList','var') && (~exist('groupSessions','var') || ~groupSessions)
    clampQualitySettings.sourceFiles = string(fileList(:));
elseif exist('sessionList','var')
    clampQualitySettings.sourceFiles = string(sessionList(:));
end
% Keep the generated dataset separate from the saved, combined dataset.
currentClampQuality = struct('clampQuality',clampQuality, ...
    'clampQualityTrials',clampQualityTrials, ...
    'clampQualityEventSummary',clampQualityEventSummary, ...
    'clampQualityWeighted',clampQualityWeighted,'clampQualityWeights',clampQualityWeights, ...
    'clampQualitySettings',clampQualitySettings);

%% 2. Save/append current clampQuality to clampQuality.mat
% New animals are appended. Existing animals are REPLACED by their current
% results, not pooled across datasets, so rerunning does not duplicate trials.
% Run section 1 after loading each animals.mat, then run this section.
if ~exist('currentClampQuality','var')
    error('ClampQuality:GenerateFirst', ...
        'Run section 1 (Generate clampQuality) for the loaded animals.mat first.');
end
qualitySaveFolder = pwd;
if exist('resultspath','var') && isfolder(resultspath)
    qualitySaveFolder = resultspath;
end
qualityFile = fullfile(qualitySaveFolder,'clampQuality.mat');
savedClampQuality = appendClampQuality(qualityFile,currentClampQuality);
fprintf('Saved %d animals to %s\n',height(savedClampQuality.clampQuality),qualityFile);

%% 3. Plot clamp quality figure from saved clampQuality.mat
% This section can run independently of sections 1 and 2.
qualitySaveFolder = pwd;
if exist('resultspath','var') && isfolder(resultspath)
    qualitySaveFolder = resultspath;
end
qualityFile = fullfile(qualitySaveFolder,'clampQuality.mat');
if ~isfile(qualityFile)
    error('ClampQuality:NoSavedData','Run section 2 first: %s does not exist.',qualityFile);
end
qualityPlotData = load(qualityFile);
qualityAnimals = string(qualityPlotData.clampQuality.Animal);
qualitySignal = qualityPlotData.clampQualitySettings.signal;
qualityPlotColor = [.232 .76 .58];
qualityPooledDotSize = 260;
qualityEventDotSize = 150;
qualityPlotEvents = {'Water','Airpuff','Tone'};
[~,~,~,~,~,~,qualityEventMap] = loadColors;
% Match the event colors used in the Random trace plots above.
qualityEventColors = [qualityEventMap(1,:); .2,.2,.2; qualityEventMap(100,:)];
eventSummary = table();
if isfield(qualityPlotData,'clampQualityEventSummary')
    eventSummary = qualityPlotData.clampQualityEventSummary;
end
% Recompute from saved counts/sums so older exports with event summaries
% can be plotted with the new weighting without reloading full traces.
[qualityWeighted,qualityPlotWeights,qualityEventRatios] = weightClampQuality( ...
    eventSummary,qualityAnimals,string(qualityPlotEvents));
qualityRatios = [qualityWeighted.RMSRatio,qualityWeighted.PeakRatio];
% Display percent reduction; retain saved ratios for reuse.
qualityReduction = 100 * (1 - qualityRatios);
qualityEventReduction = 100 * (1 - qualityEventRatios);

if isempty(qualityAnimals)
    warning('ClampQuality:NoData','No animals in the saved clampQuality.mat.');
else
    initializeFig(.6,.5); tiledlayout(1,2);
    qualityTitles = {'Trial-weighted RMS reduction','Trial-weighted absolute peak reduction'};
    for m = 1:2
        nexttile; hold on;
        for a = 1:numel(qualityAnimals)
            if ~isfinite(qualityReduction(a,m)); continue; end
            plotScatterBar(a,qualityReduction(a,m),orientation='horizontal', ...
                style='bar',color=qualityPlotColor,dotSize=qualityPooledDotSize, ...
                connectPairs=false,XJitter='none');
        end
        % Event points share the exact animal row; only the combined reduction has a bar.
        for e = 1:numel(qualityPlotEvents)
            for a = 1:numel(qualityAnimals)
                if ~isfinite(qualityEventReduction(a,e,m)); continue; end
                plotScatterBar(a,qualityEventReduction(a,e,m),orientation='horizontal', ...
                    style='points',color=qualityEventColors(e,:),dotSize=qualityEventDotSize, ...
                    connectPairs=false,XJitter='none');
            end
        end
        xline(0,'--','No reduction','LabelVerticalAlignment','bottom');
        missing = find(~isfinite(qualityReduction(:,m)));
        for a = missing'
            text(0,a,'N/A','VerticalAlignment','middle');
        end
        yticks(1:numel(qualityAnimals)); yticklabels(qualityAnimals);
        set(gca,'YDir','reverse','TickLabelInterpreter','none');
        ylim([0.5,numel(qualityAnimals)+0.5]);
        panelReduction = [qualityReduction(:,m); reshape(qualityEventReduction(:,:,m),[],1)];
        % Negative reductions indicate increased signal during clamping.
        xlim([min([0;panelReduction(isfinite(panelReduction))])-5,105]);
        xlabel('Reduction by clamping (%) (higher = better)');
        ylabel('Animal'); title(qualityTitles{m});
        % plotScatterBar hides its handles, so create explicit legend markers.
        qualityLegend = gobjects(4,1);
        qualityLegend(1) = plot(NaN,NaN,'o','MarkerFaceColor',qualityPlotColor, ...
            'MarkerEdgeColor',qualityPlotColor,'MarkerSize',sqrt(qualityPooledDotSize));
        for e = 1:numel(qualityPlotEvents)
            qualityLegend(e+1) = plot(NaN,NaN,'o','MarkerFaceColor',qualityEventColors(e,:), ...
                'MarkerEdgeColor',qualityEventColors(e,:),'MarkerSize',sqrt(qualityEventDotSize));
        end
        legend(qualityLegend,[{'All events (trial-weighted)'},qualityPlotEvents],'Location','best');
    end
    sgtitle([qualitySignal,' clamp quality'],'Interpreter','none');
end


return










%% 4. Paper-style spontaneous and reward transient panels (E/F)
% Run this section independently after loadNeuroDAP. It loads its own sources
% without replacing the animals variable used by sections 1-3.
% Layout reference: https://doi.org/10.64898/2026.07.24.740623
% Peak inclusion follows the paper's Fiber photometry methods: local
% minimum-to-maximum rise >=1.5% dF/F; exclude +/-2 s around reward.
% Saved traces retain NeuroDAP preprocessing, not the paper's raw-signal
% normalization/filtering. See measureSpontaneousTransients for details.
transientFolder = osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Project clamping/Results/202604-BiPOLES');
transientSettings.spontaneousFile = fullfile(transientFolder,'animals_20260722_random-toM446.mat');
transientSettings.rewardFile = fullfile(transientFolder,'animals_20260722_random-toM446.mat');
transientSettings.animals = 'All'; % or {'SL431','SL432','SL433'}
% One selection per condition (OFF, ON). Task uses substring matching;
% signal/event use exact matching. Edit these for other animals MAT files.
% Water-aligned windows start AFTER the reward exclusion and end before the
% next cue minus 2 s, preventing overlap and excluding neighboring events.
% These are sampled inter-event windows, not whole-session frequencies.
transientSettings.spontaneousSelections = struct( ...
    'task',{'Random','Random'},'signal',{'NAc-clamp','NAc-clamp'}, ...
    'event',{'Water (unclamp)','Water (clamp)'}, ...
    'window',{[2,4.5],[2,4.5]}, ...
    'rewardAligned',{true,true},'cueToRewardSeconds',{0,0}, ...
    'confirmedEventFree',{false,false});
% Random water is delivered at the cue. For delayed-reward tasks, set
% cueToRewardSeconds to the protocol's cue-to-water delay before use.
% For event-free ON/OFF recordings, for example, use task='Final',
% signal='NAc-left', event='Clamp', windows [6,9.5] and [1,4.5],
% rewardAligned=false, confirmedEventFree=true in both selections.
% Confirm clamp state and absence of cues/rewards for any custom windows.
transientSettings.rewardSignal = 'NAc-clamp';
transientSettings.rewardEvents = {'Water (unclamp)','Water (clamp)'};
transientSettings.rewardTime = -1:0.02:2;
transientSettings.rewardBaseline = [-0.5,0];
transientSettings.rewardPeakWindow = [0,1.5];
transientSettings.scale = 100; % fractional dF/F -> percent; use 1 if already %
% Not valid for z-scored input. downsampleSignal(...,dff=true) stores fractions.
transientSettings.minimumRisePercent = 1.5;
transientSettings.rewardExclusionSeconds = 2;
transientSettings.peakDefinition = 'Local maximum minus immediately preceding local minimum';
transientSettings.preprocessing = 'Existing NeuroDAP dF/F; 5 ms moving mean rounded to samples';
transientSettings.methodSource = 'https://doi.org/10.64898/2026.07.24.740623';
% The paper's 5 ms sliding average rounds to one sample at 50 Hz. No
% additional 100 ms smoothing, prominence or refractory-period criterion.
% Same palette as Setup; defaults also allow this section to run on its own.
if ~exist('clampColor','var'); clampColor = [.232 .76 .58]; end
if ~exist('unclampColor','var'); unclampColor = [165,209,178]./255; end
transientSettings.colors = [unclampColor; clampColor]; % OFF, ON
transientSettings.save = false; % opt in to MAT + PDF export below

spontaneousSource = load(transientSettings.spontaneousFile,'animals');
rewardSource = load(transientSettings.rewardFile,'animals');
transientPanels = clampTransientPanels(spontaneousSource.animals, ...
    rewardSource.animals,transientSettings);
disp(transientPanels.spontaneous);
disp(transientPanels.reward);
if transientSettings.save
    save(fullfile(transientFolder,'clampTransientPanels.mat'),'transientPanels','transientSettings');
    exportgraphics(transientPanels.figure,fullfile(transientFolder,'clampTransientPanels.pdf'), ...
        'ContentType','vector');
end

function out = clampTransientPanels(spontaneous, reward, cfg)
% Session metrics are averaged within animal; animals are plotted as dots.
[out.spontaneous,out.spontaneousSessions,out.spontaneousWindows] = ...
    measureSpontaneousTransients(spontaneous,cfg);
spNames = out.spontaneous.Animal;
spPeak = [out.spontaneous.OffPeak,out.spontaneous.OnPeak];
spRate = [out.spontaneous.OffHz,out.spontaneous.OnHz];
rwRows = contains({reward.task},'Random','IgnoreCase',true) & ...
    strcmpi({reward.name},cfg.rewardSignal);
rwNames = transientAnimalNames(reward,rwRows,cfg.animals);

rwPeak = nan(numel(rwNames),2); rwCount = zeros(numel(rwNames),2);
rwTrace = nan(numel(rwNames),numel(cfg.rewardTime),2);
for a = 1:numel(rwNames)
    for c = 1:2
        rows = find(rwRows & strcmpi({reward.animal},rwNames(a)) & ...
            strcmpi({reward.event},cfg.rewardEvents{c}));
        traces = [];
        for r = reshape(rows,1,[])
            [x,t] = transientWindow(reward(r),[cfg.rewardTime(1),cfg.rewardTime(end)]);
            if isempty(x); continue; end
            base = t>=cfg.rewardBaseline(1) & t<cfg.rewardBaseline(2);
            if ~any(base); continue; end
            x = cfg.scale*(x-mean(x(:,base),2));
            % Common grid permits animal averaging across sampling rates.
            x = interp1(t,x',cfg.rewardTime,'linear')';
            traces = [traces; x(all(isfinite(x),2),:)]; %#ok<AGROW>
        end
        rwCount(a,c) = size(traces,1);
        if isempty(traces); continue; end
        rwTrace(a,:,c) = mean(traces,1);
        peakIdx = cfg.rewardTime>=cfg.rewardPeakWindow(1) & ...
            cfg.rewardTime<=cfg.rewardPeakWindow(2);
        % Peak of each animal's trial-averaged response, not absolute peak.
        rwPeak(a,c) = max(rwTrace(a,peakIdx,c));
    end
end
out.reward = table(rwNames,rwPeak(:,1),rwPeak(:,2),rwCount(:,1),rwCount(:,2), ...
    'VariableNames',{'Animal','OffPeak','OnPeak','OffTrials','OnTrials'});
out.rewardTime = cfg.rewardTime; out.rewardAnimalTraces = rwTrace;
out.settings = cfg;
out.figure = figure('Color','w','Position',[100,100,1400,380]);
tiledlayout(1,4,'TileSpacing','compact','Padding','compact');
nexttile; out.p(1) = transientBar(spPeak,cfg.colors);
ylabel('Peak \DeltaF/F (%)'); title('Spontaneous transients');
nexttile; out.p(2) = transientBar(spRate,cfg.colors);
ylabel('Frequency (Hz)'); title('Spontaneous transients');
nexttile; hold on;
paired = all(isfinite(rwPeak),2); n = sum(paired);
h = gobjects(1,2);
for c = 1:2
    x = rwTrace(paired,:,c); mu = mean(x,1);
    if n>1
        se = std(x,0,1)/sqrt(n);
        fill([cfg.rewardTime,fliplr(cfg.rewardTime)],[mu-se,fliplr(mu+se)], ...
            cfg.colors(c,:),'FaceAlpha',0.2,'EdgeColor','none','HandleVisibility','off');
    end
    h(c) = plot(cfg.rewardTime,mu,'Color',cfg.colors(c,:),'LineWidth',2);
end
xline(0,':','HandleVisibility','off'); xlim(cfg.rewardTime([1,end]));
xlabel('Time from water (s)'); ylabel('\DeltaF/F (%)');
title(sprintf('Reward transients (n = %d)',n));
legend(h,{'Unclamp','Clamp'},'Location','northeast','Box','off'); box off;
nexttile; out.p(3) = transientBar(rwPeak,cfg.colors);
ylabel('Peak \DeltaF/F (%)'); title('Reward transients');
out.pairedAnimals = {spNames(all(isfinite(spPeak),2)), ...
    spNames(all(isfinite(spRate),2)),rwNames(paired)};
set(findall(out.figure,'Type','axes'),'FontSize',12,'LineWidth',1);
end

function names = transientAnimalNames(data,rows,selection)
names = unique(string({data(rows).animal}))';
if ~isequal(selection,'All')
    names = names(ismember(lower(names),lower(string(selection))));
end
end

function [x,t,trialIDs] = transientWindow(row,window,halfOpen)
% Half-open windows have duration N/Fs; require complete finite windows.
if nargin<3; halfOpen=false; end
t = double(row.timestamp(:)'); x = []; trialIDs = [];
if numel(t)<3 || size(row.data,2)~=numel(t) || ...
        any(~isfinite(t)) || any(diff(t)<=0)
    warning('ClampQuality:Timestamp','Invalid timestamps for %s/%s.',row.animal,row.event);
    return;
end
dt = median(diff(t));
if max(abs(diff(t)-dt))>dt*0.01 || t(1)>window(1)+dt/2 || t(end)<window(2)-dt/2
    warning('ClampQuality:Window','Missing/irregular window for %s/%s.',row.animal,row.event);
    return;
end
% Reward interpolation includes the endpoint; spontaneous frequency uses
% half-open windows regardless of the configured event name.
idx = t>=window(1)-dt*1e-6 & t<=window(2)+dt*1e-6;
if halfOpen; idx = idx & t<window(2)-dt*1e-6; end
if sum(idx)<3; return; end
keep = true(size(row.data,1),1);
if isfield(row,'trialInfo') && isfield(row.trialInfo,'trialTable') && ...
        ismember('performing',row.trialInfo.trialTable.Properties.VariableNames)
    keep = row.trialInfo.trialTable.performing==1;
    assert(numel(keep)==size(row.data,1),'Trial table/data length mismatch.');
else
    warning('ClampQuality:Performing','No performing flag for %s/%s; using all trials.',row.animal,row.event);
end
x = double(row.data(keep,idx)); t = t(idx);
trialIDs = find(keep);
valid = all(isfinite(x),2);
x = x(valid,:); trialIDs = trialIDs(valid);
end

function p = transientBar(values,colors)
% Matched animals for each panel; missing counterparts are not plotted.
values = values(all(isfinite(values),2),:); n = size(values,1);
hold on; p = NaN;
for c = 1:2
    if n==0; continue; end
    plotScatterBar(c,values(:,c),style='bar',color=colors(c,:), ...
        plotScatter=true,connectPairs=false,dotSize=70);
end
% Paired, two-sided signed-rank test; no trial-level pseudoreplication.
if n>=2
    if all(values(:,1)==values(:,2)); p=1;
    else; p=signrank(values(:,1),values(:,2)); end
end
xticks([1,2]); xticklabels({'Unclamp','Clamp'}); xlim([0.5,2.5]); box off;
if n==0; label='No paired data';
elseif isnan(p); label=sprintf('n = %d',n);
else; label=sprintf('n = %d, p = %.3g',n,p); end
xlabel(label);
end
