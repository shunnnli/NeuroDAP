function calibrationAnalysis = analyzeSessions_clampCalibration(sessionpath,options)
%ANALYZESESSIONS_CLAMPCALIBRATION Analyze and plot clamp calibration power sweeps.
% Called automatically by analyzeSessions_clamp for calibration sessions, or
% directly with a session directory.
%
% The function runs in four steps:
%   1. Load mat files: params from sync_*.mat, the blueClamp/redClamp command
%      traces from data_*.mat, and the processed photometry from
%      timeseries_*.mat. The blueClamp/redClamp/clampTarget command channels
%      and camera channels inside timeSeries are skipped.
%   2. Detect the calibration pulses in the NI clamp commands (findClampPulses).
%   3. Group the pulses by laser power. Each group keeps onsetIdx, the starting
%      NI sample of every pulse at that power.
%   4. Feed those starting samples to plotTraces, which aligns each photometry
%      channel to the pulse onsets and plots mean +/- SEM per power, then
%      summarize the late-on response and fit it against power.
%
% Traces are used in the units loadSessions stored them in (rolling-z or dF/F);
% the pre-onset mean is subtracted for the response statistics only. No second
% normalization, filtering, or ADC conversion is applied. Results are stored in
% calibrationAnalysis.signals, one entry per channel. Numerical settings are
% passed in the calibration struct. No YAML/TOML or PID configuration is written.
%
% Example:
%   result = analyzeSessions_clampCalibration(sessionpath, ...
%       redClampRange=[25 500],blueClampRange=[800 1600], ...
%       calibration=struct('preTime',1,'postTime',2,'lateTime',2));

arguments
    sessionpath (1,1) string
    options.outputName (1,1) string
    options.blueClampRange (1,2) double {mustBeFinite} = [800,1600]
    options.redClampRange (1,2) double {mustBeFinite} = [25,500]
    options.calibration (1,1) struct = struct()
    options.plotPhotometry (1,1) logical = true
    % Accepted so analyzeSessions_clamp can forward its own options. The
    % calibration analysis always recomputes: the photometry is already
    % processed, so there is nothing expensive left to cache.
    options.redo (1,1) logical = true
    options.analyzeTraces (1,1) logical = true
end

%% 1. Load mat files

sessionpath = strip(sessionpath,'right',filesep);
[projectPath,folderName] = fileparts(sessionpath);
if ~isfield(options,'outputName'); options.outputName = folderName; end
sessionName = resolveSessionName(options.outputName,folderName,projectPath);
dataFile       = fullfile(sessionpath,"data_"+options.outputName+".mat");
timeSeriesFile = fullfile(sessionpath,"timeseries_"+options.outputName+".mat");
syncFile       = fullfile(sessionpath,"sync_"+options.outputName+".mat");
analysisFile   = fullfile(sessionpath,"analysis_"+options.outputName+".mat");
behaviorFile   = fullfile(sessionpath,"behavior_"+options.outputName+".mat");
disp(strcat('**********',options.outputName,'**********'));

settings = mergeCalibrationOptions(options);

% sync_*.mat: sampling rates and the shared clocks used for alignment
params = loadSyncParams(syncFile,options,sessionName,projectPath);
behaviorFs = double(params.sync.behaviorFs);
timeNI = double(params.sync.timeNI);

% data_*.mat: the NI clamp commands that drive the sweep
[blueClamp,redClamp] = loadClampCommands(dataFile,timeNI);

% timeseries_*.mat: the processed photometry channels to analyze
photometry = loadPhotometryChannels(timeSeriesFile,params);
disp(['Finished: loaded ',num2str(numel(photometry)),' photometry channel(s): ', ...
    char(strjoin(string({photometry.name}),', '))]);

%% 2. Detect clamp pulses in the NI commands

[pulses,detection] = findClampPulses(blueClamp,redClamp,behaviorFs, ...
    blueClampRange=settings.blueClampRange,redClampRange=settings.redClampRange, ...
    targetFs=settings.targetFs,thresholdPct=settings.thresholdPct, ...
    minDuration=settings.minDuration,maxDuration=settings.maxDuration, ...
    minGap=settings.minGap,powerTolerance=settings.powerTolerance,laserTime=timeNI);
disp(['Finished: detected ',num2str(height(pulses)),' clamp pulses']);

%% 3. Group pulses by laser power
% One entry per laser channel and power level. onsetIdx holds the starting NI
% sample of every pulse in the group; those are the event indices used in
% step 4, and onsetTime_sec is the same onset on the shared sync clock.

lasers = ["red","blue"];
powerGroups = struct('channel',{},'type',{},'power_pct',{},'nPulses',{}, ...
    'pulseRows',{},'onsetIdx',{},'onsetTime_sec',{},'duration_sec',{},'complete',{});
for laser = lasers
    laserRows = find(pulses.channel == laser);
    for level = reshape(unique(pulses.level_pct(laserRows)),1,[])
        rows = laserRows(pulses.level_pct(laserRows) == level);
        powerGroups(end+1) = struct('channel',laser,'type',pulses.type(rows(1)), ...
            'power_pct',level,'nPulses',numel(rows),'pulseRows',rows, ...
            'onsetIdx',pulses.onset_idx(rows), ...
            'onsetTime_sec',pulses.onset_sync_sec(rows), ...
            'duration_sec',pulses.duration_sec(rows), ...
            'complete',pulses.complete(rows)); %#ok<AGROW>
    end
end
disp(['Finished: grouped pulses into ',num2str(numel(powerGroups)),' power levels']);
for laser = lasers
    levels = [powerGroups(groupsOfLaser(powerGroups,laser)).power_pct];
    disp(['     ',char(laser),' powers (% of DAC range): ',num2str(levels)]);
end

%% 4. Align each photometry channel to the pulse onsets and plot

[~,~,~,~,~,~,bluePurpleRed] = loadColors;
laserColor = struct('red',bluePurpleRed(end,:),'blue',bluePurpleRed(1,:));
channelSignals = cell(1,numel(photometry));
for channel = 1:numel(photometry)
    signal = photometry(channel).data;
    finalFs = photometry(channel).finalFs;
    system = photometry(channel).system;
    name = char(photometry(channel).name);
    units = photometry(channel).units;

    if options.plotPhotometry
        fig = initializeFig(0.67,0.8);
        tiledlayout(2,2,'TileSpacing','compact');
    end

    % 4.1 Aligned traces, one line per power level
    groups = powerGroups; % same grouping, plus this channel's traces below
    for laser = lasers
        groupIdx = groupsOfLaser(powerGroups,laser);
        if options.plotPhotometry; nexttile; hold on; end
        opacity = linspace(0.3,1,numel(groupIdx));
        for k = 1:numel(groupIdx)
            g = groupIdx(k);
            group = powerGroups(g);
            % The starting samples of this power level are NI indices;
            % plotTraces maps them onto this channel's clock through
            % params.sync and draws the mean +/- SEM of the aligned traces.
            timeRange = [-settings.preTime,max(group.duration_sec)+settings.postTime];
            [traces,timestamp] = plotTraces(group.onsetIdx,timeRange,signal, ...
                addOpacity(laserColor.(laser),opacity(k)),params, ...
                signalFs=finalFs,signalSystem=system,eventSystem='ni', ...
                plot=options.plotPhotometry,rmmissing=false, ...
                label=sprintf('%g%% (n=%d)',group.power_pct,group.nPulses));
            groups = summarizeGroup(groups,g,traces,timestamp,settings);
        end
        if options.plotPhotometry && ~isempty(groupIdx)
            plotEvent('',median(vertcat(powerGroups(groupIdx).duration_sec)), ...
                color=laserColor.(laser));
            xlabel('Time from laser onset (s)'); ylabel([name,' ',units]);
            legend('Location','best');
            title(sprintf('%s laser (%s)',laser,powerGroups(groupIdx(1)).type));
        end
    end

    % 4.2 Late-on response against power
    fits = struct('channel',{},'slope',{},'intercept',{},'r2',{});
    for laser = lasers
        groupIdx = groupsOfLaser(powerGroups,laser);
        % Sessions without pulses for this laser have no summarized groups.
        powers = []; responses = [];
        if ~isempty(groupIdx)
            powers = [groups(groupIdx).power_pct];
            responses = [groups(groupIdx).medianResponse];
        end
        fits(end+1) = fitPowerResponse(laser,powers,responses); %#ok<AGROW>
        if ~options.plotPhotometry; continue; end
        nexttile; hold on
        for k = 1:numel(groupIdx)
            trials = groups(groupIdx(k)).lateOnMedian;
            scatter(repmat(powers(k),size(trials)),trials,24,laserColor.(laser), ...
                'filled','MarkerFaceAlpha',0.3,'HandleVisibility','off');
        end
        if ~isempty(groupIdx)
            plot(powers,responses,'o-','Color',laserColor.(laser),'LineWidth',2, ...
                'DisplayName','Median across pulses');
            fit = fits(end);
            if isfinite(fit.slope)
                x = [min(powers),max(powers)];
                plot(x,fit.slope*x+fit.intercept,'--','Color',laserColor.(laser), ...
                    'DisplayName',sprintf('Fit: slope %.3g, R^2 %.3f',fit.slope,fit.r2));
            end
            legend('Location','best');
            if numel(powers) > 1
                pad = 0.05*(max(powers)-min(powers));
                xlim([min(powers)-pad,max(powers)+pad]);
            end
        end
        xlabel('Laser command (% of DAC range)');
        ylabel(sprintf('Late-on median (%s)',units));
        title(sprintf('%s laser: power response',laser));
        grid on; box off
    end

    channelSignals{channel} = struct('name',photometry(channel).name,'system',system, ...
        'source','timeSeries','signalUnits',units,'finalFs',finalFs, ...
        'demodulated',photometry(channel).demodulated, ...
        'demodFreq',photometry(channel).demodFreq, ...
        'detrend',photometry(channel).detrend, ...
        'detrend_type',photometry(channel).detrend_type, ...
        'alignment',photometry(channel).alignment, ...
        'options',settings,'groups',groups,'fits',fits, ...
        'validPulse',validPulses(pulses,groups));

    if options.plotPhotometry
        sgtitle(sprintf('%s: %s, %.3g s pre-onset baseline for the late-on median', ...
            sessionName,name,settings.preTime),'Interpreter','none');
        label = regexprep(name,'[^a-zA-Z0-9_-]','_');
        figureName = sprintf('Calibration_photometry_power_response_%02d_%s',channel,label);
        saveFigures(fig,figureName,sessionpath,savePNG=false,savePDF=true);
    end
end

%% 5. Save the analysis, the events, and the settings

signals = [channelSignals{:}];
% Coverage can differ by channel, so validTrial has one column per channel.
pulses.validTrial = false(height(pulses),numel(signals));
for channel = 1:numel(signals)
    pulses.validTrial(:,channel) = signals(channel).validPulse;
end
calibrationAnalysis = struct('schemaVersion',4,'source','timeSeries', ...
    'signalNames',string({photometry.name}),'originalFs',behaviorFs, ...
    'analysisFs',detection.analysisFs,'detection',detection,'options',settings, ...
    'powerGroups',powerGroups,'events',pulses,'signals',signals);
if isfile(analysisFile)
    save(analysisFile,'calibrationAnalysis','-append');
else
    save(analysisFile,'sessionName','calibrationAnalysis','-v7.3');
end

calibrationEvents = calibrationAnalysis.events;
if isfile(behaviorFile)
    save(behaviorFile,'calibrationEvents','-append');
else
    save(behaviorFile,'calibrationEvents','-v7.3');
end

% Record the current settings without requiring saved behavioral options.
fields = fieldnames(options);
for i = 1:numel(fields)
    params.analyze.(fields{i}) = options.(fields{i});
end
params.analyze.task = 'calibration';
save(syncFile,'params','-append');

disp('Finished: timeSeries calibration responses grouped by channel and laser power and saved');
end

%% Step 1 helpers: names, settings, and file loading

function sessionName = resolveSessionName(outputName,folderName,projectPath)
sessionName = outputName;
if ~contains(sessionName,{'-','_'})
    if contains(folderName,'calibration',IgnoreCase=true)
        sessionName = folderName;
    else
        % Match the parent-session naming used for split recordings.
        [~,sessionName] = fileparts(projectPath);
    end
end
end

function settings = mergeCalibrationOptions(options)
% Detection settings are forwarded to findClampPulses; response settings are
% used here. Unknown fields are rejected so typos cannot be ignored.
settings = struct('blueClampRange',options.blueClampRange, ...
    'redClampRange',options.redClampRange,'targetFs',200,'thresholdPct',5, ...
    'minDuration',0.5,'maxDuration',5,'minGap',2,'powerTolerance',2, ...
    'preTime',1,'postTime',2,'lateTime',2,'baselineSubtract',true);
provided = fieldnames(options.calibration);
unknown = provided(~ismember(provided,fieldnames(settings)));
if ~isempty(unknown)
    error('analyzeSessions_clampCalibration:UnknownCalibrationOption', ...
        'Unsupported calibration setting(s): %s. Supported: %s.', ...
        strjoin(unknown,', '),strjoin(fieldnames(settings),', '));
end
for i = 1:numel(provided)
    settings.(provided{i}) = options.calibration.(provided{i});
end
validateattributes(settings.preTime,{'numeric'},{'scalar','positive','finite'},'','preTime');
validateattributes(settings.postTime,{'numeric'},{'scalar','nonnegative','finite'},'','postTime');
validateattributes(settings.lateTime,{'numeric'},{'scalar','positive','finite'},'','lateTime');
validateattributes(settings.baselineSubtract,{'logical','numeric'},{'scalar','binary'},'','baselineSubtract');
settings.baselineSubtract = logical(settings.baselineSubtract);
end

function params = loadSyncParams(syncFile,options,sessionName,projectPath)
% loadSessions initializes sync_*.mat with session/sessionName and only saves
% params and timeSeries when it finishes, so a stub file means the session was
% interrupted. It also returns early while any sync file exists, so the reload
% has to be forced.
syncVariables = whos('-file',syncFile);
if ~ismember('params',{syncVariables.name})
    error('analyzeSessions_clampCalibration:MissingSamplingRate', ...
        ['%s has no params, so loadSessions never finished for this session ' ...
        '(it saves params and timeSeries last). It skips sessions that already ' ...
        'have a sync file, so rerun it with reloadAll=true.'],syncFile);
end
loaded = load(syncFile,'params');
params = loaded.params;
if ~isfield(params,'sync') || ~isfield(params.sync,'behaviorFs')
    error('analyzeSessions_clampCalibration:MissingSamplingRate', ...
        'Calibration requires params.sync.behaviorFs in %s. Rerun loadSessions with reloadAll=true.',syncFile);
end
validateattributes(params.sync.behaviorFs,{'numeric'}, ...
    {'scalar','finite','positive'},'','params.sync.behaviorFs');
if ~isfield(params.sync,'timeNI')
    error('analyzeSessions_clampCalibration:MissingSync', ...
        'Calibration requires params.sync.timeNI. Rerun loadSessions to synchronize the recordings.');
end

% Session metadata is needed before the analysis: getTraces resolves the event
% system from params.session.baselineSystem, and the pulse indices are always
% NI samples because the clamp commands are NI recordings.
if ~isfield(params,'session'); params.session = struct(); end
if ~isfield(params.session,'name'); params.session.name = options.outputName; end
parts = strsplit(sessionName,{'-','_'});
if numel(parts) >= 2
    if ~isfield(params.session,'date'); params.session.date = parts{1}; end
    if ~isfield(params.session,'animal'); params.session.animal = parts{2}; end
end
if ~isfield(params.session,'projectPath'); params.session.projectPath = projectPath; end
params.session.baselineSystem = 'NI';
params.session.task = 'calibration';
end

function [blueClamp,redClamp] = loadClampCommands(dataFile,timeNI)
variables = whos('-file',dataFile);
required = {'blueClamp','redClamp'};
missing = required(~ismember(required,{variables.name}));
if ~isempty(missing)
    error('analyzeSessions_clampCalibration:MissingCalibrationData', ...
        'Calibration requires %s in %s.',strjoin(missing,', '),dataFile);
end
data = load(dataFile,required{:});
blueClamp = double(data.blueClamp);
redClamp = double(data.redClamp);
if numel(timeNI) ~= numel(blueClamp)
    error('analyzeSessions_clampCalibration:InvalidSync', ...
        'params.sync.timeNI must contain one timestamp per NI sample (%d timestamps for %d samples).', ...
        numel(timeNI),numel(blueClamp));
end
end

function photometry = loadPhotometryChannels(timeSeriesFile,params)
% Returns one entry per photometry channel in timeSeries, with everything
% step 4 needs to align and label it.
if ~isfile(timeSeriesFile)
    error('analyzeSessions_clampCalibration:MissingTimeSeries', ...
        'Calibration requires the processed photometry in %s. Rerun loadSessions.',timeSeriesFile);
end
variables = whos('-file',timeSeriesFile);
if ~any(strcmp({variables.name},'timeSeries'))
    error('analyzeSessions_clampCalibration:MissingTimeSeries', ...
        ['%s has no timeSeries struct, only %s. loadSessions saves it last, so ' ...
        'rerun loadSessions with reloadAll=true to process the photometry.'], ...
        timeSeriesFile,strjoin({variables.name},', '));
end
loaded = load(timeSeriesFile,'timeSeries');
timeSeries = loaded.timeSeries;
required = {'name','data','finalFs','system'};
if ~isstruct(timeSeries) || isempty(timeSeries) || ~all(isfield(timeSeries,required))
    error('analyzeSessions_clampCalibration:InvalidTimeSeries', ...
        'timeSeries must be a struct array with name, data, finalFs, and system fields. Rerun loadSessions.');
end

% Photometry channels only: the clamp commands and camera traces are stored in
% timeSeries as well, but they are inputs to the sweep, not responses.
excluded = ["blueclamp","redclamp","clamptarget"];
names = string({timeSeries.name});
systems = string({timeSeries.system});
selected = find(ismember(lower(systems),["ni","lj","labjack"]) & ...
    ~ismember(lower(names),excluded) & ~cellfun(@isempty,{timeSeries.data}));
if isempty(selected)
    error('analyzeSessions_clampCalibration:NoPhotometrySignals', ...
        'timeSeries contains no NI or LabJack photometry channels to analyze.');
end
if any(ismember(lower(systems(selected)),["lj","labjack"])) && ...
        ~isfield(params.sync,'timePhotometry')
    error('analyzeSessions_clampCalibration:MissingSync', ...
        'LabJack channels require params.sync.timePhotometry. Rerun loadSessions to synchronize the recordings.');
end

photometry = struct('name',{},'system',{},'data',{},'finalFs',{},'units',{}, ...
    'demodulated',{},'demodFreq',{},'detrend',{},'detrend_type',{},'alignment',{});
for i = 1:numel(selected)
    entry = timeSeries(selected(i));
    system = char(systems(selected(i)));
    validateattributes(entry.finalFs,{'numeric'},{'scalar','finite','positive'},'','timeSeries.finalFs');
    % getTraces treats the stored trace as uniformly sampled at finalFs
    % starting when its own system started.
    if any(strcmpi(system,{'lj','labjack'}))
        signalStart = params.sync.timePhotometry(1);
    else
        signalStart = params.sync.timeNI(1);
    end
    photometry(i) = struct('name',names(selected(i)),'system',system, ...
        'data',double(entry.data(:))','finalFs',double(entry.finalFs), ...
        'units',signalUnits(entry),'demodulated',fieldValue(entry,'demux',false), ...
        'demodFreq',fieldValue(entry,'demux_freq',NaN), ...
        'detrend',fieldValue(entry,'detrend',false), ...
        'detrend_type',char(string(fieldValue(entry,'detrend_type','none'))), ...
        'alignment',struct('method','plotTraces shared sync clock', ...
        'startOffset_sec',signalStart-params.sync.timeNI(1)));
end
end

%% Step 3 and 4 helpers: grouping and response statistics

function idx = groupsOfLaser(powerGroups,laser)
if isempty(powerGroups)
    idx = [];
else
    idx = find([powerGroups.channel] == laser);
end
end

function groups = summarizeGroup(groups,g,traces,timestamp,settings)
% Adds this channel's aligned traces and the late-on response to power group g.
% A pulse needs a full pre-onset window and complete stimulation for a baseline
% and a response; everything else stays NaN.
group = groups(g);
nPulses = numel(group.onsetIdx);
preWindow = timestamp < 0;
baseline = nan(nPulses,1);
response = nan(size(traces));
lateOnMedian = baseline;
validPulse = false(nPulses,1);
for pulse = 1:nPulses
    onWindow = timestamp >= 0 & timestamp < group.duration_sec(pulse);
    if ~group.complete(pulse) || ~all(isfinite(traces(pulse,preWindow))) || ...
            ~all(isfinite(traces(pulse,onWindow)))
        continue
    end
    baseline(pulse) = mean(traces(pulse,preWindow));
    if settings.baselineSubtract
        response(pulse,:) = traces(pulse,:)-baseline(pulse);
    else
        response(pulse,:) = traces(pulse,:);
    end
    lateWindow = timestamp >= max(0,group.duration_sec(pulse)-settings.lateTime) & ...
        timestamp < group.duration_sec(pulse);
    lateOnMedian(pulse) = median(response(pulse,lateWindow),'omitnan');
    validPulse(pulse) = true;
end
counts = sum(isfinite(response),1);
sem = std(response,0,1,'omitnan')./sqrt(counts);
sem(counts == 0) = NaN;

groups(g).nTrials = sum(validPulse);
groups(g).validPulse = validPulse;
groups(g).time_sec = timestamp;
groups(g).traces = traces;
groups(g).baseline = baseline;
groups(g).response = response;
groups(g).mean_signal = mean(response,1,'omitnan');
groups(g).sem_signal = sem;
groups(g).lateOnMedian = lateOnMedian;
groups(g).medianResponse = median(lateOnMedian,'omitnan');
end

function fit = fitPowerResponse(laser,powers,responses)
finite = isfinite(powers) & isfinite(responses);
powers = powers(finite); responses = responses(finite);
slope = NaN; intercept = NaN; r2 = NaN;
if numel(unique(powers)) >= 2
    p = polyfit(powers,responses,1);
    slope = p(1); intercept = p(2);
    total = sum((responses-mean(responses)).^2);
    if total > 0; r2 = 1-sum((responses-polyval(p,powers)).^2)/total; end
end
fit = struct('channel',laser,'slope',slope,'intercept',intercept,'r2',r2);
end

function valid = validPulses(pulses,groups)
% Scatters each group's per-pulse validity back onto the shared pulse table.
valid = false(height(pulses),1);
for g = 1:numel(groups)
    valid(groups(g).pulseRows) = groups(g).validPulse;
end
end

function units = signalUnits(entry)
detrendType = lower(char(string(fieldValue(entry,'detrend_type','none'))));
switch detrendType
    case 'rolling-z'
        units = 'z-score';
    case 'dff'
        units = '\DeltaF/F';
    case 'none'
        units = 'V';
    otherwise
        units = detrendType;
end
end

function value = fieldValue(entry,name,default)
value = default;
if isfield(entry,name) && ~isempty(entry.(name)); value = entry.(name); end
end
