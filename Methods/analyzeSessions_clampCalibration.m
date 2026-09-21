function calibrationAnalysis = analyzeSessions_clampCalibration(sessionpath,options)
%ANALYZESESSIONS_CLAMPCALIBRATION Analyze and plot calibration power sweeps.
% Called automatically by analyzeSessions_clamp for calibration sessions, or
% directly with a session directory. Requires data_*.mat (labjack,
% blueClamp, redClamp) and the NI/LabJack sync timestamps in sync_*.mat. Behavioral
% trials, timeseries_*.mat, and clampTarget are not required.
%
% The response analysis follows the step-sweep section of PID-tuning.ipynb:
% 200 Hz sampling, 20 ms EMA, per-trial pre-onset F0, and mean +/- SEM traces
% by power, with late-on median responses and linear power-response fits.
% All recorded labjack.raw channels are analyzed separately. Modulated
% channels use demodulateSignal's non-detrended fluorescence. No Arduino ADC
% quantization or rolling-z normalization is applied to LabJack signals.
% Results are stored in calibrationAnalysis.signals, one entry per channel.
% Numerical settings are passed in the calibration struct; see
% analyzeClampCalibration. No YAML/TOML or PID configuration is written.
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
    options.redo (1,1) logical = true
    options.analyzeTraces (1,1) logical = true
    options.plotPhotometry (1,1) logical = true
end

sessionpath = strip(sessionpath,'right',filesep);
[projectPath,folderName] = fileparts(sessionpath);
if ~isfield(options,'outputName'); options.outputName = folderName; end
sessionName = options.outputName;
if ~contains(sessionName,{'-','_'})
    if contains(folderName,'calibration',IgnoreCase=true)
        sessionName = folderName;
    else
        % Match the parent-session naming used for split recordings.
        [~,sessionName] = fileparts(projectPath);
    end
end
dataFile = fullfile(sessionpath,"data_"+options.outputName+".mat");
syncFile = fullfile(sessionpath,"sync_"+options.outputName+".mat");
analysisFile = fullfile(sessionpath,"analysis_"+options.outputName+".mat");
behaviorFile = fullfile(sessionpath,"behavior_"+options.outputName+".mat");
syncData = load(syncFile,'params');
if ~isfield(syncData,'params') || ~isfield(syncData.params,'sync') || ...
        ~isfield(syncData.params.sync,'behaviorFs')
    error('analyzeSessions_clampCalibration:MissingSamplingRate', ...
        'Calibration requires params.sync.behaviorFs in %s.',syncFile);
end
params = syncData.params;

calibrationOptions = options.calibration;
calibrationOptions.blueClampRange = options.blueClampRange;
calibrationOptions.redClampRange = options.redClampRange;
cached = struct();
if isfile(analysisFile)
    variables = whos('-file',analysisFile);
    if any(strcmp({variables.name},'calibrationAnalysis'))
        cached = load(analysisFile,'calibrationAnalysis');
    end
end
% Settings changes invalidate the cache even when redo and analyzeTraces are
% false. Load only the calibration inputs, without other session variables.
if options.redo || options.analyzeTraces || ~isfield(cached,'calibrationAnalysis') || ...
        ~isfield(cached.calibrationAnalysis,'schemaVersion') || ...
        cached.calibrationAnalysis.schemaVersion ~= 2 || ...
        ~isfield(cached.calibrationAnalysis,'requestedOptions') || ...
        ~isequaln(cached.calibrationAnalysis.requestedOptions,calibrationOptions)
    variables = whos('-file',dataFile);
    required = {'labjack','blueClamp','redClamp'};
    missing = required(~ismember(required,{variables.name}));
    if ~isempty(missing)
        error('analyzeSessions_clampCalibration:MissingCalibrationData', ...
            'Calibration requires %s in %s.',strjoin(missing,', '),dataFile);
    end
    data = load(dataFile,required{:});
    calibrationAnalysis = analyzeLabjackChannels(data,params,calibrationOptions);
    calibrationAnalysis.requestedOptions = calibrationOptions;
    if isfile(analysisFile)
        save(analysisFile,'calibrationAnalysis','-append');
    else
        save(analysisFile,'sessionName','calibrationAnalysis','-v7.3');
    end
else
    calibrationAnalysis = cached.calibrationAnalysis;
end

calibrationEvents = calibrationAnalysis.events;
if isfile(behaviorFile)
    save(behaviorFile,'calibrationEvents','-append');
else
    save(behaviorFile,'calibrationEvents','-v7.3');
end

if ~isfield(params,'session'); params.session = struct(); end
if ~isfield(params.session,'name'); params.session.name = options.outputName; end
parts = strsplit(sessionName,{'-','_'});
if numel(parts) >= 2
    if ~isfield(params.session,'date'); params.session.date = parts{1}; end
    if ~isfield(params.session,'animal'); params.session.animal = parts{2}; end
end
if ~isfield(params.session,'projectPath'); params.session.projectPath = projectPath; end
if ~isfield(params.session,'baselineSystem'); params.session.baselineSystem = 'NI'; end
params.session.task = 'calibration';
% Record the current settings without requiring saved behavioral options.
fields = fieldnames(options);
for i = 1:numel(fields)
    params.analyze.(fields{i}) = options.(fields{i});
end
params.analyze.task = 'calibration';
save(syncFile,'params','-append');

if options.plotPhotometry
    for channel = 1:numel(calibrationAnalysis.signals)
        signal = calibrationAnalysis.signals(channel);
        fig = plotClampCalibration(signal);
        label = regexprep(char(signal.name),'[^a-zA-Z0-9_-]','_');
        figureName = sprintf('Calibration_photometry_power_response_%02d_%s',channel,label);
        saveFigures(fig,figureName,sessionpath,savePNG=true,savePDF=true);
    end
end
disp('Finished: LabJack calibration responses grouped by channel and laser power and saved');
end

function calibration = analyzeLabjackChannels(data,params,options)
lj = data.labjack;
required = {'raw','name','mod','samplerate'};
if ~isstruct(lj) || ~isscalar(lj) || ~all(isfield(lj,required)) || ...
        isempty(lj.raw) || ~ismatrix(lj.raw) || ...
        numel(string(lj.name)) ~= size(lj.raw,1) || numel(lj.mod) ~= size(lj.raw,1)
    error('analyzeSessions_clampCalibration:InvalidLabjack', ...
        'labjack must contain raw (channels by samples), matching name/mod entries, and samplerate. Rerun loadSessions.');
end
validateattributes(lj.samplerate,{'numeric'},{'scalar','finite','positive'});
if ~all(isfield(params.sync,{'timeNI','timePhotometry'}))
    error('analyzeSessions_clampCalibration:MissingSync', ...
        'LabJack calibration requires params.sync.timeNI and timePhotometry. Rerun loadSessions to synchronize the recordings.');
end
timeLJ = params.sync.timePhotometry;
if numel(timeLJ) ~= size(lj.raw,2) || numel(timeLJ) < 2 || ...
        ~isvector(timeLJ) || any(~isfinite(timeLJ(:))) || any(diff(timeLJ(:)) <= 0)
    error('analyzeSessions_clampCalibration:InvalidSync', ...
        'timePhotometry must contain one finite, increasing timestamp per raw LabJack sample.');
end
if any(lj.mod) && (~isfield(lj,'modFreq') || numel(lj.modFreq) ~= size(lj.raw,1))
    error('analyzeSessions_clampCalibration:MissingModulationFrequency', ...
        'Modulated LabJack channels require a matching labjack.modFreq entry.');
end

options.convertToADC = false;
options.laserTime = params.sync.timeNI;
names = string(lj.name);
signals = struct([]);
for channel = 1:size(lj.raw,1)
    fluorescence = double(lj.raw(channel,:));
    options.photometryTime = timeLJ;
    modFreq = NaN;
    units = 'V';
    if lj.mod(channel)
        modFreq = lj.modFreq(channel);
        validateattributes(modFreq,{'numeric'},{'scalar','finite','positive','<',lj.samplerate/2});
        targetFs = 200;
        if isfield(options,'targetFs'); targetFs = options.targetFs; end
        processed = demodulateSignal(fluorescence,originalFs=lj.samplerate, ...
            targetFs=min(targetFs,lj.samplerate),modFreq=modFreq,resample=false);
        fluorescence = processed.demodData_nodetrend;
        % Spectrogram timestamps describe window centers, not the first raw
        % sample. Map those centers through the same synchronized LJ clock.
        options.photometryTime = interp1(1:numel(timeLJ),timeLJ, ...
            1+processed.demodTimes*lj.samplerate,'linear');
        units = 'demodulated amplitude';
    end
    calibrationArgs = namedargs2cell(options);
    result = analyzeClampCalibration(fluorescence,data.blueClamp,data.redClamp, ...
        params.sync.behaviorFs,calibrationArgs{:});
    result.name = names(channel);
    result.system = 'LJ';
    result.source = 'labjack.raw';
    result.signalUnits = units;
    result.nativeFs = lj.samplerate;
    result.modulated = logical(lj.mod(channel));
    result.modFreq = modFreq;
    result.alignment = struct('method','shared sync timestamps', ...
        'startOffset_sec',timeLJ(1)-params.sync.timeNI(1));
    if channel == 1
        signals = result;
    else
        signals(channel) = result;
    end
end
events = signals(1).events;
% Detection is shared; baseline/recording coverage can differ by channel.
events.validTrial = false(height(events),numel(signals));
for channel = 1:numel(signals)
    events.validTrial(:,channel) = signals(channel).events.validTrial;
end
calibration = struct('schemaVersion',2,'source','LabJack','signalNames',names, ...
    'originalFs',params.sync.behaviorFs,'analysisFs',signals(1).analysisFs, ...
    'options',signals(1).options,'events',events,'signals',signals);
end
