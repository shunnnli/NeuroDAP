function calibrationAnalysis = analyzeSessions_clampCalibration(sessionpath,options)
%ANALYZESESSIONS_CLAMPCALIBRATION Analyze and plot calibration power sweeps.
% Called automatically by analyzeSessions_clamp for calibration sessions, or
% directly with a session directory. Requires data_*.mat (photometry_raw,
% blueClamp, redClamp) and sync_*.mat (params.sync.behaviorFs). Behavioral
% trials, timeseries_*.mat, and clampTarget are not required.
%
% The response analysis follows the step-sweep section of PID-tuning.ipynb:
% 200 Hz sampling, 20 ms EMA, per-trial pre-onset F0, and mean +/- SEM traces
% by power, with late-on median responses and linear power-response fits.
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
        ~isfield(cached.calibrationAnalysis,'requestedOptions') || ...
        ~isequaln(cached.calibrationAnalysis.requestedOptions,calibrationOptions)
    variables = whos('-file',dataFile);
    required = {'photometry_raw','blueClamp','redClamp'};
    missing = required(~ismember(required,{variables.name}));
    if ~isempty(missing)
        error('analyzeSessions_clampCalibration:MissingCalibrationData', ...
            'Calibration requires %s in %s.',strjoin(missing,', '),dataFile);
    end
    data = load(dataFile,required{:});
    calibrationArgs = namedargs2cell(calibrationOptions);
    calibrationAnalysis = analyzeClampCalibration(data.photometry_raw, ...
        data.blueClamp,data.redClamp,params.sync.behaviorFs,calibrationArgs{:});
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
    fig = plotClampCalibration(calibrationAnalysis);
    saveFigures(fig,'Calibration_photometry_power_response',sessionpath, ...
        savePNG=true,savePDF=true);
end
disp('Finished: calibration responses grouped by laser power and saved');
end
