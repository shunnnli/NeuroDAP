function tests = testClampCalibration
tests = functiontests(localfunctions);
end

function testPowerGroupsAndResponseNormalization(testCase)
[raw,blue,red,Fs] = makeSweeps;
result = analyzeClampCalibration(raw,blue,red,Fs,emaTau=0);
verifyEqual(testCase,[result.groups.power_pct],[20 60 25 75]);
verifyEqual(testCase,[result.groups.nTrials],[2 1 2 1]);
verifyEqual(testCase,[result.groups.medianResponse_dff],[0.1 0.3 -0.1 -0.3],'AbsTol',1e-12);
verifyEqual(testCase,result.events.onset_idx,(5:10:55)'*Fs+1);
verifyEqual(testCase,result.events.offset_idx,(10:10:60)'*Fs+1);
verifyEqual(testCase,result.events.duration_sec,5*ones(6,1));
verifyEqual(testCase,result.analysisFs,200);
verifyEqual(testCase,result.groups(1).time_sec([1 end]),[-1 7]);
verifyEqual(testCase,result.fits(1).slope,0.005,'AbsTol',1e-12);
verifyEqual(testCase,result.fits(2).slope,-0.004,'AbsTol',1e-12);
verifyEqual(testCase,[result.fits.r2],[1 1],'AbsTol',1e-12);
end

function testBoundaryTrialsAndMissingPostWindow(testCase)
Fs = 200;
raw = 500*ones(6000,1)*5/1023;
red = zeros(size(raw)); blue = red;
red(1:400) = laserVoltage(20,[25 500]); % recording starts during stimulation
red(1901:2300) = laserVoltage(20,[25 500]);
red(4901:5900) = laserVoltage(60,[25 500]); % complete pulse, short post-window
blue(5801:end) = laserVoltage(30,[800 1600]); % recording ends during stimulation
result = analyzeClampCalibration(raw,blue,red,Fs,emaTau=0);
verifyEqual(testCase,result.events.validTrial,[false;true;true;false]);
verifyEqual(testCase,[result.groups.nTrials],[1 1 0]);
verifyTrue(testCase,all(isnan(result.groups(1).dff(1,:))));
verifyTrue(testCase,isnan(result.groups(2).dff(1,end)));
verifyEqual(testCase,result.groups(2).medianResponse_dff,0,'AbsTol',1e-12);
verifyEqual(testCase,result.events.offset_idx(end),numel(raw)+1);
end

function testShortBaselineAndZeroF0(testCase)
Fs = 200;
red = zeros(6000,1); blue = red; raw = red;
red(101:500) = laserVoltage(20,[25 500]); % less than 1 s baseline
red(2001:2400) = laserVoltage(60,[25 500]); % full baseline, but F0 = 0
result = analyzeClampCalibration(raw,blue,red,Fs);
verifyFalse(testCase,any(result.events.validTrial));
verifyTrue(testCase,all(isnan([result.groups.medianResponse_dff])));
end

function testPulseFilteringAndGapMerging(testCase)
Fs = 200;
red = zeros(10000,1); blue = red;
red(1001:2000) = laserVoltage(40,[25 500]);
red(1401:1405) = 0; % a brief dropout should not split the pulse
red(3001:3020) = laserVoltage(40,[25 500]); % too short
red(5001:7000) = laserVoltage(40,[25 500]); % long block, not a sweep
result = analyzeClampCalibration(500*ones(size(red))*5/1023,blue,red,Fs);
verifyEqual(testCase,height(result.events),1);
verifyEqual(testCase,result.events.onset_idx,1001);
verifyEqual(testCase,result.events.offset_idx,2001);
verifyEqual(testCase,result.groups.power_pct,40);
end

function testNoPulsesAndInvalidInputs(testCase)
raw = ones(1000,1); laser = zeros(size(raw));
verifyWarning(testCase,@() analyzeClampCalibration(raw,laser,laser,200), ...
    'analyzeClampCalibration:NoPulses');
verifyError(testCase,@() analyzeClampCalibration(raw,laser(2:end),laser,200), ...
    'analyzeClampCalibration:LengthMismatch');
verifyError(testCase,@() analyzeClampCalibration(raw,laser,laser,200,redClampRange=[500 25]), ...
    'analyzeClampCalibration:InvalidRange');
end

function testPlotsSeparateChannelsAndPowers(testCase)
[raw,blue,red,Fs] = makeSweeps;
result = analyzeClampCalibration(raw,blue,red,Fs);
fig = plotClampCalibration(result);
cleanup = onCleanup(@() close(fig)); %#ok<NASGU>
traces = findobj(fig,'Type','line','-regexp','DisplayName','% \(n=');
verifyEqual(testCase,numel(traces),4);
verifyEqual(testCase,sort(string({traces.DisplayName})), ...
    sort(["20% (n=2)","60% (n=1)","25% (n=2)","75% (n=1)"]));
end

function testCalibrationSessionRoutingAndSavedResults(testCase)
[photometry_raw,blueClamp,redClamp,Fs] = makeSweeps;
folder = fullfile(tempname,'20260920-Test-CaLiBrAtIoN');
mkdir(folder);
cleanup = onCleanup(@() rmdir(fileparts(folder),'s')); %#ok<NASGU>
[~,sessionName] = fileparts(folder);
% Deliberately omit NI photometry, timeseries, behavior, and clampTarget.
[labjack,params] = makeLabjack(photometry_raw,Fs);
save(fullfile(folder,['data_',sessionName,'.mat']),'labjack','blueClamp','redClamp');
params.session.task = 'random';
params.analyze = struct('legacySetting',42); % no behavioral analysis options
save(fullfile(folder,['sync_',sessionName,'.mat']),'params');
analyzeSessions_clamp(folder,task='reward pairing',plotPhotometry=false,redo=false);
saved = load(fullfile(folder,['analysis_',sessionName,'.mat']));
sync = load(fullfile(folder,['sync_',sessionName,'.mat']));
behavior = load(fullfile(folder,['behavior_',sessionName,'.mat']));
verifyEqual(testCase,sync.params.session.task,'calibration');
verifyEqual(testCase,sync.params.analyze.legacySetting,42);
verifyEqual(testCase,[saved.calibrationAnalysis.signals(1).groups.power_pct],[20 60 25 75]);
verifyEqual(testCase,saved.calibrationAnalysis.signalNames,["green" "red"]);
verifyEqual(testCase,size(behavior.calibrationEvents.validTrial),[6 2]);
verifyEqual(testCase,height(behavior.calibrationEvents),6);
verifyFalse(testCase,isfield(behavior,'trials'));
% Changed options must invalidate the cached analysis even with redo=false.
analyzeSessions_clamp(folder,plotPhotometry=false,redo=false,analyzeTraces=false, ...
    calibration=struct('preTime',0.5),redClampRange=[25 600]);
saved = load(fullfile(folder,['analysis_',sessionName,'.mat']));
verifyEqual(testCase,saved.calibrationAnalysis.options.preTime,0.5);
verifyEqual(testCase,saved.calibrationAnalysis.options.redClampRange,[25 600]);
end

function testStandaloneCalibrationPreservesFilesAndCache(testCase)
[photometry_raw,blueClamp,redClamp,Fs] = makeSweeps;
folder = fullfile(tempname,'20260921-Test-Cal');
mkdir(folder);
cleanup = onCleanup(@() rmdir(fileparts(folder),'s')); %#ok<NASGU>
[~,sessionName] = fileparts(folder);
dataFile = fullfile(folder,['data_',sessionName,'.mat']);
[labjack,params] = makeLabjack(photometry_raw,Fs);
save(dataFile,'labjack','blueClamp','redClamp');
preserved = 123;
save(fullfile(folder,['sync_',sessionName,'.mat']),'params','preserved');
save(fullfile(folder,['behavior_',sessionName,'.mat']),'preserved');
save(fullfile(folder,['analysis_',sessionName,'.mat']),'preserved');
result = analyzeSessions_clampCalibration(folder,plotPhotometry=false, ...
    calibration=struct('emaTau',0));
verifyEqual(testCase,[result.signals(1).groups.medianResponse_dff], ...
    [0.1 0.3 -0.1 -0.3],'AbsTol',1e-12);
verifyEqual(testCase,[result.signals(2).groups.medianResponse_dff], ...
    [-0.1 -0.3 0.1 0.3],'AbsTol',1e-12);
verifyFalse(testCase,result.signals(1).options.convertToADC);
verifyFalse(testCase,isfield(result.signals(1).groups,'raw_adc'));
for prefix = ["sync_","behavior_","analysis_"]
    saved = load(fullfile(folder,prefix+sessionName+".mat"));
    verifyEqual(testCase,saved.preserved,123);
end
% A valid cache should not reload raw inputs; a settings change should.
delete(dataFile);
cached = analyzeSessions_clampCalibration(folder,plotPhotometry=false, ...
    calibration=struct('emaTau',0),redo=false,analyzeTraces=false);
verifyEqual(testCase,cached,result);
% A saved calibration task also routes through the main entry point when
% the directory itself does not contain "calibration".
analyzeSessions_clamp(folder,plotPhotometry=false,calibration=struct('emaTau',0), ...
    redo=false,analyzeTraces=false);
end

function testCalibrationFolderWithOutputAlias(testCase)
[photometry_raw,blueClamp,redClamp,Fs] = makeSweeps;
folder = fullfile(tempname,'Calibration');
mkdir(folder);
cleanup = onCleanup(@() rmdir(fileparts(folder),'s')); %#ok<NASGU>
[labjack,params] = makeLabjack(photometry_raw,Fs);
save(fullfile(folder,'data_recording.mat'),'labjack','blueClamp','redClamp');
save(fullfile(folder,'sync_recording.mat'),'params');
% Name-based routing takes precedence even when outputName has no task label.
analyzeSessions_clamp(folder,outputName='recording',task='random',plotPhotometry=false);
saved = load(fullfile(folder,'analysis_recording.mat'));
verifyEqual(testCase,height(saved.calibrationAnalysis.events),6);
verifyEqual(testCase,string(saved.sessionName),"Calibration");
end

function testLabjackClockOffsetDriftAndDifferentSamplingRates(testCase)
[raw,blue,red,Fs] = makeSweeps;
timeNI = -3+(0:numel(raw)-1)/Fs;
nativeTimeLJ = (0:24000)/400;
timeLJ = 0.5+nativeTimeLJ*1.025+0.1*sin(nativeTimeLJ/8);
fluorescence = interp1(timeNI,raw,timeLJ,'previous',raw(1));
result = analyzeClampCalibration(fluorescence,blue,red,Fs, ...
    photometryTime=timeLJ,laserTime=timeNI,convertToADC=false,emaTau=0);
verifyEqual(testCase,result.events.onset_idx,(5:10:55)'*Fs+1);
verifyEqual(testCase,result.events.onset_sync_sec,(2:10:52)');
verifyEqual(testCase,[result.groups.medianResponse_dff], ...
    [0.1 0.3 -0.1 -0.3],'AbsTol',1e-12);
% Inspect early response timing too: a start-offset-only alignment would
% drift off by over a second by the last pulse.
for g = 1:numel(result.groups)
    group = result.groups(g);
    atOnset = find(group.time_sec >= 0.1,1);
    verifyEqual(testCase,group.mean_dff(atOnset),group.medianResponse_dff,'AbsTol',1e-12);
end
verifyError(testCase,@() analyzeClampCalibration(fluorescence,blue,red,Fs, ...
    photometryTime=timeLJ,convertToADC=false),'analyzeClampCalibration:MissingSync');
verifyError(testCase,@() analyzeClampCalibration(fluorescence,blue,red,Fs, ...
    photometryTime=fliplr(timeLJ),laserTime=timeNI,convertToADC=false), ...
    'analyzeClampCalibration:InvalidSync');
end

function testLabjackCoverageAndFiltering(testCase)
[raw,blue,red,Fs] = makeSweeps;
timeNI = (0:numel(raw)-1)/Fs;
% LabJack starts halfway through the first baseline and ends during the last
% pulse. NaN padding must not propagate through the EMA into valid trials.
selected = 4501:58001;
result = analyzeClampCalibration(raw(selected),blue,red,Fs, ...
    photometryTime=timeNI(selected),laserTime=timeNI,convertToADC=false);
verifyEqual(testCase,result.events.validTrial,[false;true;true;true;true;false]);
verifyEqual(testCase,[result.groups.nTrials],[1 1 2 0]);
verifyEqual(testCase,result.groups(2).medianResponse_dff,0.3,'AbsTol',1e-10);
end

function testSessionDemodulationAndLegacyCache(testCase)
Fs = 1000;
nativeFs = 2000;
blueClamp = zeros(15*Fs,1); redClamp = blueClamp;
redClamp(5001:10000) = laserVoltage(40,[25 500]);
t = (0:15*nativeFs-1)/nativeFs;
envelope = 1+0.25*(t >= 5 & t < 10);
labjack = struct('raw',envelope.*cos(2*pi*200*t),'name',{{'modulated'}}, ...
    'samplerate',nativeFs,'mod',true,'modFreq',200);
params.sync = struct('behaviorFs',Fs,'labjackFs',nativeFs, ...
    'timeNI',(0:numel(redClamp)-1)/Fs,'timePhotometry',t);
folder = fullfile(tempname,'20260921-Test-Calibration');
mkdir(folder);
cleanup = onCleanup(@() rmdir(fileparts(folder),'s')); %#ok<NASGU>
[~,sessionName] = fileparts(folder);
save(fullfile(folder,['data_',sessionName,'.mat']),'labjack','blueClamp','redClamp');
syncFile = fullfile(folder,['sync_',sessionName,'.mat']);
save(syncFile,'params');
% Old NI-only results cannot be reused, even if all requested options match.
calibrationAnalysis = struct('source','photometry_raw (NI)', ...
    'requestedOptions',struct('blueClampRange',[800 1600],'redClampRange',[25 500]));
save(fullfile(folder,['analysis_',sessionName,'.mat']),'calibrationAnalysis');
result = analyzeSessions_clampCalibration(folder,plotPhotometry=false,redo=false,analyzeTraces=false);
verifyEqual(testCase,result.source,'LabJack');
verifyTrue(testCase,result.signals.modulated);
verifyEqual(testCase,result.signals.groups.medianResponse_dff,0.25,'AbsTol',1e-6);
% Missing sync data must produce an actionable error, never assume the
% devices started together or fall back to NI photometry.
params.sync = rmfield(params.sync,'timePhotometry');
save(syncFile,'params');
verifyError(testCase,@() analyzeSessions_clampCalibration(folder,plotPhotometry=false), ...
    'analyzeSessions_clampCalibration:MissingSync');
end

function [raw,blue,red,Fs] = makeSweeps
Fs = 1000;
raw = 500*ones(65*Fs,1);
blue = zeros(size(raw)); red = blue;
for trial = 1:6
    idx = (5+(trial-1)*10)*Fs+(1:5*Fs);
    powers = [20 20.4 60 25 25.5 75];
    responses = [50 50 150 -50 -50 -150];
    if trial <= 3
        red(idx) = laserVoltage(powers(trial),[25 500]);
    else
        blue(idx) = laserVoltage(powers(trial),[800 1600]);
    end
    raw(idx) = raw(idx)+responses(trial);
end
raw = raw*5/1023;
end

function voltage = laserVoltage(power,range)
voltage = (range(1)+(range(2)-range(1))*power/100)*5/4095;
end

function [labjack,params] = makeLabjack(raw,Fs)
labjack.raw = [raw(:)';2*raw(1)-raw(:)'];
labjack.name = {'green','red'};
labjack.mod = [false false];
labjack.samplerate = Fs;
params.sync.behaviorFs = Fs;
params.sync.labjackFs = Fs;
params.sync.timeNI = (0:numel(raw)-1)/Fs;
params.sync.timePhotometry = params.sync.timeNI;
end
