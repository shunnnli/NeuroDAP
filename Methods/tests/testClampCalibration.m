function tests = testClampCalibration
tests = functiontests(localfunctions);
end

%% findClampPulses: detection and power grouping

function testPulseGroupsAndIndices(testCase)
[blue,red,Fs] = makeSweeps;
[events,detection] = findClampPulses(blue,red,Fs);
verifyEqual(testCase,events.channel,["red";"red";"red";"blue";"blue";"blue"]);
verifyEqual(testCase,events.level_pct,[20;20;60;25;25;75]);
verifyEqual(testCase,events.onset_idx,(5:10:55)'*Fs+1);
verifyEqual(testCase,events.offset_idx,(10:10:60)'*Fs+1);
verifyEqual(testCase,events.duration_sec,5*ones(6,1));
verifyEqual(testCase,events.onset_sec,(5:10:55)');
verifyEqual(testCase,detection.analysisFs,200);
verifyEqual(testCase,detection.stride,5);
verifyFalse(testCase,ismember('onset_sync_sec',events.Properties.VariableNames));
end

function testBoundaryPulsesAndSyncClock(testCase)
Fs = 200;
red = zeros(6000,1); blue = red;
red(1:400) = laserVoltage(20,[25 500]); % recording starts during stimulation
red(1901:2300) = laserVoltage(20,[25 500]);
red(4901:5900) = laserVoltage(60,[25 500]);
blue(5801:end) = laserVoltage(30,[800 1600]); % recording ends during stimulation
timeNI = -3+(0:numel(red)-1)/Fs;
events = findClampPulses(blue,red,Fs,laserTime=timeNI);
verifyEqual(testCase,events.complete,[false;true;true;false]);
verifyEqual(testCase,events.offset_idx(end),numel(red)+1);
verifyEqual(testCase,events.onset_sync_sec,timeNI(events.onset_idx)','AbsTol',1e-12);
end

function testPulseFilteringAndGapMerging(testCase)
Fs = 200;
red = zeros(10000,1); blue = red;
red(1001:2000) = laserVoltage(40,[25 500]);
red(1401:1405) = 0; % a brief dropout should not split the pulse
red(3001:3020) = laserVoltage(40,[25 500]); % too short
red(5001:7000) = laserVoltage(40,[25 500]); % long block, not a sweep
events = findClampPulses(blue,red,Fs);
verifyEqual(testCase,height(events),1);
verifyEqual(testCase,events.onset_idx,1001);
verifyEqual(testCase,events.offset_idx,2001);
verifyEqual(testCase,events.level_pct,40);
end

function testNoPulsesAndInvalidInputs(testCase)
laser = zeros(1000,1);
verifyWarning(testCase,@() findClampPulses(laser,laser,200),'findClampPulses:NoPulses');
verifyError(testCase,@() findClampPulses(laser(2:end),laser,200),'findClampPulses:LengthMismatch');
verifyError(testCase,@() findClampPulses(laser,laser,200,redClampRange=[500 25]), ...
    'findClampPulses:InvalidRange');
verifyError(testCase,@() findClampPulses(laser,laser,200,minDuration=2,maxDuration=1), ...
    'findClampPulses:InvalidDuration');
verifyError(testCase,@() findClampPulses(laser,laser,200,laserTime=(1:10)/200), ...
    'findClampPulses:InvalidSync');
end

%% Session analysis of timeSeries photometry

function testCalibrationSessionRoutingAndSavedResults(testCase)
folder = makeSession(testCase,'20260920-Test-CaLiBrAtIoN');
[~,sessionName] = fileparts(folder);
analyzeSessions_clamp(folder,task='reward pairing',plotPhotometry=false,redo=false);
saved = load(fullfile(folder,['analysis_',sessionName,'.mat']));
sync = load(fullfile(folder,['sync_',sessionName,'.mat']));
behavior = load(fullfile(folder,['behavior_',sessionName,'.mat']));
verifyEqual(testCase,sync.params.session.task,'calibration');
verifyEqual(testCase,sync.params.analyze.legacySetting,42);
verifyEqual(testCase,saved.calibrationAnalysis.source,'timeSeries');
verifyEqual(testCase,saved.calibrationAnalysis.signalNames,["green" "red"]);
verifyEqual(testCase,[saved.calibrationAnalysis.signals(1).groups.power_pct],[20 60 25 75]);
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

function testResponsesPerChannelAndPowerLevel(testCase)
folder = makeSession(testCase,'20260921-Test-Cal');
result = analyzeSessions_clampCalibration(folder,plotPhotometry=false);
verifyEqual(testCase,[result.signals(1).groups.medianResponse],[0.1 0.3 -0.1 -0.3],'AbsTol',1e-12);
verifyEqual(testCase,[result.signals(2).groups.medianResponse],[-0.1 -0.3 0.1 0.3],'AbsTol',1e-12);
verifyEqual(testCase,[result.signals(1).groups.nTrials],[2 1 2 1]);
verifyEqual(testCase,result.signals(1).groups(1).time_sec([1 end]),[-1 7]);
verifyEqual(testCase,result.signals(1).fits(1).slope,0.005,'AbsTol',1e-12);
verifyEqual(testCase,result.signals(1).fits(2).slope,-0.004,'AbsTol',1e-12);
verifyEqual(testCase,[result.signals(1).fits.r2],[1 1],'AbsTol',1e-12);
verifyEqual(testCase,result.signals(1).source,'timeSeries');
verifyEqual(testCase,result.signals(1).system,'LJ');
verifyEqual(testCase,result.signals(1).signalUnits,'\DeltaF/F');
verifyEqual(testCase,result.signals(1).finalFs,100);
% Baselines come from the stored trace, not a second normalization.
verifyEqual(testCase,result.signals(1).groups(1).baseline,[0;0],'AbsTol',1e-12);
verifyTrue(testCase,all(result.events.validTrial(:)));
end

function testStandaloneCalibrationPreservesOtherVariables(testCase)
folder = makeSession(testCase,'20260921-Test-Preserve');
[~,sessionName] = fileparts(folder);
preserved = 123;
for prefix = ["sync_","behavior_","analysis_"]
    file = fullfile(folder,prefix+sessionName+".mat");
    if isfile(file); save(file,'preserved','-append'); else; save(file,'preserved'); end
end
result = analyzeSessions_clampCalibration(folder,plotPhotometry=false);
for prefix = ["sync_","behavior_","analysis_"]
    saved = load(fullfile(folder,prefix+sessionName+".mat"));
    verifyEqual(testCase,saved.preserved,123);
end
saved = load(fullfile(folder,['analysis_',sessionName,'.mat']),'calibrationAnalysis');
verifyEqual(testCase,saved.calibrationAnalysis,result);
% Analysis always recomputes, so the inputs are always required.
delete(fullfile(folder,['timeseries_',sessionName,'.mat']));
verifyError(testCase,@() analyzeSessions_clampCalibration(folder,plotPhotometry=false, ...
    redo=false,analyzeTraces=false), ...
    'analyzeSessions_clampCalibration:MissingTimeSeries');
end

function testPowerGroupsKeepPulseOnsets(testCase)
folder = makeSession(testCase,'20260921-Test-Onsets');
result = analyzeSessions_clampCalibration(folder,plotPhotometry=false);
groups = result.powerGroups;
verifyEqual(testCase,[groups.channel],["red" "red" "blue" "blue"]);
verifyEqual(testCase,[groups.power_pct],[20 60 25 75]);
verifyEqual(testCase,[groups.nPulses],[2 1 2 1]);
% Each group keeps the starting NI sample of its pulses, which is what
% plotTraces receives.
Fs = result.originalFs;
verifyEqual(testCase,groups(1).onsetIdx,[5;15]*Fs+1);
verifyEqual(testCase,groups(2).onsetIdx,25*Fs+1);
verifyEqual(testCase,groups(3).onsetIdx,[35;45]*Fs+1);
verifyEqual(testCase,groups(4).onsetIdx,55*Fs+1);
verifyEqual(testCase,groups(1).onsetTime_sec,[5;15],'AbsTol',1e-12);
verifyEqual(testCase,result.signals(1).groups(1).onsetIdx,groups(1).onsetIdx);
end

function testNonPhotometryChannelsAreSkipped(testCase)
folder = makeSession(testCase,'20260921-Test-Skip',extraChannels=true);
result = analyzeSessions_clampCalibration(folder,plotPhotometry=false);
verifyEqual(testCase,result.signalNames,["green" "red" "PMT"]);
verifyEqual(testCase,{result.signals.system},{'LJ','LJ','NI'});
% The NI channel carries the same response, aligned on its own clock.
verifyEqual(testCase,[result.signals(3).groups.medianResponse],[0.1 0.3 -0.1 -0.3],'AbsTol',1e-12);
end

function testLabjackStartOffsetAlignment(testCase)
% The LabJack recording starts 2 s after NI; responses must still land at
% laser onset rather than 2 s late.
folder = makeSession(testCase,'20260921-Test-Offset',photometryStart=2);
result = analyzeSessions_clampCalibration(folder,plotPhotometry=false);
group = result.signals(1).groups(1);
atOnset = find(group.time_sec >= 0.1,1);
verifyEqual(testCase,group.mean_signal(atOnset),group.medianResponse,'AbsTol',1e-12);
verifyEqual(testCase,result.signals(1).alignment.startOffset_sec,2,'AbsTol',1e-12);
end

function testSessionWithoutPulses(testCase)
folder = makeSession(testCase,'20260921-Test-Empty');
[~,sessionName] = fileparts(folder);
dataFile = fullfile(folder,['data_',sessionName,'.mat']);
loaded = load(dataFile);
blueClamp = zeros(size(loaded.blueClamp)); redClamp = blueClamp;
save(dataFile,'blueClamp','redClamp');
result = verifyWarning(testCase, ...
    @() analyzeSessions_clampCalibration(folder,plotPhotometry=false), ...
    'findClampPulses:NoPulses');
verifyEmpty(testCase,result.signals(1).groups);
verifyEqual(testCase,height(result.events),0);
verifyTrue(testCase,all(isnan([result.signals(1).fits.slope])));
end

function testMissingInputsAndOptions(testCase)
folder = makeSession(testCase,'20260921-Test-Missing');
[~,sessionName] = fileparts(folder);
verifyError(testCase,@() analyzeSessions_clampCalibration(folder, ...
    plotPhotometry=false,calibration=struct('emaTau',0)), ...
    'analyzeSessions_clampCalibration:UnknownCalibrationOption');
% Missing sync data must produce an actionable error, never assume the
% devices started together.
syncFile = fullfile(folder,['sync_',sessionName,'.mat']);
loaded = load(syncFile);
params = loaded.params;
params.sync = rmfield(params.sync,'timePhotometry');
save(syncFile,'params','-append');
verifyError(testCase,@() analyzeSessions_clampCalibration(folder,plotPhotometry=false), ...
    'analyzeSessions_clampCalibration:MissingSync');
delete(fullfile(folder,['timeseries_',sessionName,'.mat']));
verifyError(testCase,@() analyzeSessions_clampCalibration(folder,plotPhotometry=false), ...
    'analyzeSessions_clampCalibration:MissingTimeSeries');
end

function testPlotsOneFigurePerChannelWithEveryPower(testCase)
folder = makeSession(testCase,'20260921-Test-Plot');
close all
addTeardown(testCase,@() close('all'));
analyzeSessions_clampCalibration(folder,plotPhotometry=true);
figures = findall(groot,'Type','figure');
verifyEqual(testCase,numel(figures),2); % one per photometry channel
for i = 1:numel(figures)
    traces = findobj(figures(i),'Type','line','-regexp','DisplayName','% \(n=');
    verifyEqual(testCase,sort(string({traces.DisplayName})), ...
        sort(["20% (n=2)","60% (n=1)","25% (n=2)","75% (n=1)"]));
    % Both power-response panels carry their own fitted line.
    fitted = findobj(figures(i),'Type','line','-regexp','DisplayName','^Fit: ');
    verifyEqual(testCase,numel(fitted),2);
end
verifyTrue(testCase,isfile(fullfile(folder,'Calibration_photometry_power_response_01_green.pdf')));
verifyTrue(testCase,isfile(fullfile(folder,'Calibration_photometry_power_response_02_red.pdf')));
end

%% Fixtures

function folder = makeSession(testCase,name,options)
arguments
    testCase
    name (1,1) string
    options.photometryStart (1,1) double = 0
    options.extraChannels (1,1) logical = false
end
[blueClamp,redClamp,Fs] = makeSweeps;
folder = char(fullfile(tempname,name));
mkdir(folder);
addTeardown(testCase,@() rmdir(fileparts(folder),'s'));
[~,sessionName] = fileparts(folder);
% Deliberately omit labjack, NI photometry, behavior, and clampTarget.
save(fullfile(folder,['data_',sessionName,'.mat']),'blueClamp','redClamp');

finalFs = 100;
responses = [0.1 0.1 0.3 -0.1 -0.1 -0.3];
% LabJack channels start photometryStart seconds after NI, so their own
% responses sit that much earlier in the stored trace.
timeSeries = makeChannel('green','LJ',responses,finalFs,options.photometryStart);
timeSeries(2) = makeChannel('red','LJ',-responses,finalFs,options.photometryStart);
if options.extraChannels
    timeSeries(3) = makeChannel('PMT','NI',responses,finalFs,0);
    timeSeries(4) = makeChannel('blueClamp','NI',zeros(1,6),finalFs,0);
    timeSeries(5) = makeChannel('redClamp','NI',zeros(1,6),finalFs,0);
    timeSeries(6) = makeChannel('eyeArea','Cam',zeros(1,6),finalFs,0);
end
save(fullfile(folder,['timeseries_',sessionName,'.mat']),'timeSeries');

labjackFs = 2000;
params.sync = struct('behaviorFs',Fs,'labjackFs',labjackFs,'photometryFs',finalFs, ...
    'timeNI',(0:numel(blueClamp)-1)/Fs, ...
    'timePhotometry',options.photometryStart+(0:65*labjackFs-1)/labjackFs);
params.session.task = 'random';
params.analyze = struct('legacySetting',42); % no behavioral analysis options
save(fullfile(folder,['sync_',sessionName,'.mat']),'params');
end

function [blue,red,Fs] = makeSweeps
Fs = 1000;
blue = zeros(65*Fs,1); red = blue;
powers = [20 20.4 60 25 25.5 75];
for trial = 1:6
    idx = (5+(trial-1)*10)*Fs+(1:5*Fs);
    if trial <= 3
        red(idx) = laserVoltage(powers(trial),[25 500]);
    else
        blue(idx) = laserVoltage(powers(trial),[800 1600]);
    end
end
end

function voltage = laserVoltage(power,range)
voltage = (range(1)+(range(2)-range(1))*power/100)*5/4095;
end

function entry = makeChannel(name,system,amplitudes,finalFs,startSec)
% Mirrors the fields loadSessions stores for a processed photometry channel.
% startSec is when this system started relative to NI, so the pulse response
% lands at (NI onset - startSec) in the stored trace.
data = zeros(1,65*finalFs);
for trial = 1:numel(amplitudes)
    idx = round((5+(trial-1)*10-startSec)*finalFs)+(1:5*finalFs);
    data(idx) = amplitudes(trial);
end
entry = struct('name',name,'data',data,'finalFs',finalFs,'system',system, ...
    'time_offset',NaN,'demux',false,'demux_freq',NaN,'detrend',true, ...
    'detrend_type','dff','detrend_window',180,'options',struct());
end
