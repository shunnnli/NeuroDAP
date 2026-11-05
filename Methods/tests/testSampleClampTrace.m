function tests = testSampleClampTrace
tests = functiontests(localfunctions);
end

function testNativeNIReferenceForDownsampledTraces(testCase)
[signals,params] = exampleSignals();
sample = getSampleClampTrace(signals,params,5000,[-1,2],removeClampBaseline=false);
verifyEqual(testCase,sample.timestampDA,-1:0.1:2,'AbsTol',1e-12);
verifyEqual(testCase,sample.DA,signals(1).data(40:70));
verifyEqual(testCase,sample.timestampClamp,-1:0.02:2,'AbsTol',1e-12);
verifyEqual(testCase,sample.clampIndex,sample.redClamp-sample.blueClamp,'AbsTol',1e-12);
verifyEqual(testCase,sample.referenceTime_sec,4.999,'AbsTol',1e-12);
verifyEqual(testCase,sample.sessionTimeDA_sec,sample.referenceTime_sec+sample.timestampDA);
verifyEqual(testCase,sample.DAUnits,'z-score');
end

function testLabjackOffsetUsesSharedSyncClock(testCase)
[signals,params] = exampleSignals();
signals(1).system = 'LJ';
params.sync.labjackFs = 100;
% NI starts at -2 s, LJ at -1 s: LJ recording begins one second later.
params.sync.timePhotometry = -1+(0:1999)/100;
sample = getSampleClampTrace(signals,params,5001,[-1,2],removeClampBaseline=false);
verifyEqual(testCase,sample.DA,signals(1).data(30:60));
verifyEqual(testCase,sample.referenceTime_sec,5,'AbsTol',1e-12);
end

function testBaselineIncludedOutsideDisplayedWindow(testCase)
[signals,params] = exampleSignals();
% Clamp offsets differ in the pre-event baseline and should cancel separately.
signals(2).data(:) = 10;
signals(3).data(:) = 30;
signals(2).data(250:end) = 15;
signals(3).data(250:end) = 50;
sample = getSampleClampTrace(signals,params,5000,[0,1]);
verifyEqual(testCase,sample.clampIndex,15*ones(size(sample.clampIndex)),'AbsTol',1e-12);
verifyEqual(testCase,sample.blueClamp,15*ones(size(sample.blueClamp)));
verifyEqual(testCase,sample.redClamp,50*ones(size(sample.redClamp)));
verifyEqual(testCase,sample.timestampClamp,0:0.02:1,'AbsTol',1e-12);
end

function testRejectsEdgeWindowInsteadOfReturningEmptyExample(testCase)
[signals,params] = exampleSignals();
verifyError(testCase,@() getSampleClampTrace(signals,params,100,[-1,2]), ...
    'getSampleClampTrace:WindowOutsideRecording');
verifyError(testCase,@() getSampleClampTrace(signals,params,19500,[-1,2]), ...
    'getSampleClampTrace:WindowOutsideRecording');
end

function testRejectsMissingCommandsAndMismatchedSampling(testCase)
[signals,params] = exampleSignals();
verifyError(testCase,@() getSampleClampTrace(signals(1:2),params,5000,[-1,2]), ...
    'getSampleClampTrace:MissingClamp');
signals(3).finalFs = 25;
verifyError(testCase,@() getSampleClampTrace(signals,params,5000,[-1,2]), ...
    'getSampleClampTrace:ClampSamplingMismatch');
end

function [signals,params] = exampleSignals()
params.session.baselineSystem = 'NI';
params.sync.behaviorFs = 1000;
params.sync.timeNI = -2+(0:19999)/1000;
signals = struct('name',{'DA','blueClamp','redClamp'}, ...
    'data',{1:200,10+sin((1:1000)/50),30+cos((1:1000)/50)}, ...
    'finalFs',{10,50,50},'system',{'NI','NI','NI'}, ...
    'detrend_type',{'rolling-z','',''});
end
