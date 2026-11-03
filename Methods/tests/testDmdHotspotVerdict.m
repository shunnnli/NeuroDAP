function tests = testDmdHotspotVerdict

tests = functiontests(localfunctions);

end


%% Reanalysis test (the notebook's HOTSPOT_TYPE = "reanalysis")

function testFlatTraceIsNotHotspot(testCase)

[timeMs,traces] = localSyntheticSpot(0);
out = dmdHotspotVerdict('reanalysis',traces,timeMs,localCriteria());

verifyFalse(testCase,out.isHotspot);

end


function testLargeSustainedResponseIsHotspot(testCase)

[timeMs,traces] = localSyntheticSpot(-40);
out = dmdHotspotVerdict('reanalysis',traces,timeMs,localCriteria());

verifyTrue(testCase,out.isHotspot);

end


function testAmplitudeFloorRejectsTinyButSignificant(testCase)

% Clears 3*SD comfortably but stays under the 2 pA floor.
[timeMs,traces] = localSyntheticSpot(-1);

c = localCriteria(); c.minAmplitudePa = 2;
verifyFalse(testCase,dmdHotspotVerdict('reanalysis',traces,timeMs,c).isHotspot);

c.minAmplitudePa = 0;
verifyTrue(testCase,dmdHotspotVerdict('reanalysis',traces,timeMs,c).isHotspot);

end


function testSdComesFromBaselineWindowOnly(testCase)

[timeMs,traces] = localSyntheticSpot(-40);
out = dmdHotspotVerdict('reanalysis',traces,timeMs,localCriteria());

verifyLessThan(testCase,out.baselineSd,5);   % the response must not inflate it

end


function testNaNBaselineSamplesAreIgnored(testCase)

% Review Focus 5: a dropped sample must not take the SD with it.
[timeMs,traces] = localSyntheticSpot(-40);
clean = dmdHotspotVerdict('reanalysis',traces,timeMs,localCriteria()).baselineSd;
traces(1,10) = NaN;

verifyEqual(testCase, ...
    dmdHotspotVerdict('reanalysis',traces,timeMs,localCriteria()).baselineSd, ...
    clean,'RelTol',0.05);

end


%% Shared criteria

function testDefaultCriteriaCarrySpecValues(testCase)

c = dmdDefaultCriteria();

verifyEqual(testCase,c.baselineMs,[-30 -1]);
verifyEqual(testCase,c.reanalysisWindowMs,[0 50]);
verifyEqual(testCase,c.minAmplitudePa,2);
verifyEqual(testCase,c.maxRs,30);
verifyEqual(testCase,c.minChargeZ,0.5);

end


%% Helpers

function c = localCriteria()

% Delegates so the tests and loadSlicesDMD cannot drift apart.
c = dmdDefaultCriteria();

end


function [timeMs,traces] = localSyntheticSpot(amplitudePa)

% Five sweeps on the stored window, Gaussian noise plus a step over 5-25 ms.
rng(0,'twister');
timeMs = -30:0.1:100;
traces = 0.3*randn(5,numel(timeMs));
responseMask = timeMs >= 5 & timeMs <= 25;
traces(:,responseMask) = traces(:,responseMask) + amplitudePa;

end
