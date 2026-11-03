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


%% Detection metrics (the notebook's HOTSPOT_TYPE = "evoked")

function testResponsiveSpotScoresAboveFlatOnes(testCase)

timeMs = -30:0.1:100; c = localCriteria();
spots = arrayfun(@(k) localNoiseSpot(timeMs,k),1:8,'UniformOutput',false);
ref = dmdHotspotVerdict('reference',spots,timeMs,c);

flat = dmdHotspotVerdict('score',spots{1},timeMs,c,ref);
loud = dmdHotspotVerdict('score',localRespondingSpot(timeMs,-40),timeMs,c,ref);

verifyGreaterThan(testCase,loud.chargeExcessZ,flat.chargeExcessZ);
verifyGreaterThan(testCase,loud.chargeExcessZ,c.minChargeZ);

end


function testSpreadIsNaNBelowMinimumSpots(testCase)

% Review Focus 3: below 5 spots a MAD is noise, so no spot is a hotspot.
timeMs = -30:0.1:100; c = localCriteria();
spots = arrayfun(@(k) localNoiseSpot(timeMs,k),1:3,'UniformOutput',false);
ref = dmdHotspotVerdict('reference',spots,timeMs,c);

verifyTrue(testCase,isnan(ref.chargeSpread));
verifyTrue(testCase,isnan(dmdHotspotVerdict('score',spots{1},timeMs,c,ref).chargeExcessZ));

end


function testEmptyTracesScoreNaN(testCase)

% Review Focus 1: a spot whose sweeps were all dropped by the Rs filter.
timeMs = -30:0.1:100; c = localCriteria();
spots = arrayfun(@(k) localNoiseSpot(timeMs,k),1:8,'UniformOutput',false);
ref = dmdHotspotVerdict('reference',spots,timeMs,c);

out = dmdHotspotVerdict('score',zeros(0,numel(timeMs)),timeMs,c,ref);

verifyTrue(testCase,isnan(out.chargeExcessZ));
verifyEqual(testCase,out.nSweeps,0);

end


function testSweepsAreBaselineCorrectedBeforeAveraging(testCase)

timeMs = -30:0.1:100; c = localCriteria();
base = localRespondingSpot(timeMs,-40);
offset = base + (1:size(base,1))'*25;   % per-sweep DC offsets
spots = arrayfun(@(k) localNoiseSpot(timeMs,k),1:8,'UniformOutput',false);
ref = dmdHotspotVerdict('reference',spots,timeMs,c);

verifyEqual(testCase, ...
    dmdHotspotVerdict('score',offset,timeMs,c,ref).chargeExcessZ, ...
    dmdHotspotVerdict('score',base,timeMs,c,ref).chargeExcessZ,'RelTol',1e-6);

end


%% Shared criteria

function testDefaultCriteriaCarrySpecValues(testCase)

c = dmdDefaultCriteria();

verifyEqual(testCase,c.baselineMs,[-30 -1]);
verifyEqual(testCase,c.reanalysisWindowMs,[0 50]);
verifyEqual(testCase,c.minAmplitudePa,2);
verifyEqual(testCase,c.maxRs,30);
verifyEqual(testCase,c.minChargeZ,0.5);
verifyEqual(testCase,c.minSweeps,2);

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


function traces = localNoiseSpot(timeMs,seed)

% One signal-free spot; the seed makes each of them a different draw.
rng(seed,'twister');
traces = 0.3*randn(5,numel(timeMs));

end


function traces = localRespondingSpot(timeMs,amplitudePa)

rng(99,'twister');
traces = 0.3*randn(5,numel(timeMs));
responseMask = timeMs >= 5 & timeMs <= 25;
traces(:,responseMask) = traces(:,responseMask) + amplitudePa;

end
