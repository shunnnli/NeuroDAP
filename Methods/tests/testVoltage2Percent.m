function tests = testVoltage2Percent
tests = functiontests(localfunctions);
end

%% voltage2percent: nominal range

function testNominalRangeUnchanged(testCase)
% A command topping out at the calibration full scale (500 counts) read
% against the nominal [25,600] range is the 83% bug, and stays that way
% while autoMax is off (the default).
command = countsToVoltage([25 200 500]);
pct = voltage2percent(command,[25 600]);
verifyEqual(testCase,pct,[0 30.4348 82.6087],'AbsTol',1e-4);
[~,usedRange] = voltage2percent(command,[25 600]);
verifyEqual(testCase,usedRange,[25 600]);
end

%% voltage2percent: autoMax

function testAutoMaxMapsSessionMaxTo100(testCase)
command = countsToVoltage([25 200 500]);
[pct,usedRange] = voltage2percent(command,[25 600],autoMax=true);
verifyEqual(testCase,pct,[0 36.8421 100],'AbsTol',1e-4);
verifyEqual(testCase,usedRange,[25 500],'AbsTol',1e-9);
end

function testAutoMaxWorksForBriefFullScale(testCase)
% The full scale need only appear briefly within the session.
counts = 25*ones(1,10000); counts(4001:4010) = 500;
[pct,usedRange] = voltage2percent(countsToVoltage(counts),[25 600],autoMax=true);
verifyEqual(testCase,max(pct),100,'AbsTol',1e-9);
verifyEqual(testCase,usedRange(2),500,'AbsTol',1e-9);
end

function testAutoMaxRecoversClippedBlueRange(testCase)
% Blue commands above the nominal 1500 are clipped to 100% without autoMax,
% so different levels become indistinguishable.
counts = [800*ones(1,100), 1550*ones(1,100), 1600*ones(1,100)];
pctNominal = voltage2percent(countsToVoltage(counts),[800 1500]);
verifyEqual(testCase,max(pctNominal),100,'AbsTol',1e-9);
verifyEqual(testCase,pctNominal(150),pctNominal(250),'AbsTol',1e-9); % both clipped
[pct,usedRange] = voltage2percent(countsToVoltage(counts),[800 1500],autoMax=true);
verifyEqual(testCase,usedRange(2),1600,'AbsTol',1e-9);
verifyEqual(testCase,pct([50 150 250]),[0 93.75 100],'AbsTol',1e-4);
end

function testAutoMaxIgnoresSingleSpike(testCase)
% One stray sample above the plateau must not set the full scale.
counts = 25*ones(1,10000); counts(2001:3000) = 500; counts(7777) = 590;
[~,usedRange] = voltage2percent(countsToVoltage(counts),[25 600],autoMax=true);
verifyEqual(testCase,usedRange(2),500,'AbsTol',1e-9);
% Unless enough samples sit at that level to count as a real command.
counts(7777:7781) = 590;
[~,usedRange] = voltage2percent(countsToVoltage(counts),[25 600],autoMax=true);
verifyEqual(testCase,usedRange(2),590,'AbsTol',1e-9);
end

function testAutoMaxFallsBackWhenClampNeverOn(testCase)
% Nothing but the off level: keep the nominal range instead of mapping
% baseline noise onto 0-100%.
counts = 25 + rand(1,1000);
verifyWarning(testCase,@() voltage2percent(countsToVoltage(counts),[25 600],autoMax=true),...
    'voltage2percent:AutoMaxFailed');
warning('off','voltage2percent:AutoMaxFailed');
cleanup = onCleanup(@() warning('on','voltage2percent:AutoMaxFailed'));
[pct,usedRange] = voltage2percent(countsToVoltage(counts),[25 600],autoMax=true);
verifyEqual(testCase,usedRange,[25 600]);
verifyLessThan(testCase,max(pct),1);
end

function testAutoMaxHandlesNaNAndShape(testCase)
counts = [25 NaN 500; 200 500 NaN];
[pct,usedRange] = voltage2percent(countsToVoltage(counts),[25 600],autoMax=true,...
                                  autoMaxMinSamples=2);
verifyEqual(testCase,size(pct),size(counts));
verifyEqual(testCase,usedRange(2),500,'AbsTol',1e-9);
verifyTrue(testCase,all(isnan(pct([3 6]))));
verifyEqual(testCase,pct([1 2]),[0 36.8421],'AbsTol',1e-4);
end

%% Helper

function voltage = countsToVoltage(counts)
voltage = counts*5/4095;
end
