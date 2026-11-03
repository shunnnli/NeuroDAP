function tests = testDmdFitControlWindow

tests = functiontests(localfunctions);

end


function testTakesTheLastAnalysisLengthOfALongerControlWindow(testCase)

ctrl = repmat(1:400,3,1);

[data,time] = dmdFitControlWindow(ctrl,301,10);

verifyEqual(testCase,size(data),[3 301]);
verifyEqual(testCase,data(1,:),100:400);
verifyEqual(testCase,time(end),0);
verifyEqual(testCase,time(1),-30);

end


function testKeepsTheWholeControlWindowWhenItIsShorterThanTheAnalysisWindow(testCase)

% The real case since baselineGuard trimmed the stored control window:
% 291 control samples against a 301-sample analysis window. The old
% expression asked for column -9 and threw.
ctrl = repmat(1:291,3,1);

[data,time] = dmdFitControlWindow(ctrl,301,10);

verifyEqual(testCase,size(data),[3 291]);
verifyEqual(testCase,numel(time),291);
verifyEqual(testCase,time(end),0);
verifyEqual(testCase,time(1),-29);

end


function testTimeIsSpacedBySampleRateNotStretchedToFit(testCase)

[~,time] = dmdFitControlWindow(zeros(1,291),301,10);

verifyEqual(testCase,diff(time(1:2)),0.1,AbsTol=1e-12);

end


function testEmptyControlDataGivesEmptyOutputs(testCase)

[data,time] = dmdFitControlWindow(zeros(0,0),301,10);

verifyEmpty(testCase,data);
verifyEmpty(testCase,time);

end


function testWindowCanEndBeforeTheStimulusWhenAGuardIsGiven(testCase)

% The stored control window stops baselineGuard short of the pulse, so its
% last sample is at -guard ms, not at 0. Drawing it as if it ended at 0 put
% the whole baseline trace 1 ms late and closed up the guard interval.
[~,time] = dmdFitControlWindow(zeros(1,291),301,10,-1);

verifyEqual(testCase,time(end),-1,AbsTol=1e-12);
verifyEqual(testCase,time(1),-30,AbsTol=1e-12);

end


function testTheGuardShiftsWithoutChangingSpacing(testCase)

[~,time] = dmdFitControlWindow(zeros(1,291),301,10,-1);

verifyEqual(testCase,diff(time(1:2)),0.1,AbsTol=1e-12);
verifyEqual(testCase,numel(time),291);

end


function testOmittingTheGuardStillEndsAtTheStimulus(testCase)

% Data written before baselineGuard existed ran up to the pulse.
[~,withoutGuard] = dmdFitControlWindow(zeros(1,291),301,10);
[~,withZero]     = dmdFitControlWindow(zeros(1,291),301,10,0);

verifyEqual(testCase,withoutGuard(end),0);
verifyEqual(testCase,withoutGuard,withZero);

end
