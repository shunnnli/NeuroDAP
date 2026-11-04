function tests = testDmdSweepMeanStats

tests = functiontests(localfunctions);

end


function testPeakIsAveragedAroundItsOwnSample(testCase)

timeMs = -30:0.1:100;
trace = zeros(size(timeMs));
trace(timeMs>=10 & timeMs<=12) = -30;
trace(abs(timeMs-11)<0.05) = -100;          % one noisy sample
s = dmdSweepMeanStats(trace,timeMs,[0 50],1);

verifyGreaterThan(testCase,s.min,-60);      % the 1 ms average blunts it
verifyEqual(testCase,s.minTime,11,'AbsTol',0.6);

end


function testEmptyTracesGiveNaN(testCase)

% Review Focus 1: a spot whose sweeps were all dropped by the Rs filter.
s = dmdSweepMeanStats(zeros(0,1301),-30:0.1:100,[0 50],1);

verifyTrue(testCase,isnan(s.min) && isnan(s.max) && isnan(s.auc));

end


function testBoundaryPeakAverageUsesFullTrace(testCase)

timeMs = -30:0.1:100;
trace = zeros(size(timeMs));
trace(timeMs >= 29.5 & timeMs <= 30.5) = 12;
trace(timeMs == 30) = 24;
s = dmdSweepMeanStats(trace,timeMs,[0 30],1);
mask = timeMs >= s.maxTime-0.5 & timeMs <= s.maxTime+0.5;
verifyEqual(testCase,s.max,mean(trace(mask)),'AbsTol',1e-12);

end


function testBaselineWindowIsOverridable(testCase)

% loadSlicesDMD's control window follows controlWindowLength and
% baselineGuard, so the levelling window cannot be hardcoded here.
timeMs = -30:0.1:100;
trace = zeros(size(timeMs));
trace(timeMs >= -30 & timeMs < -10) = 5;
trace(timeMs >= -10 & timeMs <= -1) = 20;

sDefault = dmdSweepMeanStats(trace,timeMs,[0 50],1);              % [-30 -1] -> level 5
sNarrow  = dmdSweepMeanStats(trace,timeMs,[0 50],1,[-10 -1]);     % level 20

% The window sits at 0 either way, so the two areas are in the ratio of the
% two levels. Asserting the ratio keeps this independent of how many samples
% the inclusive mask happens to catch at its edges.
verifyLessThan(testCase,sDefault.auc,0);
verifyEqual(testCase,sNarrow.auc/sDefault.auc,4,'RelTol',1e-9);

end


function testPeakAndAucAreInvariantToPerSweepBaselineOffsets(testCase)

timeMs = -30:0.1:100;
trace = -10*exp(-((timeMs-10)/2).^2) + 5*exp(-((timeMs-20)/2).^2);
traces = [trace;2*trace];
expected = dmdSweepMeanStats(traces,timeMs,[0 30],1);
actual = dmdSweepMeanStats(traces+[100;-70],timeMs,[0 30],1);

verifyEqual(testCase,actual.max,expected.max,'AbsTol',1e-10);
verifyEqual(testCase,actual.min,expected.min,'AbsTol',1e-10);
verifyEqual(testCase,actual.auc,expected.auc,'AbsTol',1e-10);
verifyEqual(testCase,actual.max+actual.min,expected.max+expected.min,'AbsTol',1e-10);

end
