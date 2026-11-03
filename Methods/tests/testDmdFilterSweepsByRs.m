function tests = testDmdFilterSweepsByRs

tests = functiontests(localfunctions);

end


function testJudgesTheMinimumOfBothEstimates(testCase)

keep = dmdFilterSweepsByRs([50 20 28],[25 40 29],30);

verifyEqual(testCase,keep,[true true true]);

end


function testDropsWhenBothEstimatesExceedLimit(testCase)

verifyEqual(testCase,dmdFilterSweepsByRs(50,45,30),false);

end


function testNonPositiveHeaderIsMissing(testCase)

% Review Focus 2: a header written as a flat 0 is not a measurement, so it
% must not count as an unbeatably small Rs.
verifyEqual(testCase,dmdFilterSweepsByRs(45,0,30),false);
verifyEqual(testCase,dmdFilterSweepsByRs(20,-3,30),true);

end


function testNaNComputedFallsBackToHeader(testCase)

verifyEqual(testCase,dmdFilterSweepsByRs(NaN,25,30),true);

end


function testSweepWithNoEstimateIsKept(testCase)

% pyNeuroDAP slice.py apply_qc_limits: "A sweep whose metric is missing even
% after the header fallback passes: the limits remove sweeps that are known to
% be bad, not sweeps of unknown quality." Reconstructed searches store a QC
% struct with no Rs at all, so dropping these would blank a whole search.
verifyEqual(testCase,dmdFilterSweepsByRs(NaN,NaN,30),true);
verifyEqual(testCase,dmdFilterSweepsByRs(NaN,0,30),true);
verifyEqual(testCase,dmdFilterSweepsByRs([NaN 50 20],[NaN 45 25],30),[true false true]);

end
