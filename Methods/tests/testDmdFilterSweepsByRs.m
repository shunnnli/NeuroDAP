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
verifyEqual(testCase,dmdFilterSweepsByRs(NaN,NaN,30),false);

end
