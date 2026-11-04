function tests = testDmdSelectSweepsByRs
tests = functiontests(localfunctions);
end

function testOnePassingSweepExpandsTo35Inclusive(testCase)
rs = [20;31;35;35.1;40];
verifyEqual(testCase,dmdSelectSweepsByRs(rs,rs,30,35),[true;true;true;false;false]);
end

function testZeroPassingSweepsDoesNotRelax(testCase)
rs = [31;34;40];
verifyEqual(testCase,dmdSelectSweepsByRs(rs,rs,30,35),false(3,1));
end

function testTwoPassingSweepsDoesNotRelax(testCase)
rs = [20;25;34];
verifyEqual(testCase,dmdSelectSweepsByRs(rs,rs,30,35),[true;true;false]);
end

function testUsesMinimumOfComputedAndHeaderRs(testCase)
verifyEqual(testCase,dmdSelectSweepsByRs([20;40;40],[22;34;36],30,35),[true;true;false]);
end
