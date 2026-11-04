function tests = testDmdControlMeanStats
tests = functiontests(localfunctions);
end

function testSeparateExtremaUseTimeAxisNotFlattenedIndices(testCase)
a = zeros(1,291); a(100)=10; a(200)=-20;
s = dmdControlMeanStats([a+100; a-70],10000,-1);
verifyEqual(testCase,s.maxTime,-20.1,'AbsTol',1e-12);
verifyEqual(testCase,s.minTime,-10.1,'AbsTol',1e-12);
verifyGreaterThan(testCase,s.max,0);
verifyLessThan(testCase,s.min,0);
end

function testExtremaAreFoundAfterAveragingSweeps(testCase)
a = zeros(1,291); a(100)=100;
s = dmdControlMeanStats([a;-a],10000,-1);
verifyEqual(testCase,s.max,0);
verifyEqual(testCase,s.min,0);
end

function testNoUsableSweepsHaveNoLatency(testCase)
s = dmdControlMeanStats(zeros(0,291),10000,-1);
verifyTrue(testCase,isnan(s.minTime) && isnan(s.maxTime));
end
