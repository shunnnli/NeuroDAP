function tests = testDmdWindowDensity

tests = functiontests(localfunctions);

end


function testConstantOffsetContributesNothing(testCase)

timeMs = 0:0.1:45;
trace = 7*ones(size(timeMs));
[q,g] = dmdWindowDensity(trace,timeMs,[0 45],2);

verifyEqual(testCase,q,0,'AbsTol',1e-12);
verifyEqual(testCase,g,0,'AbsTol',1e-12);

end


function testChargeDensityIsPerMillisecond(testCase)

timeMs = 0:0.1:45;
trace = zeros(size(timeMs));
trace(timeMs>=10 & timeMs<20) = -20;
[q,~] = dmdWindowDensity(trace,timeMs,[0 45],2);

verifyGreaterThan(testCase,q,0);   % rectified, so positive

end


function testNaNSamplesAreIgnored(testCase)

% Review Focus 5: a dropped sample must not take the whole metric with it.
timeMs = 0:0.1:45;
trace = zeros(size(timeMs));
trace(timeMs>=10 & timeMs<20) = -20;
clean = dmdWindowDensity(trace,timeMs,[0 45],2);
trace(5) = NaN;

verifyEqual(testCase,dmdWindowDensity(trace,timeMs,[0 45],2),clean,'RelTol',0.01);

end
