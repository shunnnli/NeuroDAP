function tests = testDmdTailSamples

tests = functiontests(localfunctions);

end


function testTakesTheLastSamplesWhenPlentyAreAvailable(testCase)

idx = dmdTailSamples(291,201);

verifyEqual(testCase,idx,91:291);

end


function testKeepsWhatExistsWhenFewerSamplesAreAvailable(testCase)

% The crash this guards: a 1-sample control window asked for its last 201
% produced the index -199, which MATLAB rejects outright.
idx = dmdTailSamples(1,201);

verifyEqual(testCase,idx,1);

end


function testExactFitReturnsEverything(testCase)

verifyEqual(testCase,dmdTailSamples(201,201),1:201);

end


function testNoSamplesGivesAnEmptyIndex(testCase)

verifyEmpty(testCase,dmdTailSamples(0,201));

end


function testSingleSampleRequestTakesOnlyTheLast(testCase)

verifyEqual(testCase,dmdTailSamples(291,1),291);

end
