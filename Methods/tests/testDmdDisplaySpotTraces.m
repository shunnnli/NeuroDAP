function tests = testDmdDisplaySpotTraces
tests = functiontests(localfunctions);
end

function testExcludedSpotKeepsOriginalSweepsForGrayDisplay(testCase)
raw = {[1 2;3 4];[5 6;7 8];[]};
shown = dmdDisplaySpotTraces(raw,{[true;false];[false;false];false(0,1)});
verifyEqual(testCase,shown{1},raw{1}(1,:));
verifyEqual(testCase,shown{2},raw{2});
verifyEmpty(testCase,shown{3});
end
