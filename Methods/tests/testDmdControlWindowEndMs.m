function tests = testDmdControlWindowEndMs

tests = functiontests(localfunctions);

end


function testStoredGuardPutsTheWindowEndBeforeTheStimulus(testCase)

verifyEqual(testCase,dmdControlWindowEndMs(struct('baselineGuard',1)),-1);

end


function testAFractionalGuardIsHonoured(testCase)

verifyEqual(testCase,dmdControlWindowEndMs(struct('baselineGuard',0.5)),-0.5);

end


function testAZeroGuardEndsAtTheStimulus(testCase)

verifyEqual(testCase,dmdControlWindowEndMs(struct('baselineGuard',0)),0);

end


function testFilesWrittenBeforeTheGuardExistedEndAtTheStimulus(testCase)

% Older cells_DMD files have no baselineGuard field at all. Their control
% window really did run up to the pulse, so 0 is the right answer for them
% and their figures keep rendering as they always have.
verifyEqual(testCase,dmdControlWindowEndMs(struct('outputFs',10000)),0);

end
