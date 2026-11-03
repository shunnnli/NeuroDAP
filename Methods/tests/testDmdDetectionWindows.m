function tests = testDmdDetectionWindows

tests = functiontests(localfunctions);

end


function testControlWindowsTileAfterGap(testCase)

timeMs = linspace(-30,100,1301);
[w,c] = dmdDetectionWindows(timeMs,[5 50],5);

verifyEqual(testCase,w,[5 50]);
verifyEqual(testCase,c,[55 100],'AbsTol',1e-9);

end


function testSecondWindowAppearsWhenRoomAllows(testCase)

timeMs = linspace(-30,200,2301);
[~,c] = dmdDetectionWindows(timeMs,[5 50],5);

verifyEqual(testCase,size(c,1),3);
verifyEqual(testCase,c(2,:),[100 145],'AbsTol',1e-9);

end


function testTooShortTraceErrors(testCase)

% Review Focus 4: a stored window that cannot hold a matched control window
% must say so rather than score against a truncated one.
timeMs = linspace(-10,50,601);

verifyError(testCase,@() dmdDetectionWindows(timeMs,[5 50],5), ...
            'dmdDetectionWindows:noControlWindow');

end
