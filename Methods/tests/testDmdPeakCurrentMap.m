function tests = testDmdPeakCurrentMap
tests = functiontests(localfunctions);
end

function testUsesMeanTraceAndRejectsBadSweep(testCase)
t = -30:0.1:100;
a = zeros(size(t));
a(t>=9 & t<=11) = 8;
a(t>=19 & t<=21) = -12;
a(t==10) = 19;  % unique peak centres; 11 samples in each 1 ms window
a(t==20) = -23;
% DC offsets must cancel, and an excluded large response must not enter.
traces = [a+100; a+200; 100*a];
map = dmdPeakCurrentMap({traces},{[true;true;false]},[1 2 1 2], ...
    [2 3],t,[0 30],[-30 -1]);
verifyEqual(testCase,map(:,1:2),-4*ones(2),'AbsTol',1e-12);
verifyTrue(testCase,all(isnan(map(:,3))));
end

function testOneSweepHasAmplitudeButNoUsableSweepIsMissing(testCase)
t = -30:0.1:100;
a = zeros(size(t));
a(t>=9 & t<=11) = -10;
a(t==10) = -21;
map = dmdPeakCurrentMap({a;a},{true;false},[1 1 1 1;2 2 1 1], ...
    [1 3],t,[0 30],[-30 -1]);
verifyEqual(testCase,map(1),-11,'AbsTol',1e-12);
verifyTrue(testCase,all(isnan(map(2:3))));
end

function testAveragesSweepsBeforeFindingPeaks(testCase)
t = -30:0.1:100;
a = zeros(size(t));
a(t>=9 & t<=11) = 20;
map = dmdPeakCurrentMap({[a;-a]},{[true;true]},[1 1 1 1], ...
    [1 1],t,[0 30],[-30 -1]);
verifyEqual(testCase,map,0,'AbsTol',1e-12);
end

function testPolarityReversesSignedCurrent(testCase)
t = -30:0.1:100;
a = zeros(size(t));
a(t>=9 & t<=11) = -10;
a(t==10) = -21;
map = dmdPeakCurrentMap({a;-a},{true;true},[1 1 1 1;2 2 1 1], ...
    [1 2],t,[0 30],[-30 -1]);
verifyEqual(testCase,map,[-11 11],'AbsTol',1e-12);
end

function testAucUsesEqualWeightsAndLocalBaselines(testCase)
t = -30:0.1:100;
a = zeros(size(t));
a(t>=0 & t<=30) = 10;
traces = [a+100;3*a-70;100*a];
map = dmdPeakCurrentMap({traces},{[true;true;false]},[1 1 1 1], ...
    [1 1],t,[0 30],[-30 -1],"auc");
% Retained response mean is 20 pA over 301 samples at 10 kHz.
verifyEqual(testCase,map,20*301/10000,'AbsTol',1e-12);
end
