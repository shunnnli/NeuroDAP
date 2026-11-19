function tests = testSpontaneousClampActivity
tests = functiontests(localfunctions);
end

function testUsesEitherCommandAndExcludesUnknown(testCase)
[animals,~] = fixture;
animals(2).data(:)=0; animals(3).data(:)=0;
animals(2).data(2)=1; % threshold equality is OFF
animals(2).data(3)=2;
animals(3).data(4)=2;
animals(3).data(5)=NaN;
[active,blue,red]=getTrialClampActivity(animals,1,1);
verifyEqual(testCase,active(1:4),[0,0,1,1]);
verifyTrue(testCase,isnan(active(5)));
verifyEqual(testCase,[blue,red],[2,3]);
end

function testThreeSecondExclusionAndStateBoundaries(testCase)
[animals,cfg]=fixture;
% Both OFF and ON activity deliberately live in a '(clamp)' labeled row.
% A pre-3 s peak and a peak crossing OFF->ON must not be counted.
t=animals(1).timestamp; x=zeros(size(t));
for pair=[2.7,2.8; 3.1,3.2; 3.45,3.5; 3.7,3.8]'
    [~,lo]=min(abs(t-pair(1))); [~,hi]=min(abs(t-pair(2)));
    x(lo)=-0.01; x(hi)=0.1;
end
animals(1).data=x;
[summary,~,windows]=measureSpontaneousTransients(animals,cfg);
verifyEqual(testCase,[summary.OffCount,summary.OnCount],[1,1]);
verifyTrue(testCase,all(windows.StartTime>3));
verifyTrue(testCase,all(windows.EndTime<=4.5+1e-10));
verifyEqual(testCase,sum(windows.Seconds),sum(t>3+1e-9 & t<4.5-1e-9)*0.05,'AbsTol',1e-10);
verifyTrue(testCase,all(windows.EndTime(windows.Condition==1)<=3.5+1e-10));
verifyTrue(testCase,all(windows.StartTime(windows.Condition==2)>=3.5-1e-10));
end

function testRejectsMisalignedTrialMetadata(testCase)
[animals,~]=fixture;
animals(2).trialInfo.trialNumber=99;
verifyError(testCase,@()getTrialClampActivity(animals,1,1),'ClampQuality:CommandTrialMismatch');
end

function [animals,cfg]=fixture
t=0:0.05:4.5;
row=struct('animal','test','task','Random','event','Water (clamp)', ...
    'name','NAc-clamp','timestamp',t,'data',zeros(size(t)));
row.trialInfo=struct('trialNumber',1,'trialTable', ...
    table(1,1,10,0,'VariableNames',{'performing','isReward','ITI','CueTime'}));
row.options.sessionList={'test-session'};
animals=repmat(row,1,3);
animals(2).name='blueClamp'; animals(2).data(t>=3.5)=2;
animals(3).name='redClamp';
cfg=struct('animals','All','scale',100,'minimumRisePercent',3, ...
    'rewardExclusionSeconds',3,'nextCueGuardSeconds',2,'clampActiveThresholdPercent',1);
cfg.spontaneousSelections=struct('task','Random','signal','NAc-clamp', ...
    'event','Water (clamp)','window',[3,4.5],'rewardAligned',true, ...
    'cueToRewardSeconds',0,'confirmedEventFree',false);
end
