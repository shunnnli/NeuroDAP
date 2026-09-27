function tests = testRescaleClampTraces
tests = functiontests(localfunctions);
end

%% Summary struct

function testRescalesSummaryTracesAndStages(testCase)
% Red command topping out at 500 counts read against [25,600] peaks at 82.6%
summary = makeSummary(82.6087);
[summary,info] = rescaleClampTraces(summary,verbose=false);

verifyEqual(testCase,max(summary(1).data,[],'all'),100,'AbsTol',1e-9);
verifyEqual(testCase,max(summary(2).data,[],'all'),100,'AbsTol',1e-9);
% Zero stays zero, intermediate levels keep their ratio to full scale
verifyEqual(testCase,summary(1).data(1,1),0,'AbsTol',1e-9);
verifyEqual(testCase,summary(1).data(2,1),50,'AbsTol',1e-9);
% Stage stats derived from the trace are scaled by the same factor
verifyEqual(testCase,summary(1).stageAvg.data,[0;50;100],'AbsTol',1e-9);
verifyEqual(testCase,summary(1).stageArea.data,2*[0;50;100],'AbsTol',1e-9);
% Photometry rows are untouched
verifyEqual(testCase,summary(3).data,0.5*ones(3,10));
verifyEqual(testCase,summary(3).stageAvg.data,0.5*ones(3,1));

verifyEqual(testCase,height(info),2); % one row per session-channel
verifyEqual(testCase,info.status,["rescaled";"rescaled"]);
verifyEqual(testCase,info.detectedMax,[82.6087;82.6087],'AbsTol',1e-4);
verifyEqual(testCase,info.scaleFactor,[1.2105;1.2105],'AbsTol',1e-4);
verifyEqual(testCase,sort(info.name),["blueClamp";"redClamp"]);
end

function testMaxIsPooledAcrossEventsOfTheSameSession(testCase)
% Full scale may appear in only one event of a session; every row of that
% session & channel must still get the same scale factor.
summary = makeSummary(82.6087);
summary(4) = summary(1); summary(4).event = 'Water'; % this event reaches full scale
summary(1).data = [0;20;40].*ones(3,10); % this one peaks at 40%
[summary,info] = rescaleClampTraces(summary,verbose=false);

verifyEqual(testCase,info.nRows(info.name == "redClamp"),2);
verifyEqual(testCase,max(summary(1).data,[],'all'),48.4211,'AbsTol',1e-4);
verifyEqual(testCase,max(summary(4).data,[],'all'),100,'AbsTol',1e-9);
end

function testSeparateSessionsScaleIndependently(testCase)
summary = makeSummary(82.6087);
session2 = summary(1); session2.session = 's2'; session2.date = '20260716';
session2.data = [0;25;50].*ones(3,10); % this session peaks at 50%
summary(end+1) = session2;
[summary,info] = rescaleClampTraces(summary,verbose=false);

verifyEqual(testCase,height(info),3);
verifyEqual(testCase,max(summary(1).data,[],'all'),100,'AbsTol',1e-9);
verifyEqual(testCase,max(summary(end).data,[],'all'),100,'AbsTol',1e-9);
verifyEqual(testCase,summary(end).data(2,1),50,'AbsTol',1e-9);
end

%% Date filtering

function testDateBeforeSelectsSessions(testCase)
summary = makeSummary(82.6087);
late = summary(1); late.date = '20261001'; late.session = 's2';
summary(end+1) = late;
[summary,info] = rescaleClampTraces(summary,dateBefore=20260927,verbose=false);

verifyEqual(testCase,height(info),2); % only the 20260715 session
verifyEqual(testCase,max(summary(1).data,[],'all'),100,'AbsTol',1e-9);
verifyEqual(testCase,max(summary(end).data,[],'all'),82.6087,'AbsTol',1e-4);
end

%% Animals struct

function testAnimalsStructWithoutDateField(testCase)
animals = makeAnimals(82.6087);
verifyWarning(testCase,@() rescaleClampTraces(animals,dateBefore=20260927,verbose=false),...
    'rescaleClampTraces:NoDateField');
warning('off','rescaleClampTraces:NoDateField');
cleanup = onCleanup(@() warning('on','rescaleClampTraces:NoDateField')); %#ok<NASGU>

[animals,info] = rescaleClampTraces(animals,dateBefore=20260927,verbose=false);
verifyEqual(testCase,max(animals(1).data,[],'all'),100,'AbsTol',1e-9);
verifyEqual(testCase,animals(1).stageMax.data,[0;50;100],'AbsTol',1e-9);
verifyEqual(testCase,info.status,"rescaled");
verifyEqual(testCase,info.group,"SL431_Random-clamp_redClamp");
end

%% Edge cases

function testAlreadyFullScaleAndClampOffAreLeftAlone(testCase)
summary = makeSummary(100); % correctly scaled already (or clipped)
[scaled,info] = rescaleClampTraces(summary,verbose=false);
verifyEqual(testCase,scaled(1).data,summary(1).data);
verifyEqual(testCase,unique(info.status),"unchanged: already at targetPct");
verifyEqual(testCase,unique(info.scaleFactor),1);

summary = makeSummary(2); % clamp never turned on
[scaled,info] = rescaleClampTraces(summary,verbose=false);
verifyEqual(testCase,scaled(1).data,summary(1).data);
verifyEqual(testCase,unique(info.status),"skipped: peak below minPct");
end

function testSpikeDoesNotSetFullScale(testCase)
summary = makeSummary(82.6087);
summary(1).data(1,5) = 95; % single stray sample
[~,info] = rescaleClampTraces(summary,verbose=false);
verifyEqual(testCase,info.detectedMax(info.name == "redClamp"),82.6087,'AbsTol',1e-4);
end

function testIdempotentAndNoClampRows(testCase)
summary = makeSummary(82.6087);
once = rescaleClampTraces(summary,verbose=false);
twice = rescaleClampTraces(once,verbose=false);
verifyEqual(testCase,twice(1).data,once(1).data,'AbsTol',1e-12);

[unchanged,info] = rescaleClampTraces(summary(3),verbose=false); % photometry only
verifyEqual(testCase,unchanged,summary(3));
verifyEqual(testCase,height(info),0);
end

%% Helpers

function summary = makeSummary(peakPct)
% Three trials at 0 / half / full scale, for both clamp channels plus one
% photometry row that must not be touched.
trace = [0; peakPct/2; peakPct].*ones(3,10);
mk = @(name,data,stageData) struct('animal','SL431','date','20260715','session','s1',...
    'task','Random-clamp','event','Tone','name',name,'system','NI','data',data,...
    'stageAvg',struct('data',stageData),'stageMax',struct('data',stageData),...
    'stageArea',struct('data',2*stageData),'finalFs',50);
stage = [0; peakPct/2; peakPct];
summary = [mk('redClamp',trace,stage), mk('blueClamp',trace,stage), ...
           mk('NAc-clamp',0.5*ones(3,10),0.5*ones(3,1))];
end

function animals = makeAnimals(peakPct)
% animals struct: no date or session field
trace = [0; peakPct/2; peakPct].*ones(3,10);
animals = struct('animal','SL431','task','Random-clamp','event','Tone','name','redClamp',...
    'system','NI','data',trace,'stageMax',struct('data',[0;peakPct/2;peakPct]),'finalFs',50);
end
