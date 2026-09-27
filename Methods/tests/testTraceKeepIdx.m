function tests = testTraceKeepIdx
tests = functiontests(localfunctions);
end

%% getTraces

function testKeepIdxMarksDroppedEdgeEvents(testCase)
% Events whose window runs off the start or the end of the recording are
% dropped from traces, and keepIdx must say which ones.
Fs = 50; signal = 1:(60*Fs); % 60 s recording
events = [5 20 40 58]*Fs;    % first and last cannot fill a [-15,15] window
[traces,~,keepIdx] = getTraces(events,[-15,15],signal,params=struct(),signalFs=Fs,syncFs=Fs,sameSystem=true);

verifyEqual(testCase,keepIdx,[false;true;true;false]);
verifyEqual(testCase,size(traces,1),2);
% The surviving rows are the ones keepIdx marks, in order
verifyEqual(testCase,traces(1,1),(20-15)*Fs,'AbsTol',1);
verifyEqual(testCase,traces(2,1),(40-15)*Fs,'AbsTol',1);
end

function testKeepIdxMarksInteriorNaNs(testCase)
% rmmissing drops a row for any missing sample, not just at the edges.
Fs = 50; signal = ones(1,60*Fs);
signal(30*Fs) = NaN; % lands inside the second event's window
events = [20 31 45]*Fs;
[traces,~,keepIdx] = getTraces(events,[-2,2],signal,params=struct(),signalFs=Fs,syncFs=Fs,sameSystem=true);

verifyEqual(testCase,keepIdx,[true;false;true]);
verifyEqual(testCase,size(traces,1),2);
verifyFalse(testCase,any(isnan(traces),'all'));
end

function testKeepIdxAllTrueWhenNothingDropped(testCase)
Fs = 50; signal = ones(1,60*Fs);
events = [20 30 40]*Fs;
[traces,~,keepIdx] = getTraces(events,[-2,2],signal,params=struct(),signalFs=Fs,syncFs=Fs,sameSystem=true);
verifyEqual(testCase,keepIdx,true(3,1));
verifyEqual(testCase,size(traces,1),3);

% rmmissing=false keeps the NaN rows, so nothing is reported as dropped
signal(30*Fs) = NaN;
[traces,~,keepIdx] = getTraces(events,[-2,2],signal,params=struct(),signalFs=Fs,syncFs=Fs,...
                               sameSystem=true,rmmissing=false);
verifyEqual(testCase,keepIdx,true(3,1));
verifyEqual(testCase,size(traces,1),3);
end

%% plotTraces

function testPlotTracesPassesKeepIdxThrough(testCase)
Fs = 50; signal = 1:(60*Fs);
events = [5 20 40 58]*Fs;
params.sync.behaviorFs = Fs; params.sync.timeNI = (0:numel(signal)-1)/Fs;
[traces,~,keepIdx] = plotTraces(events,[-15,15],signal,params,...
                                signalFs=Fs,signalSystem='ni',eventSystem='ni',plot=false);
verifyEqual(testCase,keepIdx,[false;true;true;false]);
verifyEqual(testCase,size(traces,1),2);
end

function testPlotTracesKeepIdxWhenOnlyPlotting(testCase)
% Plot-only calls extract nothing, so every input row is kept.
fig = figure('Visible','off');
cleanup = onCleanup(@() close(fig)); %#ok<NASGU>
[traces,~,keepIdx] = plotTraces(ones(4,10),1:10);
verifyEqual(testCase,keepIdx,true(4,1));
verifyEqual(testCase,size(traces,1),4);
end

%% analyzeTraces

function testTrialInfoIsSlicedToKeptTrials(testCase)
% The end-to-end guarantee: rows of data, trialNumber and trialTable all
% describe the same trials, even when an edge trial is dropped.
Fs = 50; nSec = 60;
timeSeries = struct('name','redClamp','data',ones(1,nSec*Fs),'finalFs',Fs,'system','NI');
params.session.animal = 'SL479'; params.session.date = '20260924';
params.session.name = 'testSession'; params.session.task = 'Random-clamp';
params.session.baselineSystem = 'ni';
params.sync.behaviorFs = Fs; params.sync.timeNI = (0:nSec*Fs-1)/Fs;
events = {[5 20 40 58]*Fs};      % first and last do not fit a [-15,15] window
trialNumber = {[11 12 13 14]'};  % trial numbers need not start at 1
trialTable = table((1:20)',(1:20)','VariableNames',{'TrialNumber','performing'});

analysis = analyzeTraces(timeSeries,zeros(1,nSec*Fs),events,{'Airpuff'},params,...
    timeRange=[-15,15],stageTime=[-2,0;0,2],nboot=2,save=false,...
    trialNumber=trialNumber,trialTable=trialTable);

clampRow = analysis(strcmp({analysis.name},'redClamp'));
verifyEqual(testCase,size(clampRow.data,1),2);
verifyEqual(testCase,clampRow.trialInfo.trialNumber,[12;13]);
verifyEqual(testCase,height(clampRow.trialInfo.trialTable),2);
verifyEqual(testCase,clampRow.trialInfo.trialTable.TrialNumber,[12;13]);
verifyEqual(testCase,size(clampRow.stageAvg.data,1),2);

% Lick traces are built by getLicks, which drops nothing, so that row keeps
% all four trials and stays self-consistent.
lickRow = analysis(strcmp({analysis.name},'Lick'));
verifyEqual(testCase,size(lickRow.data.lickRate,1),4);
verifyEqual(testCase,lickRow.trialInfo.trialNumber,[11;12;13;14]);
verifyEqual(testCase,height(lickRow.trialInfo.trialTable),4);
end
