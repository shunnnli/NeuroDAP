function tests = testInputParamsTable
tests = functiontests(localfunctions);
end

%% Helpers

function [Prompt,Formats,DefAns] = fixture(nSessions,varargin)
% A miniature version of inputSessionParams: one of each column kind.
names = {'Session','reloadAll','rollingWindowTime','Paradigm'};

Formats(1,1).enable = 'inactive';
Formats(2,1).type = 'check';
Formats(3,1).type = 'edit';
Formats(4,1).style = 'popupmenu';
Formats(4,1).type = 'list';
Formats(4,1).items = {'random','reward pairing','punish pairing','RTPP'};

D = cell(numel(names),nSessions);
for s = 1:nSessions
    D{1,s} = sprintf('202401%02d-SL001',s);
    D{2,s} = mod(s,2)==0;
    D{3,s} = num2str(180+s);
    D{4,s} = mod(s-1,4)+1;
end
DefAns = cell2struct(D,names,1);
Prompt = repmat(names',1,2);

for k = 1:2:numel(varargin)  % override a Formats entry: fixture(n,'items',{...})
    Formats(4,1).(varargin{k}) = varargin{k+1};
end
end

%% Round trip

function testDefaultsRoundTripUnchanged(testCase)
% Building the table and converting straight back must reproduce the input.
[Prompt,Formats,DefAns] = fixture(5);
[T,meta] = inputParamsTable_build(Prompt,Formats,DefAns);
Answer = inputParamsTable_toAnswer(T,meta);
verifyEqual(testCase,Answer,DefAns(:));
end

function testColumnTypesMatchInputsdlgContract(testCase)
% edit stays char (callers run str2double/eval), check is logical, and the
% drop-down comes back as a numeric index, not its text.
[Prompt,Formats,DefAns] = fixture(3);
[T,meta] = inputParamsTable_build(Prompt,Formats,DefAns);
Answer = inputParamsTable_toAnswer(T,meta);

verifyClass(testCase,Answer(1).Session,'char');
verifyClass(testCase,Answer(1).reloadAll,'logical');
verifyClass(testCase,Answer(1).rollingWindowTime,'char');
verifyClass(testCase,Answer(1).Paradigm,'double');
verifyEqual(testCase,Answer(2).Paradigm,2);
verifyEqual(testCase,str2double(Answer(2).rollingWindowTime),182);
end

function testTableColumnClassesDriveUitableEditors(testCase)
% uitable picks the editor from the table variable type, so these classes
% are what make the check box and the drop-down appear.
[Prompt,Formats,DefAns] = fixture(4);
[T,meta] = inputParamsTable_build(Prompt,Formats,DefAns);
verifyClass(testCase,T.reloadAll,'logical');
verifyClass(testCase,T.Paradigm,'categorical');
verifyClass(testCase,T.rollingWindowTime,'string');
verifyEqual(testCase,categories(T.Paradigm),Formats(4).items(:));
verifyEqual(testCase,[meta.editable],logical([0 1 1 1]));
end

%% Editing

function testEditedDropdownReturnsNewIndex(testCase)
% The bug this guards: returning the item text would break
% taskOptions{sessionParams(s).Paradigm} in the callers.
[Prompt,Formats,DefAns] = fixture(3);
[T,meta] = inputParamsTable_build(Prompt,Formats,DefAns);
T.Paradigm(2) = 'RTPP';
Answer = inputParamsTable_toAnswer(T,meta);
verifyEqual(testCase,Answer(2).Paradigm,4);
verifyEqual(testCase,Answer(1).Paradigm,1);
end

function testEditedCheckAndEditCellsPropagate(testCase)
[Prompt,Formats,DefAns] = fixture(3);
[T,meta] = inputParamsTable_build(Prompt,Formats,DefAns);
T.reloadAll(1) = true;
T.rollingWindowTime(3) = "42";
Answer = inputParamsTable_toAnswer(T,meta);
verifyTrue(testCase,Answer(1).reloadAll);
verifyEqual(testCase,Answer(3).rollingWindowTime,'42');
end

%% Fill behaviour (the table operation behind the fill buttons)

function testFillWholeColumn(testCase)
[Prompt,Formats,DefAns] = fixture(5);
[T,meta] = inputParamsTable_build(Prompt,Formats,DefAns);
T.Paradigm((1:height(T))') = T{2,4};
Answer = inputParamsTable_toAnswer(T,meta);
verifyEqual(testCase,[Answer.Paradigm],[2 2 2 2 2]);
end

function testFillSelectedRowsLeavesOthersAlone(testCase)
[Prompt,Formats,DefAns] = fixture(5);
[T,meta] = inputParamsTable_build(Prompt,Formats,DefAns);
before = inputParamsTable_toAnswer(T,meta);
T.reloadAll([2;4]) = T{1,2};
Answer = inputParamsTable_toAnswer(T,meta);
verifyEqual(testCase,[Answer([2 4]).reloadAll],[before(1).reloadAll before(1).reloadAll]);
verifyEqual(testCase,[Answer([1 3 5]).reloadAll],[before([1 3 5]).reloadAll]);
end

%% Edge cases

function testSingleSessionWorks(testCase)
% inputsdlg took a separate "no tiling" branch at one session, so this is
% the obvious regression point.
[Prompt,Formats,DefAns] = fixture(1);
[T,meta] = inputParamsTable_build(Prompt,Formats,DefAns);
Answer = inputParamsTable_toAnswer(T,meta);
verifySize(testCase,Answer,[1 1]);
verifyEqual(testCase,Answer,DefAns(:));
end

function testManySessionsAllSurvive(testCase)
[Prompt,Formats,DefAns] = fixture(30);
[T,meta] = inputParamsTable_build(Prompt,Formats,DefAns);
Answer = inputParamsTable_toAnswer(T,meta);
verifySize(testCase,Answer,[30 1]);
verifyEqual(testCase,Answer(30).Session,'20240130-SL001');
end

function testNoneTypeColumnIsOmitted(testCase)
[Prompt,Formats,DefAns] = fixture(2);
Formats(3,1).type = 'none';
[T,meta] = inputParamsTable_build(Prompt,Formats,DefAns);
verifyFalse(testCase,ismember('rollingWindowTime',T.Properties.VariableNames));
Answer = inputParamsTable_toAnswer(T,meta);
verifyFalse(testCase,isfield(Answer,'rollingWindowTime'));
verifyTrue(testCase,isfield(Answer,'Paradigm'));
end

function testNumericEditDefaultIsCoercedToChar(testCase)
% A wrapper that forgets num2str must not break the eval/str2double callers.
[Prompt,Formats,DefAns] = fixture(2);
[DefAns.rollingWindowTime] = deal(180);
[T,meta] = inputParamsTable_build(Prompt,Formats,DefAns);
Answer = inputParamsTable_toAnswer(T,meta);
verifyClass(testCase,Answer(1).rollingWindowTime,'char');
verifyEqual(testCase,str2double(Answer(1).rollingWindowTime),180);
end

function testEvalStyleEditDefaultSurvives(testCase)
% Shun_loadSessionData.m:22 runs eval() on recordLJ.
[Prompt,Formats,DefAns] = fixture(2);
[DefAns.rollingWindowTime] = deal('[1 1 0 0]');
[T,meta] = inputParamsTable_build(Prompt,Formats,DefAns);
Answer = inputParamsTable_toAnswer(T,meta);
verifyEqual(testCase,eval(Answer(1).rollingWindowTime),[1 1 0 0]);
end

function testListDefaultMayBeItemText(testCase)
[Prompt,Formats,DefAns] = fixture(2);
[DefAns.Paradigm] = deal('punish pairing');
[T,meta] = inputParamsTable_build(Prompt,Formats,DefAns);
Answer = inputParamsTable_toAnswer(T,meta);
verifyEqual(testCase,Answer(1).Paradigm,3);
end

%% Input validation

function testDuplicateItemsAreRejected(testCase)
% Duplicates would make the text->index lookup ambiguous on output.
[Prompt,Formats,DefAns] = fixture(2,'items',{'random','random','RTPP','x'});
verifyError(testCase,@()inputParamsTable_build(Prompt,Formats,DefAns),...
    'inputParamsTable:DuplicateItems');
end

function testOutOfRangeListDefaultIsRejected(testCase)
[Prompt,Formats,DefAns] = fixture(2);
[DefAns.Paradigm] = deal(9);
verifyError(testCase,@()inputParamsTable_build(Prompt,Formats,DefAns),...
    'inputParamsTable:IndexOutOfRange');
end

function testUnsupportedFormatTypeIsRejected(testCase)
[Prompt,Formats,DefAns] = fixture(2);
Formats(3,1).type = 'slider';
verifyError(testCase,@()inputParamsTable_build(Prompt,Formats,DefAns),...
    'inputParamsTable:UnsupportedFormat');
end

function testPromptFormatsLengthMismatchIsRejected(testCase)
[Prompt,Formats,DefAns] = fixture(2);
verifyError(testCase,@()inputParamsTable_build(Prompt(1:3,:),Formats,DefAns),...
    'inputParamsTable:SizeMismatch');
end

%% Real wrapper shapes

function testInputSessionParamsFormatsShape(testCase)
% Exercise the actual 15-column Formats/Prompt that inputSessionParams builds.
names = {'Session','Paradigm','redStim','Pavlovian','ReactionTime','minLicks',...
         'OptoTriggered','OptoInverted','RedPulseFreq','RedPulseDuration',...
         'RedStimDuration','BluePulseFreq','BluePulseDuration','BlueStimDuration',...
         'IncludeOtherStim'};

Formats(1,1).enable = 'inactive';
Formats(2,1).style = 'popupmenu';
Formats(2,1).type = 'list';
Formats(2,1).items = {'random','reward pairing','punish pairing','RTPP'};
for k = [3 4 7 8 15]; Formats(k,1).type = 'check'; end
for k = [5 6 9 10 11 12 13 14]; Formats(k,1).type = 'edit'; end

n = 12;
D = cell(numel(names),n);
for s = 1:n
    D{1,s} = sprintf('sess%02d',s);
    D{2,s} = 2;
    for k = [3 4 7 8 15]; D{k,s} = true; end
    for k = [5 6 9 10 11 12 13 14]; D{k,s} = num2str(k*10); end
end
DefAns = cell2struct(D,names,1);

[T,meta] = inputParamsTable_build(repmat(names',1,2),Formats,DefAns);
verifySize(testCase,T,[12 15]);
verifyEqual(testCase,numel(meta),15);

Answer = inputParamsTable_toAnswer(T,meta);
verifySize(testCase,Answer,[12 1]);

% The exact expressions Shun_loadSessionData.m runs on the result.
taskOptions = {'random','reward pairing','punish pairing','RTPP'};
verifyEqual(testCase,taskOptions{Answer(1).Paradigm},'reward pairing');
verifyEqual(testCase,str2double(Answer(1).ReactionTime),50);

% Edit fields are char, not string. Shun_loadSessionData.m used to gate the
% str2double on isstring(), which is false for char, so ReactionTime and
% minLicks silently stayed text. Keep that trap documented.
verifyFalse(testCase,isstring(Answer(1).ReactionTime));
verifyFalse(testCase,isstring(Answer(1).minLicks));
verifyEqual(testCase,str2double(Answer(1).minLicks),60);
end
