function tests = testLabjackAIN9
tests = functiontests(localfunctions);
end

function testOptionalAIN9AndSharedReference(testCase)
% Exercise all selections and both modes of each shared DAC independently.
for modeMask = 0:3
    for mask = 0:15
        record = logical(bitget(mask,1:4));
        [folder,cleanup,labjack,scans] = makeRecording(record,logical(bitget(modeMask,1:2))); %#ok<ASGLU>
        concatLabjack(folder,record=double(record),saveDuplicate=false);
        result = load(fullfile(folder,'data_labjack.mat'));
        rows = labjack.rawScanIdx(record);
        refs = labjack.modScanIdx(record);
        verifyEqual(testCase,result.labjack.raw,scans(rows,:));
        verifyEqual(testCase,result.labjack.modulation,scans(refs,:));
        verifyEqual(testCase,result.sync_labjack,scans(end,:));
        verifyEqual(testCase,result.labjack.mod,labjack.mod(record));
        verifyEqual(testCase,result.labjack.modFreq,labjack.modFreq(record));
        verifyEqual(testCase,result.labjack.LEDpowers,labjack.LEDpowers(record));
        verifyEqual(testCase,result.labjack.rawScanIdx,rows);
        verifyEqual(testCase,result.labjack.modScanIdx,refs);
        verifyEqual(testCase,result.labjack.channelDAC,labjack.channelDAC(record));
        verifyEqual(testCase,result.labjack.display,labjack.display(record));
        verifyEqual(testCase,result.labjack.sourceChannelIdx,find(record));
        clear cleanup
    end
end
end

function testLegacyThreeChannelRecording(testCase)
[folder,cleanup,labjack,scans] = makeRecording([true true true false]); %#ok<ASGLU>
labjack = rmfield(labjack,{'rawScanIdx','modScanIdx','numAddressesIn','syncScanIdx','channelDAC','display'});
labjack.record = [true true true];
labjack.name = labjack.name(1:3);
labjack.mod = labjack.mod(1:3);
labjack.modFreq = labjack.modFreq(1:3);
save(fullfile(folder,'Photometry','info.mat'),'labjack');
concatLabjack(folder,record=[1 1 1],saveDuplicate=false);
result = load(fullfile(folder,'data_labjack.mat'));
verifyEqual(testCase,result.labjack.raw,scans([1 2 5],:));
verifyEqual(testCase,result.labjack.modulation,scans([3 4 6],:));
verifyEqual(testCase,result.sync_labjack,scans(7,:));
end

function testSelectingOnlyAIN9DuringConcatenation(testCase)
[folder,cleanup,labjack,scans] = makeRecording([true true true true]); %#ok<ASGLU>
concatLabjack(folder,record=[0 0 0 1],followOriginal=false,saveDuplicate=false);
result = load(fullfile(folder,'data_labjack.mat'));
verifyEqual(testCase,result.labjack.raw,scans(7,:));
verifyEqual(testCase,result.labjack.modulation,scans(4,:));
verifyEqual(testCase,result.labjack.mod,true);
verifyEqual(testCase,result.labjack.modFreq,250);
verifyEqual(testCase,result.labjack.name,{'AIN9'});
end

function testShortFinalFileAndNumericFileOrder(testCase)
[folder,cleanup,~,scans] = makeRecording([true true true true]); %#ok<ASGLU>
middle = scans + 1000;
tail = scans(:,1:3) + 2000;
temp = middle(:)';
save(fullfile(folder,'Photometry','Raw_9999.mat'),'temp');
temp = tail(:)';
save(fullfile(folder,'Photometry','Raw_10000.mat'),'temp');
concatLabjack(folder,record=[1 1 1 1],saveDuplicate=false);
result = load(fullfile(folder,'data_labjack.mat'));
expected = [scans middle tail];
verifyEqual(testCase,result.labjack.raw,expected([1 2 5 7],:));
verifyEqual(testCase,result.labjack.modulation,expected([3 4 6 4],:));
verifyEqual(testCase,result.sync_labjack,expected(end,:));
verifyEqual(testCase,result.labjack.totalLen,numel(expected));
end

function testDefaultSelectionPreservesRecordedAIN9(testCase)
[folder,cleanup,~,scans] = makeRecording([false false false true]); %#ok<ASGLU>
concatLabjack(folder,saveDuplicate=false);
result = load(fullfile(folder,'data_labjack.mat'));
verifyEqual(testCase,result.labjack.raw,scans(7,:));
verifyEqual(testCase,result.labjack.modulation,scans(4,:));
verifyEqual(testCase,result.labjack.channelDAC,1);
verifyEqual(testCase,result.labjack.sourceChannelIdx,4);
end

function testAIN9CannotBeRecoveredIfNotAcquired(testCase)
[folder,cleanup] = makeRecording([true true false false]); %#ok<ASGLU>
verifyError(testCase,@() concatLabjack(folder,record=[0 0 0 1], ...
    followOriginal=false,save=false),'concatLabjack:InputNotAcquired');
end

function testGUISharedModeAndOptionalInput(testCase)
% An old three-channel preset must still open, defaulting AIN9 to unselected.
folder = tempname;
mkdir(folder);
cleanup = onCleanup(@() rmdir(folder,'s')); %#ok<NASGU>
configFile = fullfile(folder,'config.json');
fid = fopen(configFile,'w');
presets(1) = struct('animal','legacy','spikeGLX',false, ...
    'record',[false false false],'freqMod',[false true true]);
presets(2) = struct('animal','AIN9','spikeGLX',false, ...
    'record',[false false false false],'freqMod',[false false false true]);
fprintf(fid,'%s',jsonencode(presets));
fclose(fid);
% Initialize graphics before timers can interrupt library loading, then wait
% until the dialog has completed construction and entered uiwait.
warmup = figure('Visible','off');
delete(warmup);
t = timer('ExecutionMode','fixedSpacing','Period',0.2,'TimerFcn',@editDialog);
timerCleanup = onCleanup(@() delete(t)); %#ok<NASGU>
start(t);
[labjack,~,livePlot] = inputLabjackRecordingConfig(2000,configFile);
verifyNotEmpty(testCase,labjack);
verifyEqual(testCase,labjack.record,[false false false true]);
verifyEqual(testCase,labjack.mod,[true true true true]);
verifyEqual(testCase,livePlot.channelIdx,7);

    function editDialog(~,~)
        fig = findall(groot,'Type','figure','Name','LabJack recording config');
        if isempty(fig) || ~strcmp(get(fig,'WaitStatus'),'waiting'); return; end
        stop(t);
        closeCleanup = onCleanup(@() closeDialog(fig)); %#ok<NASGU>
        control = @(tag) findobj(fig,'Tag',tag);
        verifyEqual(testCase,get(control('record4'),'Value'),0);
        for c = 1:4
            verifyEqual(testCase,get(control(sprintf('freqMod%d',c)),'Enable'),'on');
            verifyEqual(testCase,get(control(sprintf('freqMod%d',c)),'Value'),1);
        end
        popup = findobj(fig,'Style','popupmenu');
        set(popup,'Value',2);
        invoke(popup);
        verifyEqual(testCase,get(control('freqMod1'),'Value'),0);
        verifyEqual(testCase,get(control('freqMod3'),'Value'),0);
        verifyEqual(testCase,get(control('freqMod2'),'Value'),1);
        verifyEqual(testCase,get(control('freqMod4'),'Value'),1);
        for c = 1:4
            partner = [3 4 1 2];
            for value = [0 1]
                set(control(sprintf('freqMod%d',c)),'Value',value);
                invoke(control(sprintf('freqMod%d',c)));
                verifyEqual(testCase,get(control(sprintf('freqMod%d',partner(c))),'Value'),value);
            end
        end
        set(control('display4'),'Value',1);
        invoke(control('display4'));
        verifyEqual(testCase,get(control('record4'),'Value'),1);
        set(control('record4'),'Value',0);
        invoke(control('record4'));
        verifyEqual(testCase,get(control('display4'),'Value'),0);
        set(control('display4'),'Value',1);
        invoke(control('display4'));
        set(control('sessionName'),'String','test-session');
        invoke(findobj(fig,'Style','pushbutton','String','OK'));
    end
end

function invoke(control)
callback = get(control,'Callback');
callback(control,[]);
end

function closeDialog(fig)
% Release uiwait even if a GUI assertion fails; leave accepted state intact.
if isgraphics(fig); uiresume(fig); end
end

function [folder,cleanup,labjack,scans] = makeRecording(record,modes)
if nargin < 2; modes = [false true]; end
folder = tempname;
mkdir(fullfile(folder,'Photometry'));
cleanup = onCleanup(@() rmdir(folder,'s'));
labjack.record = record;
labjack.nSignals = sum(record);
labjack.name = {'first','second','PMT','AIN9'};
labjack.mod = modes([1 2 1 2]);
labjack.modFreq = [200 250 200 250];
labjack.LEDpowers = [0.8 3 0.8 3];
labjack.LEDpowersMin = [0.3 0.2 0.3 0.2];
labjack.channelDAC = [0 1 0 1];
labjack.display = record;
labjack.samplerate = 20;
labjack.numAddressesIn = 7 + double(record(4));
labjack.syncScanIdx = labjack.numAddressesIn;
labjack.rawScanIdx = [1 2 5 nan];
if record(4); labjack.rawScanIdx(4) = 7; end
labjack.modScanIdx = [3 4 6 4];
if modes(1); labjack.modScanIdx(3) = 3; end
scans = (1:labjack.numAddressesIn)'*100 + (1:labjack.samplerate);
temp = reshape(scans,1,[]);
save(fullfile(folder,'Photometry','info.mat'),'labjack');
save(fullfile(folder,'Photometry','Raw_1001.mat'),'temp');
end
