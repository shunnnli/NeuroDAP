function tests = testLoadSessionsLabjack
tests = functiontests(localfunctions);
end

function testDefaultPhotometrySelection(testCase)
% Exercise loadSessions itself with cached photometry, as on subsequent loads.
% Include both the new logical masks and older numeric three-channel masks.
selections = {logical([1 1 0 1]), logical([1 1 0 0]), [0 1 0], [1 1 0]};
expectMismatch = [true false true false];
for k = 1:numel(selections)
    folder = tempname;
    sessionName = '20260921-test_g0';
    sessionPath = fullfile(folder,sessionName);
    mkdir(fullfile(sessionPath,'Photometry'));
    cleanup = onCleanup(@() removeSession(folder)); %#ok<NASGU>
    record = selections{k};
    selected = logical(record);
    names = {'first','second','third','AIN9'};
    labjack.record = record;
    labjack.nSignals = sum(record);
    labjack.name = names(selected);
    labjack.mod = false(1,sum(record));
    labjack.modFreq = zeros(1,sum(record));
    labjack.samplerate = 100;
    labjack.options = struct();
    t = (0:399)/labjack.samplerate;
    allRaw = (1:4)' + sin(2*pi*(1:4)'*t);
    labjack.raw = allRaw(selected,:);
    labjack.modulation = zeros(size(labjack.raw));
    sync_labjack = double(mod(0:399,100) >= 25 & mod(0:399,100) < 50);
    save(fullfile(sessionPath,'data_labjack.mat'),'labjack','sync_labjack');

    % Keep the default recordLJ and followOriginal values. Disable only the
    % optional photometry summary plot; processing and syncing run normally.
    output = evalc('loadSessions(sessionPath,plotPhotometry=false);');
    verifyEqual(testCase,contains(output,'labjack.record differs from recordLJ'),expectMismatch(k));
    if expectMismatch(k)
        verifyTrue(testCase,contains(output,['labjack.record: ' num2str(record)]));
        verifyTrue(testCase,contains(output,'options.recordLJ: 1  1  0'));
    end
    result = load(fullfile(sessionPath,['timeseries_' sessionName '.mat']));
    verifyEqual(testCase,{result.timeSeries.name},names(selected));
    verifyEqual(testCase,numel(result.timeSeries),sum(record));
    saved = load(fullfile(sessionPath,['data_' sessionName '.mat']),'labjack');
    verifyEqual(testCase,saved.labjack.record,record);
    verifyEqual(testCase,saved.labjack.raw,allRaw(selected,:));
    clear cleanup
end
end

function removeSession(folder)
close all;
rmdir(folder,'s');
end
