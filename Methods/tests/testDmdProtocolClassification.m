function tests = testDmdProtocolClassification

tests = functiontests(localfunctions);

end


function testCycleNameNormalization(testCase)

cycleNames = ["randomSearch", "Random_Search.cyc", "random search", ...
              "Chrimson_InputMapping_FullField_Nov25.cyc"];

verifyEqual(testCase,isRandomSearchCycle(cycleNames),[true true true false]);

end


function testSweepLevelProtocolCells(testCase)

protocols = {struct('cycle','FullField.cyc'); ...
             struct('cycle','randomSearch.cyc'); ...
             struct('cycle','FullField.cyc')};

cycles = getSweepCycleNames(protocols,3);
verifyEqual(testCase,cycles,["FullField.cyc"; "randomSearch.cyc"; "FullField.cyc"]);
verifyEqual(testCase,isRandomSearchCycle(cycles),[false; true; false]);

end


function testMergedProtocolStruct(testCase)

protocol = struct();
protocol.cycle = {'FullField.cyc'; 'Random_Search.cyc'; 'FullField.cyc'};

cycles = getSweepCycleNames(protocol,3);
verifyEqual(testCase,cycles,["FullField.cyc"; "Random_Search.cyc"; "FullField.cyc"]);
verifyEqual(testCase,isRandomSearchCycle(cycles),[false; true; false]);

end


function testMissingProtocolPreservesAlignment(testCase)

protocols = {struct('cycle','FullField.cyc'); []; ...
             struct('cycle','randomSearch.cyc')};

cycles = getSweepCycleNames(protocols,3);
verifySize(testCase,cycles,[3 1]);
verifyEqual(testCase,cycles(2),"");
verifyEqual(testCase,isRandomSearchCycle(cycles),[false; false; true]);

end


function testMixedEpochSelectsOnlySearchSweeps(testCase)

mixedProtocol = struct();
mixedProtocol.cycle = {'FullField.cyc'; 'randomSearch.cyc'; 'FullField.cyc'};
fullFieldProtocol = struct();
fullFieldProtocol.cycle = {'FullField.cyc'; 'FullField.cyc'};

protocolRows = {mixedProtocol; fullFieldProtocol};
sweepNames = {{'AD0_1'; 'AD0_2'; 'AD0_3'}; {'AD0_4'; 'AD0_5'}};
dmdMasks = cell(2,1);
for row = 1:2
    cycles = getSweepCycleNames(protocolRows{row},numel(sweepNames{row}));
    dmdMasks{row} = isRandomSearchCycle(cycles);
end

verifyEqual(testCase,cellfun(@any,dmdMasks),[true; false]);
verifyEqual(testCase,sweepNames{1}(dmdMasks{1}),{'AD0_2'});

end
