function tests = testDmdRefreshSearchDetection
tests = functiontests(localfunctions);
end

function testOldMasksAndScoresAreReplacedWithCurrentFallback(testCase)
c = localCell([20;31;35;36;40]);
c = dmdRefreshSearchDetection(c,1);
r = c.("Response map"){1};
verifyEqual(testCase,r.nSweepsKept{1}{1},repmat(3,8,1));
verifyEqual(testCase,r.keptRows{1}{1}{1},[true;true;true;false;false]);
verifyTrue(testCase,all(isfinite(r.chargeExcessZ{1}{1})));
verifyTrue(testCase,isfield(r,'repMap'));
end

function testZeroSurvivorsRemainExcludedButTracesArePreserved(testCase)
c = localCell([31;32;35;36;40]);
before = c.("Response map"){1}.currentMap;
c = dmdRefreshSearchDetection(c,1);
r = c.("Response map"){1};
verifyEqual(testCase,r.nSweepsKept{1}{1},zeros(8,1));
verifyTrue(testCase,all(isnan(r.chargeExcessZ{1}{1})));
verifyEqual(testCase,r.currentMap,before);
end

function testStoredSweepIdentitySupportsUnequalRepetitions(testCase)
c = localCell([20;31;35;36;40]);
r = c.("Response map"){1};
r.currentMap{1}{1}{1} = r.currentMap{1}{1}{1}([2 5],:);
names = "sweep_"+string((1:5)');
rm = repmat({names},8,1); rm{1}=names([2 5]);
r.repMap = {{rm}}; c.("Response map"){1}=r;
c = dmdRefreshSearchDetection(c,1);
verifyEqual(testCase,c.("Response map"){1}.nSweepsKept{1}{1}(1),0);
end

function testLegacyUnequalRepetitionsUseSourceRows(testCase)
c = localCell([20;31;35;36;40]);
r = c.("Response map"){1};
r.currentMap{1}{1}{1} = r.currentMap{1}{1}{1}([2 5],:);
sequence=[]; names=strings(0,1);
for rep=1:5
    for spot=1:8
        if spot==1 && ~ismember(rep,[2 5]); continue; end
        sequence(end+1,1)=spot;
        names(end+1,1)="sweep_"+rep;
    end
end
r.spotSequence={{sequence}}; r.depths={4}; c.("Response map"){1}=r;
root=tempname; mkdir(fullfile(root,'cell1'));
cleanup=onCleanup(@() rmdir(root,'s')); %#ok<NASGU>
spotsAtDepth=table(names,'VariableNames',{'Sweep'});
save(fullfile(root,'cell1','search_depth4.mat'),'spotsAtDepth');
c.Session=string(root); c.Cell=1; c.Epochs={{'search'}};
c = dmdRefreshSearchDetection(c,1);
verifyEqual(testCase,c.("Response map"){1}.nSweepsKept{1}{1}(1),0);
end

function c = localCell(rs)
rng(3); spots = arrayfun(@(~) randn(5,1301),1:8,'UniformOutput',false)';
opt = dmdDefaultCriteria(); opt.minSweeps=2; opt.fallbackMaxRs=30;
opt.outputFs=10000; opt.timeRange=[-30 100];
opt.controlWindowLength=30; opt.baselineGuard=1; opt.analysisWindowLength=30;
r.currentMap = {{spots}};
r.keptRows = {{repmat({[true;false;false;false;false]},8,1)}};
r.nSweepsKept = {{ones(8,1)}};
r.chargeExcessZ = {{nan(8,1)}};
qc=struct('Sweep',"sweep_"+string((1:5)'), 'Rs',rs,'Rs_headerString',rs);
c=table; c.Options={opt}; c.("Response map")={r};
c.QC={{{qc}}}; c.Stats={struct('filtered',{{}})};
end
