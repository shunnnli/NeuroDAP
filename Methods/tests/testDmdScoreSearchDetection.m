function tests = testDmdScoreSearchDetection

tests = functiontests(localfunctions);

end


function testRespondingSpotScoresAboveFlatOnes(testCase)

[cm,rm,qc,opt] = localSearch(localFlatSpots(8));
cm{1}{3} = localRespondingSpot();           % one real response among eight

d = dmdScoreSearchDetection(cm,rm,qc,opt);

verifyGreaterThan(testCase,d.chargeExcessZ{1}(3),opt.minChargeZ);
verifyGreaterThan(testCase,d.chargeExcessZ{1}(3),max(d.chargeExcessZ{1}([1 2 4:8])));

end


function testHighRsSweepsAreDroppedBeforeAveraging(testCase)

[cm,rm,qc,opt] = localSearch(localFlatSpots(8));
qc{1}.Rs(3) = 90;                            % one bad patch, shared by every spot
qc{1}.Rs_headerString(3) = 85;

d = dmdScoreSearchDetection(cm,rm,qc,opt);

verifyEqual(testCase,d.nSweepsKept{1},repmat(4,8,1));

end


function testSparseDepthBorrowsTheSearchReference(testCase)

% Depth 1 has too few spots for a MAD, so it falls back to the whole search.
[cm,rm,qc,opt] = localSearch(localFlatSpots(2));
[cm2,rm2,qc2,~] = localSearch(localFlatSpots(8));
cm = [cm; cm2]; rm = [rm; rm2]; qc = [qc; qc2];

d = dmdScoreSearchDetection(cm,rm,qc,opt);

verifyTrue(testCase,all(isfinite(d.chargeExcessZ{1})));

end


function testSpotWithNoSurvivingSweepsScoresNaN(testCase)

% Review Focus 1: every sweep of one spot dropped by the Rs filter.
[cm,rm,qc,opt] = localSearch(localFlatSpots(8));
cm{1}{1} = cm{1}{1}(1,:);
rm{1}{1} = "sweep_bad";
qc{1}.Sweep(end+1) = "sweep_bad";
qc{1}.Rs(end+1) = 90;
qc{1}.Rs_headerString(end+1) = 85;

d = dmdScoreSearchDetection(cm,rm,qc,opt);

verifyEqual(testCase,d.nSweepsKept{1}(1),0);
verifyTrue(testCase,isnan(d.chargeExcessZ{1}(1)));
verifyTrue(testCase,isnan(d.filtered{1}(1).auc));

end


%% Helpers

function timeMs = localTimeMs()

timeMs = -30:0.1:100;

end


function spots = localFlatSpots(nSpots)

timeMs = localTimeMs();
spots = cell(nSpots,1);
for k = 1:nSpots
    rng(k,'twister');
    spots{k} = 0.3*randn(5,numel(timeMs));
end

end


function traces = localRespondingSpot()

timeMs = localTimeMs();
rng(99,'twister');
traces = 0.3*randn(5,numel(timeMs));
responseMask = timeMs >= 5 & timeMs <= 25;
traces(:,responseMask) = traces(:,responseMask) - 40;

end


function [cm,rm,qc,opt] = localSearch(spots)

nSweeps = 5;
sweepNames = "sweep_" + string(1:nSweeps)';

cm = {spots};
rm = {repmat({sweepNames},numel(spots),1)};
qc = {struct('Sweep',sweepNames,'Rs',repmat(12,nSweeps,1), ...
             'Rs_headerString',repmat(14,nSweeps,1))};

opt = dmdDefaultCriteria();
opt.timeRange = [-30 100];
opt.controlWindowLength = 30;
opt.baselineGuard = 1;
opt.analysisWindowLength = 30;

end
