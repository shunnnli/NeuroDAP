function tests = testDmdSpotPlotState

tests = functiontests(localfunctions);

end


function testHotspotWithSurvivingSweeps(testCase)

verifyEqual(testCase,dmdSpotPlotState(true,10,8),1);

end


function testScoredButNotAHotspot(testCase)

verifyEqual(testCase,dmdSpotPlotState(false,10,8),0);

end


function testTestedButNoSweepSurvivedQC(testCase)

% The spot was stimulated and recorded, but nothing passed the quality
% limits, so there is no verdict to trust -- not "not a hotspot", unknown.
verifyEqual(testCase,dmdSpotPlotState(false,10,0),-1);

end


function testQcExclusionOutranksAStaleVerdict(testCase)

% If a verdict somehow survives alongside zero kept sweeps, the QC answer
% wins: a metric computed from no sweeps cannot be believed.
verifyEqual(testCase,dmdSpotPlotState(true,10,0),-1);

end


function testNeverStimulatedTileIsNotMarkedExcluded(testCase)

% An unsampled grid position has nothing to draw and nothing to exclude.
verifyEqual(testCase,dmdSpotPlotState(false,0,0),0);

end


function testVectorisedOverSpots(testCase)

state = dmdSpotPlotState([true false false false],[10 10 10 0],[8 8 0 0]);

verifyEqual(testCase,state,[1;0;-1;0]);

end
