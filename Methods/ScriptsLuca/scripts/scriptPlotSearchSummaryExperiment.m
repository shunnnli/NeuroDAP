disp(['*-*-*-* Running: scriptPlotSearchSummaryExperiment *-*-*-*'])

for iDepth = 1:5
    
selectedSearchesTableDepth = selectedSearchesTable(cell2mat(selectedSearchesTable.runDepth) == iDepth, :);

initializeFig(0.9,0.9);
postTextSuptitle = ['_','Depth',num2str(iDepth)];

nRows = 4;
nColumns = 10;
masterLayout = tiledlayout(nRows,nColumns);
masterLayout.TileSpacing = 'compact';
masterLayout.Padding = 'compact';

colorBlue = [85, 161, 254]./255;
colorRed = [255, 50, 58]./255;
colorGrey = [192, 192, 192]./255;

colorBlueOverlay = [85, 161, 254]./255 * 0.6;
colorRedOverlay = [255, 50, 58]./255 * 0.6;
colorGreyOverlay = [192, 192, 192]./255 * 0.6;

groupByCell = 1;
isOverallExperimentPlot = 1;

if strcmp(responseType,'excitatory'); color = colorBlue; colorOverlay = colorBlueOverlay;
elseif strcmp(responseType,'inhibitory'); color = colorRed; colorOverlay = colorRedOverlay; end    

%% Plot hotspotCount traces
nexttile(masterLayout,1,[1 3]); axis off;

queryFeature = 'hotspotCount';
plotData = generatePlotDataStruct;
plotData.queryRunFeature.cat1 = queryFeature;
plotData.selectedSearchesTable.cat1 = selectedSearchesTableDepth;

options = generateOptionsStruct(groupByCell,isOverallExperimentPlot);
options.plotColor.cat1 = color;
options.plotColorOverlay.cat1 = colorOverlay;
options.yLabelText = queryFeature;
options.xAxisText = 'Radius [pixel]';
options.plotTitle = 'RadialDistribution - Count';

generateLinePlotSearch(plotData,options)

%% Plot hotspotCount MaxValue
nexttile(masterLayout,4,[1 1]); axis off;

queryFeature = 'hotspotCountMaxValue';
plotData = generatePlotDataStruct;
plotData.queryRunFeature.cat1 = queryFeature;
plotData.selectedSearchesTable.cat1 = selectedSearchesTableDepth;

options = generateOptionsStruct(groupByCell,isOverallExperimentPlot);
options.plotColor.cat1 = color;
options.plotColorOverlay.cat1 = colorOverlay;
options.xLabelText.cat1 = [' '];
options.yLabelText = queryFeature;
options.plotTitle = 'MaxValue';

generateBarplotSearch(plotData,options)

%% Plot hotspotCount MaxRadius
nexttile(masterLayout,5,[1 1]); axis off;

queryFeature = 'hotspotCountMaxRadius';
plotData = generatePlotDataStruct;
plotData.queryRunFeature.cat1 = queryFeature;
plotData.selectedSearchesTable.cat1 = selectedSearchesTableDepth;

options = generateOptionsStruct(groupByCell,isOverallExperimentPlot);
options.plotColor.cat1 = color;
options.plotColorOverlay.cat1 = colorOverlay;
options.xLabelText.cat1 = [' '];
options.yLabelText = queryFeature;
options.plotTitle = 'MaxRadius';
options.yLimits = [0,600];

generateBarplotSearch(plotData,options)

%% Plot hotspotCount MeanRadius
nexttile(masterLayout,6,[1 1]); axis off;

queryFeature = 'hotspotCountMeanRadius';
plotData = generatePlotDataStruct;
plotData.queryRunFeature.cat1 = queryFeature;
plotData.selectedSearchesTable.cat1 = selectedSearchesTableDepth;

options = generateOptionsStruct(groupByCell,isOverallExperimentPlot);
options.plotColor.cat1 = color;
options.plotColorOverlay.cat1 = colorOverlay;
options.xLabelText.cat1 = [' '];
options.yLabelText = queryFeature;
options.plotTitle = 'MeanRadius';
options.yLimits = [0,600];

generateBarplotSearch(plotData,options)

%% Plot hotspotCount StdRadius
nexttile(masterLayout,7,[1 1]); axis off;

queryFeature = 'hotspotCountStdRadius';
plotData = generatePlotDataStruct;
plotData.queryRunFeature.cat1 = queryFeature;
plotData.selectedSearchesTable.cat1 = selectedSearchesTableDepth;

options = generateOptionsStruct(groupByCell,isOverallExperimentPlot);
options.plotColor.cat1 = color;
options.plotColorOverlay.cat1 = colorOverlay;
options.xLabelText.cat1 = [' '];
options.yLabelText = queryFeature;
options.plotTitle = 'StdRadius';
options.yLimits = [0,600];

generateBarplotSearch(plotData,options)

%% Plot hotspotCount MedianRadius
nexttile(masterLayout,8,[1 1]); axis off;

queryFeature = 'hotspotCountMedianRadius';
plotData = generatePlotDataStruct;
plotData.queryRunFeature.cat1 = queryFeature;
plotData.selectedSearchesTable.cat1 = selectedSearchesTableDepth;

options = generateOptionsStruct(groupByCell,isOverallExperimentPlot);
options.plotColor.cat1 = color;
options.plotColorOverlay.cat1 = colorOverlay;
options.xLabelText.cat1 = [' '];
options.yLabelText = queryFeature;
options.plotTitle = 'MedianRadius';
options.yLimits = [0,600];

generateBarplotSearch(plotData,options)

%% Plot hotspotCount ModeRadius
nexttile(masterLayout,9,[1 1]); axis off;

queryFeature = 'hotspotCountModeRadius';
plotData = generatePlotDataStruct;
plotData.queryRunFeature.cat1 = queryFeature;
plotData.selectedSearchesTable.cat1 = selectedSearchesTableDepth;

options = generateOptionsStruct(groupByCell,isOverallExperimentPlot);
options.plotColor.cat1 = color;
options.plotColorOverlay.cat1 = colorOverlay;
options.xLabelText.cat1 = [' '];
options.yLabelText = queryFeature;
options.plotTitle = 'ModeRadius';
options.yLimits = [0,600];

generateBarplotSearch(plotData,options)

%% Plot hotspotCount Prctile
nexttile(masterLayout,10,[1 1]); axis off;

queryFeature = 'hotspotCountPrctile';
plotData = generatePlotDataStruct;
plotData.queryRunFeature.cat1 = queryFeature;
plotData.selectedSearchesTable.cat1 = selectedSearchesTableDepth;

options = generateOptionsStruct(groupByCell,isOverallExperimentPlot);
options.plotColor.cat1 = color;
options.plotColorOverlay.cat1 = colorOverlay;
options.yLabelText = queryFeature;
options.plotTitle = 'PrctileRadius';
options.yLimits = [0,600];

generateBarplotSearch(plotData,options)
xticklabels({'25%','50%','75%'})

%% Plot hotspotCountNorm traces
nexttile(masterLayout,1 + nColumns,[1 3]); axis off;

queryFeature = 'hotspotCountNorm';
plotData = generatePlotDataStruct;
plotData.queryRunFeature.cat1 = queryFeature;
plotData.selectedSearchesTable.cat1 = selectedSearchesTableDepth;

options = generateOptionsStruct(groupByCell,isOverallExperimentPlot);
options.plotColor.cat1 = color;
options.plotColorOverlay.cat1 = colorOverlay;
options.yLabelText = queryFeature;
options.xAxisText = 'Radius [pixel]';
options.plotTitle = 'RadialDistribution - CountNorm';

generateLinePlotSearch(plotData,options)

%% Plot hotspotCountNorm MaxValue
nexttile(masterLayout,4 + nColumns,[1 1]); axis off;

queryFeature = 'hotspotCountNormMaxValue';
plotData = generatePlotDataStruct;
plotData.queryRunFeature.cat1 = queryFeature;
plotData.selectedSearchesTable.cat1 = selectedSearchesTableDepth;

options = generateOptionsStruct(groupByCell,isOverallExperimentPlot);
options.plotColor.cat1 = color;
options.plotColorOverlay.cat1 = colorOverlay;
options.xLabelText.cat1 = [' '];
options.yLabelText = queryFeature;
options.plotTitle = 'MaxValue';

generateBarplotSearch(plotData,options)

%% Plot hotspotCountNorm MaxRadius
nexttile(masterLayout,5 + nColumns,[1 1]); axis off;

queryFeature = 'hotspotCountNormMaxRadius';
plotData = generatePlotDataStruct;
plotData.queryRunFeature.cat1 = queryFeature;
plotData.selectedSearchesTable.cat1 = selectedSearchesTableDepth;

options = generateOptionsStruct(groupByCell,isOverallExperimentPlot);
options.plotColor.cat1 = color;
options.plotColorOverlay.cat1 = colorOverlay;
options.xLabelText.cat1 = [' '];
options.yLabelText = queryFeature;
options.yLimits = [0,600];
options.plotTitle = 'MaxRadius';

generateBarplotSearch(plotData,options)

%% Plot hotspotCountNorm MeanRadius
nexttile(masterLayout,6 + nColumns,[1 1]); axis off;

queryFeature = 'hotspotCountNormMeanRadius';
plotData = generatePlotDataStruct;
plotData.queryRunFeature.cat1 = queryFeature;
plotData.selectedSearchesTable.cat1 = selectedSearchesTableDepth;

options = generateOptionsStruct(groupByCell,isOverallExperimentPlot);
options.plotColor.cat1 = color;
options.plotColorOverlay.cat1 = colorOverlay;
options.xLabelText.cat1 = [' '];
options.yLabelText = queryFeature;
options.yLimits = [0,600];
options.plotTitle = 'MeanRadius';

generateBarplotSearch(plotData,options)

%% Plot hotspotCountNorm StdRadius
nexttile(masterLayout,7 + nColumns,[1 1]); axis off;

queryFeature = 'hotspotCountNormStdRadius';
plotData = generatePlotDataStruct;
plotData.queryRunFeature.cat1 = queryFeature;
plotData.selectedSearchesTable.cat1 = selectedSearchesTableDepth;

options = generateOptionsStruct(groupByCell,isOverallExperimentPlot);
options.plotColor.cat1 = color;
options.plotColorOverlay.cat1 = colorOverlay;
options.xLabelText.cat1 = [' '];
options.yLabelText = queryFeature;
options.yLimits = [0,600];
options.plotTitle = 'StdRadius';

generateBarplotSearch(plotData,options)

%% Plot hotspotCountNorm MedianRadius
nexttile(masterLayout,8 + nColumns,[1 1]); axis off;

queryFeature = 'hotspotCountNormMedianRadius';
plotData = generatePlotDataStruct;
plotData.queryRunFeature.cat1 = queryFeature;
plotData.selectedSearchesTable.cat1 = selectedSearchesTableDepth;

options = generateOptionsStruct(groupByCell,isOverallExperimentPlot);
options.plotColor.cat1 = color;
options.plotColorOverlay.cat1 = colorOverlay;
options.xLabelText.cat1 = [' '];
options.yLabelText = queryFeature;
options.yLimits = [0,600];
options.plotTitle = 'MedianRadius';

generateBarplotSearch(plotData,options)

%% Plot hotspotCountNorm ModeRadius
nexttile(masterLayout,9 + nColumns,[1 1]); axis off;

queryFeature = 'hotspotCountNormModeRadius';
plotData = generatePlotDataStruct;
plotData.queryRunFeature.cat1 = queryFeature;
plotData.selectedSearchesTable.cat1 = selectedSearchesTableDepth;

options = generateOptionsStruct(groupByCell,isOverallExperimentPlot);
options.plotColor.cat1 = color;
options.plotColorOverlay.cat1 = colorOverlay;
options.xLabelText.cat1 = [' '];
options.yLabelText = queryFeature;
options.yLimits = [0,600];
options.plotTitle = 'ModeRadius';

generateBarplotSearch(plotData,options)

%% Plot hotspotCount Prctile
nexttile(masterLayout,10 + nColumns,[1 1]); axis off;

queryFeature = 'hotspotCountNormPrctile';
plotData = generatePlotDataStruct;
plotData.queryRunFeature.cat1 = queryFeature;
plotData.selectedSearchesTable.cat1 = selectedSearchesTableDepth;

options = generateOptionsStruct(groupByCell,isOverallExperimentPlot);
options.plotColor.cat1 = color;
options.plotColorOverlay.cat1 = colorOverlay;
options.yLabelText = queryFeature;
options.yLimits = [0,600];
options.plotTitle = 'PrctileRadius';

generateBarplotSearch(plotData,options)
xticklabels({'25%','50%','75%'})

%% Plot radialCharge traces
nexttile(masterLayout,1 + 3*nColumns,[1 3]); axis off;

queryFeature = 'radialCharge';
plotData = generatePlotDataStruct;
plotData.queryRunFeature.cat1 = queryFeature;
plotData.selectedSearchesTable.cat1 = selectedSearchesTableDepth;

options = generateOptionsStruct(groupByCell,isOverallExperimentPlot);
options.plotColor.cat1 = color;
options.plotColorOverlay.cat1 = colorOverlay;
options.yLabelText = queryFeature;
options.xAxisText = 'Radius [pixel]';
options.plotTitle = 'RadialDistribution - Charge';

generateLinePlotSearch(plotData,options)

%% Plot radialCharge MaxValue
nexttile(masterLayout,4 + 3*nColumns,[1 1]); axis off;

queryFeature = 'radialChargeMaxValue';
plotData = generatePlotDataStruct;
plotData.queryRunFeature.cat1 = queryFeature;
plotData.selectedSearchesTable.cat1 = selectedSearchesTableDepth;

options = generateOptionsStruct(groupByCell,isOverallExperimentPlot);
options.plotColor.cat1 = color;
options.plotColorOverlay.cat1 = colorOverlay;
options.xLabelText.cat1 = [' '];
options.yLabelText = queryFeature;
options.plotTitle = 'MaxValue';

generateBarplotSearch(plotData,options)

%% Plot radialCharge MaxRadius
nexttile(masterLayout,5 + 3*nColumns,[1 1]); axis off;

queryFeature = 'radialChargeMaxRadius';
plotData = generatePlotDataStruct;
plotData.queryRunFeature.cat1 = queryFeature;
plotData.selectedSearchesTable.cat1 = selectedSearchesTableDepth;

options = generateOptionsStruct(groupByCell,isOverallExperimentPlot);
options.plotColor.cat1 = color;
options.plotColorOverlay.cat1 = colorOverlay;
options.xLabelText.cat1 = [' '];
options.yLabelText = queryFeature;
options.yLimits = [0,600];
options.plotTitle = 'MaxRadius';

generateBarplotSearch(plotData,options)

%% Plot radialChargeNorm traces
nexttile(masterLayout,6 + 3*nColumns,[1 3]); axis off;

queryFeature = 'radialChargeNorm';
plotData = generatePlotDataStruct;
plotData.queryRunFeature.cat1 = queryFeature;
plotData.selectedSearchesTable.cat1 = selectedSearchesTableDepth;

options = generateOptionsStruct(groupByCell,isOverallExperimentPlot);
options.plotColor.cat1 = color;
options.plotColorOverlay.cat1 = colorOverlay;
options.yLabelText = queryFeature;
options.xAxisText = 'Radius [pixel]';
options.plotTitle = 'RadialDistribution - ChargeNorm';

generateLinePlotSearch(plotData,options)

%% Plot radialChargeNorm MaxValue
nexttile(masterLayout,9 + 3*nColumns,[1 1]); axis off;

queryFeature = 'radialChargeNormMaxValue';
plotData = generatePlotDataStruct;
plotData.queryRunFeature.cat1 = queryFeature;
plotData.selectedSearchesTable.cat1 = selectedSearchesTableDepth;

options = generateOptionsStruct(groupByCell,isOverallExperimentPlot);
options.plotColor.cat1 = color;
options.plotColorOverlay.cat1 = colorOverlay;
options.xLabelText.cat1 = [' '];
options.yLabelText = queryFeature;
options.plotTitle = 'MaxValue';

generateBarplotSearch(plotData,options)

%% Plot radialChargeNorm MaxRadius
nexttile(masterLayout,10 + 3*nColumns,[1 1]); axis off;

queryFeature = 'radialChargeNormMaxRadius';
plotData = generatePlotDataStruct;
plotData.queryRunFeature.cat1 = queryFeature;
plotData.selectedSearchesTable.cat1 = selectedSearchesTableDepth;

options = generateOptionsStruct(groupByCell,isOverallExperimentPlot);
options.plotColor.cat1 = color;
options.plotColorOverlay.cat1 = colorOverlay;
options.xLabelText.cat1 = [' '];
options.yLabelText = queryFeature;
options.yLimits = [0,600];
options.plotTitle = 'MaxRadius';

generateBarplotSearch(plotData,options)

%% Plot hotspotNumber
nexttile(masterLayout,2 + 2*nColumns,[1 1]); axis off;

queryFeature = 'hotspotNumber';
plotData = generatePlotDataStruct;
plotData.queryRunFeature.cat1 = queryFeature;
plotData.selectedSearchesTable.cat1 = selectedSearchesTableDepth;

options = generateOptionsStruct(groupByCell,isOverallExperimentPlot);
options.plotColor.cat1 = color;
options.plotColorOverlay.cat1 = colorOverlay;
options.xLabelText.cat1 = [' '];
options.yLabelText = queryFeature;
options.plotTitle = 'Number';

generateBarplotSearch(plotData,options)

%% Plot hotspotPrcSoma
nexttile(masterLayout,3 + 2*nColumns,[1 1]); axis off;

queryFeature = 'hotspotPrcSoma';
plotData = generatePlotDataStruct;
plotData.queryRunFeature.cat1 = queryFeature;
plotData.selectedSearchesTable.cat1 = selectedSearchesTableDepth;

options = generateOptionsStruct(groupByCell,isOverallExperimentPlot);
options.plotColor.cat1 = color;
options.plotColorOverlay.cat1 = colorOverlay;
options.xLabelText.cat1 = [' '];
options.yLabelText = queryFeature;
options.plotTitle = 'On soma [%]';

generateBarplotSearch(plotData,options)

%% Plot hotspotDistanceMean
nexttile(masterLayout,4 + 2*nColumns,[1 1]); axis off;

queryFeature = 'hotspotDistanceMean';
plotData = generatePlotDataStruct;
plotData.queryRunFeature.cat1 = queryFeature;
plotData.selectedSearchesTable.cat1 = selectedSearchesTableDepth;

options = generateOptionsStruct(groupByCell,isOverallExperimentPlot);
options.plotColor.cat1 = color;
options.plotColorOverlay.cat1 = colorOverlay;
options.xLabelText.cat1 = [' '];
options.yLabelText = queryFeature;
options.yLimits = [0,600];
options.plotTitle = 'DistanceMean';

generateBarplotSearch(plotData,options)

%% Plot hotspotDistanceStd
nexttile(masterLayout,5 + 2*nColumns,[1 1]); axis off;

queryFeature = 'hotspotDistanceStd';
plotData = generatePlotDataStruct;
plotData.queryRunFeature.cat1 = queryFeature;
plotData.selectedSearchesTable.cat1 = selectedSearchesTableDepth;

options = generateOptionsStruct(groupByCell,isOverallExperimentPlot);
options.plotColor.cat1 = color;
options.plotColorOverlay.cat1 = colorOverlay;
options.xLabelText.cat1 = [' '];
options.yLabelText = queryFeature;
options.yLimits = [0,600];
options.plotTitle = 'DistanceStd';

generateBarplotSearch(plotData,options)

%% Plot hotspotDistanceMedian
nexttile(masterLayout,6 + 2*nColumns,[1 1]); axis off;

queryFeature = 'hotspotDistanceMedian';
plotData = generatePlotDataStruct;
plotData.queryRunFeature.cat1 = queryFeature;
plotData.selectedSearchesTable.cat1 = selectedSearchesTableDepth;

options = generateOptionsStruct(groupByCell,isOverallExperimentPlot);
options.plotColor.cat1 = color;
options.plotColorOverlay.cat1 = colorOverlay;
options.xLabelText.cat1 = [' '];
options.yLabelText = queryFeature;
options.yLimits = [0,600];
options.plotTitle = 'DistanceMedian';

generateBarplotSearch(plotData,options)

%% Plot hotspotDistanceMode
nexttile(masterLayout,7 + 2*nColumns,[1 1]); axis off;

queryFeature = 'hotspotDistanceMode';
plotData = generatePlotDataStruct;
plotData.queryRunFeature.cat1 = queryFeature;
plotData.selectedSearchesTable.cat1 = selectedSearchesTableDepth;

options = generateOptionsStruct(groupByCell,isOverallExperimentPlot);
options.plotColor.cat1 = color;
options.plotColorOverlay.cat1 = colorOverlay;
options.xLabelText.cat1 = [' '];
options.yLabelText = queryFeature;
options.yLimits = [0,600];
options.plotTitle = 'DistanceMode';

generateBarplotSearch(plotData,options)

%% Plot hotspotDistancePrctile
nexttile(masterLayout,8 + 2*nColumns,[1 1]); axis off;

queryFeature = 'hotspotDistancePrctile';
plotData = generatePlotDataStruct;
plotData.queryRunFeature.cat1 = queryFeature;
plotData.selectedSearchesTable.cat1 = selectedSearchesTableDepth;

options = generateOptionsStruct(groupByCell,isOverallExperimentPlot);
options.plotColor.cat1 = color;
options.plotColorOverlay.cat1 = colorOverlay;
options.xLabelText.cat1 = [' '];
options.yLabelText = queryFeature;
options.yLimits = [0,600];
options.plotTitle = 'DistancePrctile';

generateBarplotSearch(plotData,options)
xticklabels({'25%','50%','75%'})

%% Title and saving

plotMainTitle = [overallExperimentName, postTextSuptitle];
sgtitle(plotMainTitle,'FontSize',14,'Interpreter', 'none')

plotName = ['summaryPlotSearch_' plotMainTitle];
savePlotPath = overallExperimentDirectory;
plotFullPath = fullfile(savePlotPath, plotName);

if savePlots == 1
    exportgraphics(gcf, fullfile(strcat(plotFullPath,'.pdf')), 'ContentType', 'vector');
    saveas(gcf, plotFullPath, 'fig')
end

end