function options = generateOptionsStruct(groupByCell,isOverallExperimentPlot)

if nargin == 1; isOverallExperimentPlot = 0; end

    options = [];
    options.plotColor = [];
    options.plotColorOverlay = [];
    options.plotTitle = [];
    options.xAxisText = [];
    options.xLabelText = [];
    options.yLabelText = [];
    options.legendText = [];
    options.yLimits = [];
    options.xLimits = [];
    options.xEvent = [];
    options.groupByCell =  groupByCell;
    options.isOverallExperimentPlot = isOverallExperimentPlot;

end