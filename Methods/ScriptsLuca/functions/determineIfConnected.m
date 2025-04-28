function cellsConnectionInfo = determineIfConnected(cellsConnectionInfo,cellName,channel,plotData,options)

    nCat = numel(fieldnames(plotData.selectedRunGroupsAnalysisTable));
    options.groupByCell = 0;
    
    allDataForStats = [];

    for iCat = 1:nCat
        
        catName = ['cat' num2str(iCat)];
        queryRunFeature = plotData.queryRunFeature.(catName);
        selectedRunGroupsAnalysisTable = plotData.selectedRunGroupsAnalysisTable.(catName);
        
        if height(selectedRunGroupsAnalysisTable) == 0; continue; end
        
        allDataToPlot = cell2mat(selectedRunGroupsAnalysisTable.(queryRunFeature));

        catDataForStats = allDataToPlot;    
        allDataForStats.(catName) = catDataForStats;
    
    end
    
    if isempty(allDataForStats) || numel(fields(allDataForStats)) < 2
        
        cellsConnectionInfo.(channel).(cellName){end+1} = -1;
    
    else
    
        [effectSize,pValue] = ttest(allDataForStats.cat1,allDataForStats.cat2);

        if strcmp('heightPulsePeak',queryRunFeature) && abs(mean(allDataForStats.cat1) - mean(allDataForStats.cat2)) > 25; heuristicCriterion = 1;
        elseif strcmp('areaPulse',queryRunFeature) && abs(mean(allDataForStats.cat1) - mean(allDataForStats.cat2)) > 0.5; heuristicCriterion = 1;
        else; heuristicCriterion = 0; end
        
        if pValue < 0.05 || heuristicCriterion == 1; cellsConnectionInfo.(channel).(cellName){end+1} = 1;
        else cellsConnectionInfo.(channel).(cellName){end+1} = 0; end

    end

end