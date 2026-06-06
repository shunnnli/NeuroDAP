function generateBarPlotAPpulse(featuresStruct, queryOptoParameter, queryOptoParameterValues, queryRunFeature, selectedRunGroups, selectedCellsIDs, overallExperimentName, savePlotPath, savePlots)

    queryVariable = [queryRunFeature 'allData'];
    allData = featuresStruct.(queryVariable);
    struct2vars(featuresStruct)

    meansAllDataCells = nan(numel(selectedCellsIDs),numel(selectedRunGroups));
    stdsAllDataCells = nan(numel(selectedCellsIDs),numel(selectedRunGroups));
    figure;

    for iRunGroup = 1:numel(selectedRunGroups)

        cellNames = unique(cellNameData{1,iRunGroup});

        for iCell = 1:numel(cellNames)

            groupCells = vertcat(cellNameData{1,iRunGroup});
            groupCellRows = find(groupCells == cellNames(iCell));

            allDataCells = vertcat(allData{1,iRunGroup}{:});
            allDataCell = allDataCells(groupCellRows);

            positionOptoParameterValue = find(sort(queryOptoParameterValues) == queryOptoParameterValues(iRunGroup));

            meansAllDataCells(find(selectedCellsIDs == cellNames(iCell)),positionOptoParameterValue) = nanmean(allDataCell);
            stdsAllDataCells(find(selectedCellsIDs == cellNames(iCell)),positionOptoParameterValue) = nanstd(allDataCell);

            hscatter = scatter(positionOptoParameterValue, nanmean(allDataCell),'k','filled');
            hold on;

        end

    end

    meanAllDataCells = nanmean(meansAllDataCells,1);
    stdAllDataCells = nanstd(meansAllDataCells,1);
    xValues = linspace(1,numel(queryOptoParameterValues),numel(queryOptoParameterValues));
    bar(xValues,meanAllDataCells,'FaceColor',[0 0.4470 0.7410],'BarWidth', 0.5)
    
    if numel(selectedCellsIDs) > 1
        
        errorbar(xValues,meanAllDataCells,stdAllDataCells,'.','Color',[0 0 0],'LineWidth', 1.2, 'CapSize', 14)
    
    end

    hScatter = findobj(gca, 'Type', 'scatter');
    uistack(hScatter, 'top');

    xticks(xValues)
    xticklabels([unique(queryOptoParameterValues)])
    xlabel(queryOptoParameter);
    ylabel(queryRunFeature);
    box on;
    title([queryOptoParameter,' vs ', queryRunFeature]);
    set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');

    plotName = ['IntP' '_' overallExperimentName '_' queryRunFeature,'_vs_', queryOptoParameter];
    plotFullPath = fullfile(savePlotPath, plotName);

    if savePlots == 1
        saveas(gcf, plotFullPath, 'pdf')
        saveas(gcf, plotFullPath, 'fig')
    end

end