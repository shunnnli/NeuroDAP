function generateLinePlot(queryOptoParameter, queryRunFeature, selectedRunGroups, selectedMembersTable, overallExperimentName, savePlotPath, savePlots)

    selectedRunGroupsAnalysisTable = selectedMembersTable(ismember(selectedMembersTable.runGroup,selectedRunGroups),:);

    queryOptoParameterValues = zeros(numel(selectedRunGroups),1);
    queryOptoParameterArray = [];
    queryRunFeatureData = {};
    cellNameData = {};

    groupsCounter = 1;

    for iRunGroup = 1:numel(selectedRunGroups)

        rowsRunGroup = find(selectedRunGroupsAnalysisTable.runGroup(:) == selectedRunGroups(iRunGroup));
        queryOptoParameterValues(groupsCounter) = selectedRunGroupsAnalysisTable.optoParameters{rowsRunGroup(1)}.(queryOptoParameter);

        queryOptoParameterArray = [queryOptoParameterArray; repmat(queryOptoParameterValues(groupsCounter),numel(rowsRunGroup),1)];
        queryRunFeatureData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTable.(queryRunFeature)(rowsRunGroup));
        cellNameData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTable.overallCellName(rowsRunGroup));

        groupsCounter = groupsCounter + 1;

    end

    allData = vertcat(queryRunFeatureData{:});
    cellNameDataArray = vertcat(cellNameData{:});

    meansQueryRunFeatureDataAllCells = nan(numel(unique(cellNameDataArray)),numel(unique(queryOptoParameterValues)));

    figure;

    for iQueryOptoParameter = 1:length(queryRunFeatureData)

        cellNames = unique(cellNameData{iQueryOptoParameter});

        for iCell = 1:numel(cellNames)

            cellRows = find(cellNameData{iQueryOptoParameter} == cellNames(iCell)); 
            meanQueryRunFeatureDataCell = nanmean(queryRunFeatureData{iQueryOptoParameter}(cellRows));
            positionOptoParameterValue = find(sort(queryOptoParameterValues) == queryOptoParameterValues(iQueryOptoParameter));
            scatter(positionOptoParameterValue, meanQueryRunFeatureDataCell,'k','filled');
            hold on;

            meansQueryRunFeatureDataAllCells(iCell,positionOptoParameterValue) = meanQueryRunFeatureDataCell;

        end

    end

    meanQueryRunFeatureDataAllCells = nanmean(meansQueryRunFeatureDataAllCells,1);
    stdQueryRunFeatureDataAllCells = nanstd(meansQueryRunFeatureDataAllCells,1);
    xValues = linspace(1,numel(queryOptoParameterValues),numel(queryOptoParameterValues));
    plot(xValues, meanQueryRunFeatureDataAllCells,'Color',[0 0.4470 0.7410],'LineWidth', 2);
    errorbar(xValues,meanQueryRunFeatureDataAllCells,stdQueryRunFeatureDataAllCells,'.','Color',[0 0 0],'LineWidth', 1.2, 'CapSize', 14)

    hScatter = findobj(gca, 'Type', 'scatter');
    uistack(hScatter, 'top');

    xticks(xValues)
    xticklabels([unique(queryOptoParameterValues)])
    xlabel(queryOptoParameter);
    ylabel(queryRunFeature);
    box on;
    title([queryOptoParameter,' vs ', queryRunFeature]);
    set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');

    plotName = ['IntP'  '_' overallExperimentName '_' queryRunFeature,'_vs_', queryOptoParameter];
    plotFullPath = fullfile(savePlotPath, plotName);
    if savePlots == 1
        saveas(gcf, plotFullPath, 'pdf')
        saveas(gcf, plotFullPath, 'fig')
    end

end