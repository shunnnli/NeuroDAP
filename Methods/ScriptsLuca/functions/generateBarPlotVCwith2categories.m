function generateBarPlotVCwith2categories(cat1_queryOptoParameter, cat2_queryOptoParameter, queryRunFeature, queryOptoParameterValuesOverall, cat1_selectedRunGroupsAnalysisTable, cat2_selectedRunGroupsAnalysisTable, cat1_selectedRunGroups, cat2_selectedRunGroups, cat1_selectedCellsIDs, cat2_selectedCellsIDs, groupByCell, plotBaseName, savePlotPath, savePlots, nSigns)
    
    if nargin == 14 || nSigns == 1
        
        nSigns = 1;
        
    end
    
    for iSign = 1:nSigns
        
        nNaNs = numel(find(isnan(queryOptoParameterValuesOverall)));
        cat1_nSelectedRunGroups = numel(cat1_selectedRunGroups);
        cat2_nSelectedRunGroups = numel(cat2_selectedRunGroups);

        queryOptoParameterValues = nan(cat1_nSelectedRunGroups,1);

        cat1_queryRunFeatureData = {};
        cat1_cellNameData = {};
        cat1_traces = {};
        cat1_groupsCounter = 1;

        for iRunGroup = 1:cat1_nSelectedRunGroups

            rowsRunGroup = find(cat1_selectedRunGroupsAnalysisTable.runGroup(:) == cat1_selectedRunGroups(iRunGroup));

            if ~isempty(rowsRunGroup)

                queryOptoParameterValues(cat1_groupsCounter) = cat1_selectedRunGroupsAnalysisTable.optoParameters{rowsRunGroup(1)}.(cat1_queryOptoParameter);
                cat1_queryRunFeatureData{cat1_groupsCounter} = cell2mat(cat1_selectedRunGroupsAnalysisTable.(queryRunFeature)(rowsRunGroup));
                cat1_queryRunFeatureData{cat1_groupsCounter} = cat1_queryRunFeatureData{cat1_groupsCounter}(:,iSign);
                cat1_cellNameData{cat1_groupsCounter} = cell2mat(cat1_selectedRunGroupsAnalysisTable.overallCellName(rowsRunGroup));
                cat1_traces{cat1_groupsCounter} = cat1_selectedRunGroupsAnalysisTable.trace(rowsRunGroup);
                traceDuration = size(cat1_selectedRunGroupsAnalysisTable.trace{1},2);

            end

            cat1_groupsCounter = cat1_groupsCounter + 1;

        end

        cat1_cellNameDataArray = vertcat(cat1_cellNameData{:});

        if groupByCell == 1

            meansQueryRunFeatureDataAllCells = nan(numel(unique(cat1_cellNameDataArray)),numel(queryOptoParameterValuesOverall));
            meansTracesAllCells = nan(numel(unique(cat1_cellNameDataArray)),numel(queryOptoParameterValuesOverall),traceDuration);
            counterDifferentRunGroupsSameOptoParameter = zeros(numel(unique(cat1_cellNameDataArray)),numel(queryOptoParameterValuesOverall));

            for iQueryOptoParameter = 1:length(cat1_queryRunFeatureData)

                cellNames = unique(cat1_cellNameData{iQueryOptoParameter});

                for iCell = 1:numel(cellNames)

                    cellRows = find(cat1_cellNameData{iQueryOptoParameter} == cellNames(iCell)); 
                    meanQueryRunFeatureDataCell = nanmean(cat1_queryRunFeatureData{iQueryOptoParameter}(cellRows));
                    meanTracesCell = nanmean(cell2mat(cat1_traces{iQueryOptoParameter}(cellRows)),1);

                    if ~isnan(queryOptoParameterValues(iQueryOptoParameter))

                        positionOptoParameterValue = find(queryOptoParameterValuesOverall == queryOptoParameterValues(iQueryOptoParameter));

                    else

                        positionOptoParameterValue = find(isnan(queryOptoParameterValuesOverall));

                    end

                    if counterDifferentRunGroupsSameOptoParameter(find(cat1_selectedCellsIDs == cellNames(iCell)),positionOptoParameterValue) == 0

                        meansQueryRunFeatureDataAllCells(find(cat1_selectedCellsIDs == cellNames(iCell)),positionOptoParameterValue) = meanQueryRunFeatureDataCell;
                        meansTracesAllCells(find(cat1_selectedCellsIDs == cellNames(iCell)),positionOptoParameterValue,:) = meanTracesCell;
                        counterDifferentRunGroupsSameOptoParameter(find(cat1_selectedCellsIDs == cellNames(iCell)),positionOptoParameterValue) = 1;

                    else

                        counterValue = counterDifferentRunGroupsSameOptoParameter(find(cat1_selectedCellsIDs == cellNames(iCell)),positionOptoParameterValue);
                        meansQueryRunFeatureDataAllCells(find(cat1_selectedCellsIDs == cellNames(iCell)),positionOptoParameterValue) = sum([counterValue.*meansQueryRunFeatureDataAllCells(find(cat1_selectedCellsIDs == cellNames(iCell)),positionOptoParameterValue); meanQueryRunFeatureDataCell],1)./(1+counterValue);
                        meansTracesAllCells(find(cat1_selectedCellsIDs == cellNames(iCell)),positionOptoParameterValue,:) = sum([counterValue.*squeeze(meansTracesAllCells(find(cat1_selectedCellsIDs == cellNames(iCell)),positionOptoParameterValue,:))'; meanTracesCell],1)./(1+counterValue);
                        counterDifferentRunGroupsSameOptoParameter(find(cat1_selectedCellsIDs == cellNames(iCell)),positionOptoParameterValue) = counterValue + 1;

                    end

                end

            end

            cat1_allDataToPlot = meansQueryRunFeatureDataAllCells;
            cat1_meanQueryRunFeatureDataToPlot = nanmean(meansQueryRunFeatureDataAllCells,1);
            cat1_stdQueryRunFeatureDataToPlot = nanstd(meansQueryRunFeatureDataAllCells,1);
            cat1_allTracesToPlot = meansTracesAllCells - mean(meansTracesAllCells,3);
            cat1_meanTracesToPlot = nanmean(cat1_allTracesToPlot,1);
            cat1_stdTracesToPlot = nanstd(cat1_allTracesToPlot,1);

        else

            cat1_allDataToPlot = vertcat(cat1_queryRunFeatureData{:});
            cat1_meanQueryRunFeatureDataToPlot = nanmean(cat1_allDataToPlot,1);
            cat1_stdQueryRunFeatureDataToPlot = nanstd(cat1_allDataToPlot,1);
            cat1_allTracesToPlot = cell2mat(vertcat(cat1_traces{:}));
            cat1_meanTracesToPlot = nanmean(cat1_allTracesToPlot,1);
            cat1_stdTracesToPlot = nanstd(cat1_allTracesToPlot,1);

        end

        queryOptoParameterValues = nan(cat2_nSelectedRunGroups,1);

        cat2_queryRunFeatureData = {};
        cat2_cellNameData = {};
        cat2_traces = {};
        cat2_groupsCounter = 1;

        for iRunGroup = 1:cat2_nSelectedRunGroups

            rowsRunGroup = find(cat2_selectedRunGroupsAnalysisTable.runGroup(:) == cat2_selectedRunGroups(iRunGroup));

            if ~isempty(rowsRunGroup)

                queryOptoParameterValues(cat2_groupsCounter) = cat2_selectedRunGroupsAnalysisTable.optoParameters{rowsRunGroup(1)}.(cat2_queryOptoParameter);
                cat2_queryRunFeatureData{cat2_groupsCounter} = cell2mat(cat2_selectedRunGroupsAnalysisTable.(queryRunFeature)(rowsRunGroup));
                cat2_queryRunFeatureData{cat2_groupsCounter} = cat2_queryRunFeatureData{cat2_groupsCounter}(:,iSign);
                cat2_cellNameData{cat2_groupsCounter} = cell2mat(cat2_selectedRunGroupsAnalysisTable.overallCellName(rowsRunGroup));
                cat2_traces{cat2_groupsCounter} = cat2_selectedRunGroupsAnalysisTable.trace(rowsRunGroup);
                traceDuration = size(cat2_selectedRunGroupsAnalysisTable.trace{1},2);

            end

            cat2_groupsCounter = cat2_groupsCounter + 1;

        end

        cat2_cellNameDataArray = vertcat(cat2_cellNameData{:});

        if groupByCell == 1

            meansQueryRunFeatureDataAllCells = nan(numel(unique(cat2_cellNameDataArray)),numel(queryOptoParameterValuesOverall));
            meansTracesAllCells = nan(numel(unique(cat2_cellNameDataArray)),numel(queryOptoParameterValuesOverall),traceDuration);
            counterDifferentRunGroupsSameOptoParameter = zeros(numel(unique(cat2_cellNameDataArray)),numel(queryOptoParameterValuesOverall));

            for iQueryOptoParameter = 1:length(cat2_queryRunFeatureData)

                cellNames = unique(cat2_cellNameData{iQueryOptoParameter});

                for iCell = 1:numel(cellNames)

                    cellRows = find(cat2_cellNameData{iQueryOptoParameter} == cellNames(iCell)); 
                    meanQueryRunFeatureDataCell = nanmean(cat2_queryRunFeatureData{iQueryOptoParameter}(cellRows));
                    meanTracesCell = nanmean(cell2mat(cat2_traces{iQueryOptoParameter}(cellRows)),1);

                    if ~isnan(queryOptoParameterValues(iQueryOptoParameter))

                        positionOptoParameterValue = find(queryOptoParameterValuesOverall == queryOptoParameterValues(iQueryOptoParameter));

                    else

                        positionOptoParameterValue = find(isnan(queryOptoParameterValuesOverall));

                    end

                    if counterDifferentRunGroupsSameOptoParameter(find(cat2_selectedCellsIDs == cellNames(iCell)),positionOptoParameterValue) == 0

                        meansQueryRunFeatureDataAllCells(find(cat2_selectedCellsIDs == cellNames(iCell)),positionOptoParameterValue) = meanQueryRunFeatureDataCell;
                        meansTracesAllCells(find(cat2_selectedCellsIDs == cellNames(iCell)),positionOptoParameterValue,:) = meanTracesCell;
                        counterDifferentRunGroupsSameOptoParameter(find(cat2_selectedCellsIDs == cellNames(iCell)),positionOptoParameterValue) = 1;

                    else

                        counterValue = counterDifferentRunGroupsSameOptoParameter(find(cat2_selectedCellsIDs == cellNames(iCell)),positionOptoParameterValue);
                        meansQueryRunFeatureDataAllCells(find(cat2_selectedCellsIDs == cellNames(iCell)),positionOptoParameterValue) = sum([counterValue.*meansQueryRunFeatureDataAllCells(find(cat2_selectedCellsIDs == cellNames(iCell)),positionOptoParameterValue); meanQueryRunFeatureDataCell],1)./(1+counterValue);
                        meansTracesAllCells(find(cat2_selectedCellsIDs == cellNames(iCell)),positionOptoParameterValue,:) = sum([counterValue.*squeeze(meansTracesAllCells(find(cat2_selectedCellsIDs == cellNames(iCell)),positionOptoParameterValue,:))'; meanTracesCell],1)./(1+counterValue);
                        counterDifferentRunGroupsSameOptoParameter(find(cat2_selectedCellsIDs == cellNames(iCell)),positionOptoParameterValue) = counterValue + 1;

                    end

                end

            end

            cat2_allDataToPlot = meansQueryRunFeatureDataAllCells;
            cat2_meanQueryRunFeatureDataToPlot = nanmean(meansQueryRunFeatureDataAllCells,1);
            cat2_stdQueryRunFeatureDataToPlot = nanstd(meansQueryRunFeatureDataAllCells,1);
            cat2_allTracesToPlot = meansTracesAllCells - mean(meansTracesAllCells,3);
            cat2_meanTracesToPlot = nanmean(cat2_allTracesToPlot,1);
            cat2_stdTracesToPlot = nanstd(cat2_allTracesToPlot,1);

        else

            cat2_allDataToPlot = vertcat(cat2_queryRunFeatureData{:});
            cat2_meanQueryRunFeatureDataToPlot = nanmean(cat2_allDataToPlot,1);
            cat2_stdQueryRunFeatureDataToPlot = nanstd(cat2_allDataToPlot,1);
            cat2_allTracesToPlot = cell2mat(vertcat(cat2_traces{:}));
            cat2_meanTracesToPlot = nanmean(cat2_allTracesToPlot,1);
            cat2_stdTracesToPlot = nanstd(cat2_allTracesToPlot,1);

        end

        figure
        scatter([0], cat1_allDataToPlot,'k','filled');
        hold on
        scatter([1], cat2_allDataToPlot,'k','filled');
        hold on

        bar([0,1], [cat1_meanQueryRunFeatureDataToPlot, cat2_meanQueryRunFeatureDataToPlot] ,'FaceColor',[0 0.4470 0.7410],'BarWidth', 0.5);
        errorbar([0,1],[cat1_meanQueryRunFeatureDataToPlot, cat2_meanQueryRunFeatureDataToPlot],[cat1_stdQueryRunFeatureDataToPlot, cat2_stdQueryRunFeatureDataToPlot],'.','Color',[0.3 0.3 0.3],'LineWidth', 1.2, 'CapSize', 14)

        hScatter = findobj(gca, 'Type', 'scatter');
        uistack(hScatter, 'top');

        xticks([0,1])
        xticklabels({'Category 1','Category 2'})
        xlabel(['Cells']);
        ylabel(queryRunFeature);

        box on;
        set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');

        if groupByCell == 1

            plotName = [plotBaseName '_' 'values' '_' 'iSign' num2str(iSign) '_' 'byCell'];

        else

            plotName = [plotBaseName '_' 'values' '_' 'iSign' num2str(iSign) '_' 'allRuns'];

        end

        title([plotName], 'Interpreter', 'none');
        plotFullPath = fullfile(savePlotPath, plotName);

        if savePlots == 1
            exportgraphics(gcf, fullfile(strcat(plotFullPath,'.pdf')), 'ContentType', 'vector');
            saveas(gcf, plotFullPath, 'fig')
        end

        figure;

        hold on
        plot(preprocessSignalVCv2(squeeze(cat1_meanTracesToPlot)'),'Color',[0 0.4470 0.7410],'LineWidth',2);
        hold on
        plot(preprocessSignalVCv2(squeeze(cat2_meanTracesToPlot)'),'Color',[0.4 0.4 0.4],'LineWidth',2);


        for iTrace = 1:size(cat1_allTracesToPlot,1)

            hold on
            plot(preprocessSignalVCv2(squeeze(cat1_allTracesToPlot(iTrace,:))'),'Color',[0 0.4470 0.7410],'LineWidth',0.1);

        end

        for iTrace = 1:size(cat2_allTracesToPlot,1)

            hold on
            plot(preprocessSignalVCv2(squeeze(cat2_allTracesToPlot(iTrace,:))'),'Color',[0.4 0.4 0.4],'LineWidth',0.1);

        end

        if length(squeeze(cat1_meanTracesToPlot)) == 20000
            xlim([4000,7000])
        else
            xlim([length(squeeze(cat1_meanTracesToPlot))/2-1000,length(squeeze(cat1_meanTracesToPlot))/2+2000])
        end

        xlabel('Time')
        ylabel('Current [pA]')
        box on;
        set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');
        legend('Category 1','Category 2','Location', 'northeast', 'FontSize',12)

        if groupByCell == 1

            plotName = [plotBaseName '_' 'traces' '_' 'byCell'];

        else

            plotName = [plotBaseName '_' 'traces' '_' 'allRuns'];

        end

        title([plotName], 'Interpreter', 'none');
        plotFullPath = fullfile(savePlotPath, plotName);

        if savePlots == 1 && strcmp(queryRunFeature,'heightPulsePeak')
            exportgraphics(gcf, fullfile(strcat(plotFullPath,'.pdf')), 'ContentType', 'vector');
            saveas(gcf, plotFullPath, 'fig')
        end
        
    end

end