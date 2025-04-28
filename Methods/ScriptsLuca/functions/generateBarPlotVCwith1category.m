function generateBarPlotVCwith1category(cat1_queryOptoParameter, queryRunFeature, queryOptoParameterValuesOverall, cat1_selectedRunGroupsAnalysisTable, cat1_selectedRunGroups, cat1_selectedCellsIDs, groupByCell, plotBaseName, savePlotPath, savePlots, nSigns)

    if nargin == 10 || nSigns == 1
        
        nSigns = 1;
        
    end
    
    for iSign = 1:nSigns
    
        nNaNs = numel(find(isnan(queryOptoParameterValuesOverall)));
        cat1_nSelectedRunGroups = numel(cat1_selectedRunGroups);

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
                cat1_cellNameData{cat1_groupsCounter} = cell2mat(cat1_selectedRunGroupsAnalysisTable.cellName(rowsRunGroup));
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

            [~,sortingVector] = sort(queryOptoParameterValues);
            cat1_allDataToPlot = cat1_queryRunFeatureData;
            cat1_allTracesToPlot = cat1_traces(sortingVector);

        end

        queryOptoParameterValuesSorted = sort(queryOptoParameterValues);
        xValues = linspace(1,numel(queryOptoParameterValuesOverall),numel(queryOptoParameterValuesOverall));

        figure;

        if groupByCell == 1

            scatter(xValues, cat1_allDataToPlot,'k','filled');
            hold on

            bar(xValues, [cat1_meanQueryRunFeatureDataToPlot] ,'FaceColor',[0 0.4470 0.7410],'BarWidth', 0.5);
            errorbar(xValues,[cat1_meanQueryRunFeatureDataToPlot],[cat1_stdQueryRunFeatureDataToPlot],'.','Color',[0.3 0.3 0.3],'LineWidth', 1.2, 'CapSize', 14)

        else

            for iQueryOptoParameter = 1:length(cat1_queryRunFeatureData)

                scatter(find(queryOptoParameterValuesSorted == queryOptoParameterValues(iQueryOptoParameter)), cat1_allDataToPlot{iQueryOptoParameter},'k','filled');
                hold on
                bar(find(queryOptoParameterValuesSorted == queryOptoParameterValues(iQueryOptoParameter)), nanmean(cat1_allDataToPlot{iQueryOptoParameter},1),'FaceColor',[0 0.4470 0.7410],'BarWidth', 0.5);
                errorbar(find(queryOptoParameterValuesSorted == queryOptoParameterValues(iQueryOptoParameter)), nanmean(cat1_allDataToPlot{iQueryOptoParameter},1),nanstd(cat1_allDataToPlot{iQueryOptoParameter},1),'.','Color',[0.3 0.3 0.3],'LineWidth', 1.2, 'CapSize', 14)

            end

        end

        hScatter = findobj(gca, 'Type', 'scatter');
        uistack(hScatter, 'top');

        xticks(xValues)
        xticklabels(sort(queryOptoParameterValues))
        xlabel(cat1_queryOptoParameter);
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
        cmap = winter(length(queryOptoParameterValuesOverall));

        if groupByCell == 1

            for iQueryOptoParameter = 1:length(cat1_queryRunFeatureData)

                plot(preprocessSignalVCv2(squeeze(cat1_meanTracesToPlot(1,iQueryOptoParameter,:))'),'Color',cmap(iQueryOptoParameter,:),'LineWidth',2);
                hold on

            end

            for iQueryOptoParameter = 1:length(cat1_queryRunFeatureData)

                cellNames = unique(cat1_cellNameData{iQueryOptoParameter});

                for iCell = 1:numel(cellNames)

                    plot(preprocessSignalVCv2(squeeze(cat1_allTracesToPlot(iCell,iQueryOptoParameter,:))'),'Color',cmap(iQueryOptoParameter,:),'LineWidth',0.2);
                    hold on

                end

            end

            %xlim([length(squeeze(cat1_meanTracesToPlot(1,1,:)))/2-1000,length(squeeze(cat1_meanTracesToPlot(1,1,:)))/2+2000])

        else

            for iQueryOptoParameter = 1:length(cat1_queryRunFeatureData)

                plot(preprocessSignalVCv2(nanmean(cell2mat(cat1_allTracesToPlot{iQueryOptoParameter}),1)),'Color',cmap(iQueryOptoParameter,:),'LineWidth',2);
                hold on

            end

            for iQueryOptoParameter = 1:length(cat1_queryRunFeatureData)

                for iTrace = 1:size(cat1_allTracesToPlot{iQueryOptoParameter},1)

                    plot(preprocessSignalVCv2(cell2mat(cat1_allTracesToPlot{iQueryOptoParameter}(iTrace))),'Color',cmap(iQueryOptoParameter,:),'LineWidth',0.2);
                    hold on

                end

            end

            %xlim([size(cell2mat(cat1_allTracesToPlot{1}(1)),2)/2-1000,size(cell2mat(cat1_allTracesToPlot{1}(1)),2)/2+2000])

        end

        xlim([4000,7000])
        xlabel('Time')
        ylabel('Current [pA]')
        box on;
        set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');
        legend('Value 1','Value 2','Value 3','Location', 'northeast', 'FontSize',12)

        if groupByCell == 1

            plotName = [plotBaseName '_' 'traces' '_' 'byCell'];

        else

            plotName = [plotBaseName '_' 'traces' '_' 'allRuns'];

        end

        title([plotName], 'Interpreter', 'none');
        plotFullPath = fullfile(savePlotPath, plotName);

        if savePlots == 1 && strcmp('heightPulsePeak',queryRunFeature)
            exportgraphics(gcf, fullfile(strcat(plotFullPath,'.pdf')), 'ContentType', 'vector');
            saveas(gcf, plotFullPath, 'fig')
        end
    
    end

end