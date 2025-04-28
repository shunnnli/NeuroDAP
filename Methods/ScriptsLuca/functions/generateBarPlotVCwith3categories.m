function generateBarPlotVCwith3categories(cat1_queryOptoParameter, cat2_queryOptoParameter, cat3_queryOptoParameter, queryRunFeature, queryOptoParameterValuesOverall, cat1_selectedRunGroupsAnalysisTable, cat2_selectedRunGroupsAnalysisTable, cat3_selectedRunGroupsAnalysisTable, cat1_selectedRunGroups, cat2_selectedRunGroups, cat3_selectedRunGroups, cat1_selectedCellsIDs, cat2_selectedCellsIDs, cat3_selectedCellsIDs, groupByCell, plotBaseName, savePlotPath, savePlots)

    nNaNs = numel(find(isnan(queryOptoParameterValuesOverall)));
    cat1_nSelectedRunGroups = numel(cat1_selectedRunGroups);
    cat2_nSelectedRunGroups = numel(cat2_selectedRunGroups);
    cat3_nSelectedRunGroups = numel(cat2_selectedRunGroups);

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
        cat1_allTracesToPlot = meansTracesAllCells;
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
    
    cat2_queryRunFeatureData = {};
    cat2_cellNameData = {};
    cat2_traces = {};
    cat2_groupsCounter = 1;

    for iRunGroup = 1:cat2_nSelectedRunGroups

        rowsRunGroup = find(cat2_selectedRunGroupsAnalysisTable.runGroup(:) == cat2_selectedRunGroups(iRunGroup));

        if ~isempty(rowsRunGroup)

            queryOptoParameterValues(cat2_groupsCounter) = cat2_selectedRunGroupsAnalysisTable.optoParameters{rowsRunGroup(1)}.(cat2_queryOptoParameter);
            cat2_queryRunFeatureData{cat2_groupsCounter} = cell2mat(cat2_selectedRunGroupsAnalysisTable.(queryRunFeature)(rowsRunGroup));
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
        cat2_allTracesToPlot = meansTracesAllCells;
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
    
    
    cat3_queryRunFeatureData = {};
    cat3_cellNameData = {};
    cat3_traces = {};
    cat3_groupsCounter = 1;

    for iRunGroup = 1:cat3_nSelectedRunGroups

        rowsRunGroup = find(cat3_selectedRunGroupsAnalysisTable.runGroup(:) == cat3_selectedRunGroups(iRunGroup));

        if ~isempty(rowsRunGroup)

            queryOptoParameterValues(cat3_groupsCounter) = cat3_selectedRunGroupsAnalysisTable.optoParameters{rowsRunGroup(1)}.(cat3_queryOptoParameter);
            cat3_queryRunFeatureData{cat3_groupsCounter} = cell2mat(cat3_selectedRunGroupsAnalysisTable.(queryRunFeature)(rowsRunGroup));
            cat3_cellNameData{cat3_groupsCounter} = cell2mat(cat3_selectedRunGroupsAnalysisTable.overallCellName(rowsRunGroup));
            cat3_traces{cat3_groupsCounter} = cat3_selectedRunGroupsAnalysisTable.trace(rowsRunGroup);
            traceDuration = size(cat3_selectedRunGroupsAnalysisTable.trace{1},2);
            
        end

        cat3_groupsCounter = cat3_groupsCounter + 1;

    end
    
    cat3_cellNameDataArray = vertcat(cat3_cellNameData{:});

    if groupByCell == 1

        meansQueryRunFeatureDataAllCells = nan(numel(unique(cat3_cellNameDataArray)),numel(queryOptoParameterValuesOverall));
        meansTracesAllCells = nan(numel(unique(cat3_cellNameDataArray)),numel(queryOptoParameterValuesOverall),traceDuration);
        counterDifferentRunGroupsSameOptoParameter = zeros(numel(unique(cat3_cellNameDataArray)),numel(queryOptoParameterValuesOverall));

        for iQueryOptoParameter = 1:length(cat3_queryRunFeatureData)

            cellNames = unique(cat3_cellNameData{iQueryOptoParameter});

            for iCell = 1:numel(cellNames)

                cellRows = find(cat3_cellNameData{iQueryOptoParameter} == cellNames(iCell)); 
                meanQueryRunFeatureDataCell = nanmean(cat3_queryRunFeatureData{iQueryOptoParameter}(cellRows));
                meanTracesCell = nanmean(cell2mat(cat3_traces{iQueryOptoParameter}(cellRows)),1);
                
                if ~isnan(queryOptoParameterValues(iQueryOptoParameter))

                    positionOptoParameterValue = find(queryOptoParameterValuesOverall == queryOptoParameterValues(iQueryOptoParameter));

                else

                    positionOptoParameterValue = find(isnan(queryOptoParameterValuesOverall));

                end
                
                if counterDifferentRunGroupsSameOptoParameter(find(cat3_selectedCellsIDs == cellNames(iCell)),positionOptoParameterValue) == 0

                    meansQueryRunFeatureDataAllCells(find(cat3_selectedCellsIDs == cellNames(iCell)),positionOptoParameterValue) = meanQueryRunFeatureDataCell;
                   meansTracesAllCells(find(cat3_selectedCellsIDs == cellNames(iCell)),positionOptoParameterValue,:) = meanTracesCell;
                    counterDifferentRunGroupsSameOptoParameter(find(cat3_selectedCellsIDs == cellNames(iCell)),positionOptoParameterValue) = 1;
                    
                else
                    
                    counterValue = counterDifferentRunGroupsSameOptoParameter(find(cat3_selectedCellsIDs == cellNames(iCell)),positionOptoParameterValue);
                    meansQueryRunFeatureDataAllCells(find(cat3_selectedCellsIDs == cellNames(iCell)),positionOptoParameterValue) = sum([counterValue.*meansQueryRunFeatureDataAllCells(find(cat3_selectedCellsIDs == cellNames(iCell)),positionOptoParameterValue); meanQueryRunFeatureDataCell],1)./(1+counterValue);
                   meansTracesAllCells(find(cat3_selectedCellsIDs == cellNames(iCell)),positionOptoParameterValue,:) = sum([counterValue.*squeeze(meansTracesAllCells(find(cat3_selectedCellsIDs == cellNames(iCell)),positionOptoParameterValue,:))'; meanTracesCell],1)./(1+counterValue);
                    counterDifferentRunGroupsSameOptoParameter(find(cat3_selectedCellsIDs == cellNames(iCell)),positionOptoParameterValue) = counterValue + 1;
                
                end
                     
            end

        end

        cat3_allDataToPlot = meansQueryRunFeatureDataAllCells;
        cat3_meanQueryRunFeatureDataToPlot = nanmean(meansQueryRunFeatureDataAllCells,1);
        cat3_stdQueryRunFeatureDataToPlot = nanstd(meansQueryRunFeatureDataAllCells,1);
        cat3_allTracesToPlot = meansTracesAllCells;
        cat3_meanTracesToPlot = nanmean(cat3_allTracesToPlot,1);
        cat3_stdTracesToPlot = nanstd(cat3_allTracesToPlot,1);
        
    else

        cat3_allDataToPlot = vertcat(cat3_queryRunFeatureData{:});
        cat3_meanQueryRunFeatureDataToPlot = nanmean(cat3_allDataToPlot,1);
        cat3_stdQueryRunFeatureDataToPlot = nanstd(cat3_allDataToPlot,1);
        cat3_allTracesToPlot = cell2mat(vertcat(cat3_traces{:}));
        cat3_meanTracesToPlot = nanmean(cat3_allTracesToPlot,1);
        cat3_stdTracesToPlot = nanstd(cat3_allTracesToPlot,1);
        
    end
    
    figure
    scatter([0], cat1_allDataToPlot,'k','filled');
    hold on
    scatter([1], cat2_allDataToPlot,'k','filled');
    hold on
    scatter([2], cat3_allDataToPlot,'k','filled');
    hold on
    
    bar([0,1,2], [cat1_meanQueryRunFeatureDataToPlot, cat2_meanQueryRunFeatureDataToPlot, cat3_meanQueryRunFeatureDataToPlot] ,'FaceColor',[0 0.4470 0.7410],'BarWidth', 0.5);
    errorbar([0,1,2],[cat1_meanQueryRunFeatureDataToPlot, cat2_meanQueryRunFeatureDataToPlot, cat3_meanQueryRunFeatureDataToPlot],[cat1_stdQueryRunFeatureDataToPlot, cat2_stdQueryRunFeatureDataToPlot, cat3_stdQueryRunFeatureDataToPlot],'.','Color',[0.3 0.3 0.3],'LineWidth', 1.2, 'CapSize', 14)

    hScatter = findobj(gca, 'Type', 'scatter');
    uistack(hScatter, 'top');

    xticks([0,1,2])
    xticklabels({'Category 1','Category 2','Category 3'})
    xlabel(['Cells']);
    ylabel(queryRunFeature);

    box on;
    set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');

    if groupByCell == 1
        
        plotName = [plotBaseName '_' 'values' '_' 'byCell'];
            
    else
    
        plotName = [plotBaseName '_' 'values' '_' 'allRuns'];
    
    end
    
    title([plotName], 'Interpreter', 'none');    
    plotFullPath = fullfile(savePlotPath, plotName);

    if savePlots == 1
        saveas(gcf, plotFullPath, 'pdf')
        saveas(gcf, plotFullPath, 'fig')
    end
    
    figure;
    
    hold on
    plot(preprocessSignalVCv2(squeeze(cat1_meanTracesToPlot)'),'Color',[0 0.4470 0.7410],'LineWidth',2);
    
%     errorUp = preprocessSignalVCv2(squeeze(cat1_meanTracesToPlot)' + squeeze(cat1_stdTracesToPlot)') ;
%     errorDown = preprocessSignalVCv2(squeeze(cat1_meanTracesToPlot)' - squeeze(cat1_stdTracesToPlot)') ;
%     xLinspace = linspace(1,size(errorUp,2),size(errorUp,2));
%     plotpatch = patch([xLinspace fliplr(xLinspace)], [errorDown fliplr(errorUp)], [0 0.4470 0.7410], 'FaceAlpha',0.5, 'EdgeColor','none');
%     plotpatch.Annotation.LegendInformation.IconDisplayStyle = 'off';
%     
    hold on
    plot(preprocessSignalVCv2(squeeze(cat2_meanTracesToPlot)'),'Color',[0.4 0.4 0.4],'LineWidth',2);
    
%     errorUp = preprocessSignalVCv2(squeeze(cat2_meanTracesToPlot)' + squeeze(cat2_stdTracesToPlot)') ;
%     errorDown = preprocessSignalVCv2(squeeze(cat2_meanTracesToPlot)' - squeeze(cat2_stdTracesToPlot)') ;
%     xLinspace = linspace(1,size(errorUp,2),size(errorUp,2));
%     plotpatch = patch([xLinspace fliplr(xLinspace)], [errorDown fliplr(errorUp)], [0.4 0.4 0.4], 'FaceAlpha',0.5, 'EdgeColor','none');
%     plotpatch.Annotation.LegendInformation.IconDisplayStyle = 'off';
%    
    hold on
    plot(preprocessSignalVCv2(squeeze(cat3_meanTracesToPlot)'),'Color',[0.9 0.1 0.1],'LineWidth',2);
    
%     errorUp = preprocessSignalVCv2(squeeze(cat3_meanTracesToPlot)' + squeeze(cat3_stdTracesToPlot)') ;
%     errorDown = preprocessSignalVCv2(squeeze(cat3_meanTracesToPlot)' - squeeze(cat3_stdTracesToPlot)') ;
%     xLinspace = linspace(1,size(errorUp,2),size(errorUp,2));
%     plotpatch = patch([xLinspace fliplr(xLinspace)], [errorDown fliplr(errorUp)], [0.9 0.1 0.1], 'FaceAlpha',0.5, 'EdgeColor','none');
%     plotpatch.Annotation.LegendInformation.IconDisplayStyle = 'off';
   
%     for iTrace = 1:size(cat1_allTracesToPlot,1)
%         
%         hold on
%         plot(preprocessSignalVCv2(squeeze(cat1_allTracesToPlot(iTrace,:))'),'Color',[0 0.4470 0.7410],'LineWidth',0.1);
%     
%     end
%     
%     for iTrace = 1:size(cat2_allTracesToPlot,1)
%         
%         hold on
%         plot(preprocessSignalVCv2(squeeze(cat2_allTracesToPlot(iTrace,:))'),'Color',[0.4 0.4 0.4],'LineWidth',0.1);
%     
%     end
% 
%     for iTrace = 1:size(cat3_allTracesToPlot,1)
%         
%         hold on
%         plot(preprocessSignalVCv2(squeeze(cat3_allTracesToPlot(iTrace,:))'),'Color',[0.4 0.4 0.4],'LineWidth',0.1);
%     
%     end
    
    xlim([length(squeeze(cat1_meanTracesToPlot))/2-1000,length(squeeze(cat1_meanTracesToPlot))/2+2000])
    xlabel('Time')
    ylabel('Current [pA]')
    box on;
    set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');
    legend('Category 1','Category 2','Category 3','Location', 'northeast', 'FontSize',12)

    if groupByCell == 1
        
        plotName = [plotBaseName '_' 'traces' '_' 'byCell'];
            
    else
    
        plotName = [plotBaseName '_' 'traces' '_' 'allRuns'];
    
    end
    
    title([plotName], 'Interpreter', 'none');
    plotFullPath = fullfile(savePlotPath, plotName);

    if savePlots == 1
        saveas(gcf, plotFullPath, 'pdf')
        saveas(gcf, plotFullPath, 'fig')
    end

end