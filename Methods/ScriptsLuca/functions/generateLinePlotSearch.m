function generateLinePlotSearch(plotData,options)

    nCat = numel(fieldnames(plotData.selectedSearchesTable));
    xLabelText = {};
    
    meansToPlot = nan(1,nCat);

    for iCat = 1:nCat
        
        catName = ['cat' num2str(iCat)];
        queryRunFeature = plotData.queryRunFeature.(catName);
        selectedSearchesTable = plotData.selectedSearchesTable.(catName);
        
        if height(selectedSearchesTable) == 0; continue; end
        
        nDataPerRun = length(selectedSearchesTable.(queryRunFeature){1});
        allDataToPlot = reshape(cell2mat(selectedSearchesTable.(queryRunFeature))',nDataPerRun,[])';        

        if options.isOverallExperimentPlot == 1; allCellNames = cell2mat(selectedSearchesTable.overallCellName);
        else; allCellNames = cell2mat(selectedSearchesTable.cellName); end
        
        uniqueCellNames = unique(allCellNames);
        
        if options.groupByCell == 1
            
            nCellsToPlot = numel(uniqueCellNames);
            meansCellsToPlot = {};
            
            for iCell = 1:nCellsToPlot
                
                cellName = uniqueCellNames(iCell);
                rowsCellData = find(allCellNames == cellName);
                cellDataToPlot = allDataToPlot(rowsCellData,:);                
                meanCellToPlot = nanmean(cellDataToPlot,1);
                meansCellsToPlot{end+1} = meanCellToPlot;
            
            end
            
            meansCellsToPlot = reshape(cell2mat(meansCellsToPlot)',nDataPerRun,[])';%cell2mat(meansCellsToPlot);
            meanCellsToPlot = nanmean(meansCellsToPlot,1);
            
            if height(meansCellsToPlot) > 1; stdCellsToPlot = nanstd(meansCellsToPlot,1); semCellsToPlot = nanstd(meansCellsToPlot,1)./sqrt(height(meansCellsToPlot));
            else; stdCellsToPlot = zeros(1,length(meansCellsToPlot)); semCellsToPlot = zeros(1,length(meansCellsToPlot)); end
        
        else
            
            meanRunsToPlot = nanmean(allDataToPlot,1);
            
            if height(allDataToPlot) > 1; stdRunsToPlot = nanstd(allDataToPlot,1);
            else; stdRunsToPlot = zeros(1,nDataPerRun); end
            
        end
        
        % Plotting

        if options.groupByCell == 1

%             scatter(iCat*ones(size(meansCellsToPlot,1),1)*[1:nDataPerRun], meansCellsToPlot, 20, options.plotColorOverlay.(catName),'filled');
%             hold on

%             scatter(iCat*[1:nDataPerRun],meanCellsToPlot,40,'MarkerFaceColor', options.plotColor.(catName),'MarkerEdgeColor',options.plotColor.(catName),'Marker','o');
%             errorbar(iCat*[1:nDataPerRun],meanCellsToPlot,stdCellsToPlot,'.','Color',[0.3 0.3 0.3],'LineWidth', 1.2, 'CapSize', 14)
%             hold on
            
            plot([1:nDataPerRun],meanCellsToPlot,'-','Color',options.plotColor.cat1,'LineWidth',2);
            xValues = [1:nDataPerRun];

            hold on;
            
            patch('XData',[xValues, fliplr(xValues)], 'YData',[meanCellsToPlot-semCellsToPlot, fliplr(meanCellsToPlot+semCellsToPlot)], ...
            'FaceColor',options.plotColor.(catName),'EdgeColor','none','FaceAlpha',0.25,'HandleVisibility','off');
            hold on

        else

            if nCat == 1 || nDataPerRun == 1
                
                scatter(iCat*ones(size(allDataToPlot,1),1)*[1:nDataPerRun],allDataToPlot,20,options.plotColorOverlay.(catName),'filled');
                hold on
                scatter(iCat*[1:nDataPerRun],meanRunsToPlot,40,'MarkerFaceColor',options.plotColor.(catName),'MarkerEdgeColor',options.plotColor.(catName),'Marker','o');
                errorbar(iCat*[1:nDataPerRun],meanRunsToPlot,stdRunsToPlot,'.','Color',[0.3 0.3 0.3],'LineWidth', 1.2, 'CapSize', 14)
                
            else
                
                positionShift = 0.2*(-1)^(iCat);
                scatter(ones(size(allDataToPlot,1),1)*[1:nDataPerRun] + positionShift,allDataToPlot, 20, options.plotColorOverlay.(catName),'filled');
                hold on
                scatter([1:nDataPerRun] + positionShift,meanRunsToPlot,20,'MarkerFaceColor', options.plotColor.(catName),'MarkerEdgeColor',options.plotColor.(catName));
                errorbar([1:nDataPerRun] + positionShift,meanRunsToPlot,stdRunsToPlot,'.','Color',[0.3 0.3 0.3],'LineWidth', 1.2, 'CapSize', 14)

            end
                
        end
        
        hScatter = findobj(gca, 'Type', 'scatter');
        uistack(hScatter, 'top');
        
        if ~isempty(options.xLimits)
            
            xLimMax = options.xLimits(end);
            xLimMin = options.xLimits(1);
            xlim([xLimMin,xLimMax]);
            
        end
        
        if ~isempty(options.yLimits)
            
            yLimMax = max(options.yLimits(end),0);
            yLimMin = min(options.yLimits(1),0);
            yRange = abs(yLimMax - yLimMin);
            ylim([yLimMin - 0.4*yRange,yLimMax + 0.4*yRange]);
            
        end

        if nDataPerRun == 1; xticks([1:nCat]);
        else; xticks([1:2:nDataPerRun]); xticklabels(string(round(selectedSearchesTable.radii{1}(1:2:end)))); end
        
        if ~isempty(options.xLabelText); xLabelText{end+1} = [options.xLabelText.(catName)]; xticklabels(xLabelText); end
        if ~isempty(options.xAxisText); xlabel(options.xAxisText,'FontSize', 10); end
        
        ylabel(options.yLabelText,'FontSize', 10)
        box on;
        set(gca, 'LineWidth', 0.5, 'FontSize', 10, 'FontName', 'Arial');
        title(options.plotTitle,'FontSize',12)
    
    end

end