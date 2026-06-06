function generateRatiosPlot(plotData,options)

    nCat = numel(fieldnames(plotData.selectedRunGroupsAnalysisTable));
    xLabelText = {};

    for iCat = 1:nCat
        
        catName = ['cat' num2str(iCat)];
        queryRunFeature = plotData.queryRunFeature.(catName);
        selectedRunGroupsAnalysisTable = plotData.selectedRunGroupsAnalysisTable.(catName);
        
        if height(selectedRunGroupsAnalysisTable) == 0; continue; end
        
        nDataPerRun = length(selectedRunGroupsAnalysisTable.(queryRunFeature){1});
        allDataToPlot = reshape(cell2mat(selectedRunGroupsAnalysisTable.(queryRunFeature))',nDataPerRun,[])';        
        nDataToPlot = nDataPerRun-1;
        allDataToPlot = allDataToPlot(:,2:end) ./ allDataToPlot(:,1:end-1);
        
        allCellNames = cell2mat(selectedRunGroupsAnalysisTable.cellName);
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
            
            meansCellsToPlot = cell2mat(meansCellsToPlot);
            meanCellsToPlot = nanmean(meansCellsToPlot,1);
            stdCellsToPlot = nanstd(meansCellsToPlot,1);
            
        else
            
            meanRunsToPlot = nanmean(allDataToPlot,1);
            
            if height(allDataToPlot) > 1; stdRunsToPlot = nanstd(allDataToPlot,1);
            else; stdRunsToPlot = zeros(1,nDataToPlot); end
            
        end
        
        % Plotting

        if options.groupByCell == 1

            scatter(iCat*ones(size(meansCellsToPlot,1),1)*[1:nDataToPlot], meansCellsToPlot, 20, options.plotColorOverlay.(catName),'filled');
            hold on

            bar(iCat*[1:nDataToPlot],meanCellsToPlot,'FaceColor',options.plotColor.(catName),'BarWidth', 0.5);
            errorbar(iCat*[1:nDataToPlot],meanCellsToPlot,stdCellsToPlot,'.','Color',[0.3 0.3 0.3],'LineWidth', 1.2, 'CapSize', 14)

        else

%             scatter(iCat*ones(size(allDataToPlot,1),1)*[1:nDataToPlot],allDataToPlot, 20, options.plotColorOverlay.(catName),'filled');
%             hold on
% 
%             bar(iCat*[1:nDataToPlot],meanRunsToPlot,'FaceColor',options.plotColor.(catName),'BarWidth', 0.5);
%             errorbar(iCat*[1:nDataToPlot],meanRunsToPlot,stdRunsToPlot,'.','Color',[0.3 0.3 0.3],'LineWidth', 1.2, 'CapSize', 14)

            if nCat == 1 || nDataToPlot == 1
                
                scatter(iCat*ones(size(allDataToPlot,1),1)*[1:nDataToPlot],allDataToPlot, 20, options.plotColorOverlay.(catName),'filled');
                hold on
                bar(iCat*[1:nDataToPlot],meanRunsToPlot,'FaceColor',options.plotColor.(catName),'BarWidth', 0.5);
                errorbar(iCat*[1:nDataToPlot],meanRunsToPlot,stdRunsToPlot,'.','Color',[0.3 0.3 0.3],'LineWidth', 1.2, 'CapSize', 14)
                
            else
                
                positionShift = 0.2*(-1)^(iCat);
                scatter(ones(size(allDataToPlot,1),1)*[1:nDataToPlot] + positionShift,allDataToPlot, 20, options.plotColorOverlay.(catName),'filled');
                hold on
                bar([1:nDataToPlot] + positionShift,meanRunsToPlot,'FaceColor',options.plotColor.(catName),'BarWidth', 0.3);
                errorbar([1:nDataToPlot] + positionShift,meanRunsToPlot,stdRunsToPlot,'.','Color',[0.3 0.3 0.3],'LineWidth', 1.2, 'CapSize', 14)

            end
        
        end
        
        hScatter = findobj(gca, 'Type', 'scatter');
        uistack(hScatter, 'top');
        plotLine = yline(1,'--k','LineWidth',1);
        set(plotLine, 'HandleVisibility', 'off')
        
        if ~isempty(options.yLimits)
            
            yLimMax = max(options.yLimits(end),0);
            yLimMin = min(options.yLimits(1),0);
            yRange = abs(yLimMax - yLimMin);
            ylim([yLimMin - 0.4*yRange,yLimMax + 0.4*yRange]);
            
        end

        if nDataPerRun == 1; xticks([1:nCat]);
        else; xticks([1:nDataToPlot]); xticklabels(2:nDataToPlot+1); end
        
        if ~isempty(options.xLabelText); xLabelText{end+1} = [options.xLabelText.(catName)]; xticklabels(xLabelText); end
        if ~isempty(options.xAxisText); xlabel(options.xAxisText,'FontSize', 10); end
        
        ylabel(options.yLabelText,'FontSize', 10)
        box on;
        set(gca, 'LineWidth', 0.5, 'FontSize', 10, 'FontName', 'Arial');
        title(options.plotTitle,'FontSize',12)
    
    end

end