function generateBarplotSearch(plotData,options,statsMode)

    if nargin == 2; statsMode = 'none'; end

    nCat = numel(fieldnames(plotData.selectedSearchesTable));
    xLabelText = {};

    for iCat = 1:nCat
        
        catName = ['cat' num2str(iCat)];
        queryRunFeature = plotData.queryRunFeature.(catName);
        selectedSearchesTable = plotData.selectedSearchesTable.(catName);
        
        if height(selectedSearchesTable) == 0; continue; end
        
        nDataPerRun = length(selectedSearchesTable.(queryRunFeature){1});
        allDataToPlot = reshape(cell2mat(selectedSearchesTable.(queryRunFeature))',nDataPerRun,[])';        

        if options.isOverallExperimentPlot == 1; allCellNames = cell2mat(selectedSearchesTable.overallCellName);
        else; allCellNames = cell2mat(selectedSearchesTable.mouseCell); end
        
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
            
            meansCellsToPlot = reshape(cell2mat(meansCellsToPlot)',nDataPerRun,[])';
            meanCellsToPlot = nanmean(meansCellsToPlot,1);
            
            if height(meansCellsToPlot) > 1; stdCellsToPlot = nanstd(meansCellsToPlot,1);
            else; stdCellsToPlot = zeros(1,length(meansCellsToPlot)); end
            
        elseif strcmp(statsMode,'none')
            
            meanRunsToPlot = nanmean(allDataToPlot,1);
            
            if height(allDataToPlot) > 1; stdRunsToPlot = nanstd(allDataToPlot,1);
            else; stdRunsToPlot = zeros(1,nDataPerRun); end
            
        else
            
            meanRunsToPlot = mean(allDataToPlot);
            stdRunsToPlot = std(allDataToPlot);
            medianRunsToPlot = median(allDataToPlot);
            modeRunsToPlot = mode(allDataToPlot);
            prctileRunsToPlot = prctile(allDataToPlot,[25,50,75]);
            
        end
        
        % Plotting

        if options.groupByCell == 1

            scatter(iCat*ones(size(meansCellsToPlot,1),1)*[1:nDataPerRun], meansCellsToPlot, 20, options.plotColorOverlay.(catName),'filled');
            hold on

            bar(iCat*[1:nDataPerRun],meanCellsToPlot,'FaceColor',options.plotColor.(catName),'BarWidth', 0.5);
            errorbar(iCat*[1:nDataPerRun],meanCellsToPlot,stdCellsToPlot,'.','Color',[0.3 0.3 0.3],'LineWidth', 1.2, 'CapSize', 14)

        else

            if (nCat == 1 || nDataPerRun == 1)
                
                if strcmp(statsMode,'none')
                
                    scatter(iCat*ones(size(allDataToPlot,1),1)*[1:nDataPerRun],allDataToPlot, 20, options.plotColorOverlay.(catName),'filled');
                    hold on
                    bar(iCat*[1:nDataPerRun],meanRunsToPlot,'FaceColor',options.plotColor.(catName),'BarWidth', 0.5);
                    errorbar(iCat*[1:nDataPerRun],meanRunsToPlot,stdRunsToPlot,'.','Color',[0.3 0.3 0.3],'LineWidth', 1.2, 'CapSize', 14)

                else
               
                    if strcmp(statsMode,'mean')
                        
                    scatter(iCat*ones(1,size(allDataToPlot,2)),allDataToPlot, 20, options.plotColorOverlay.(catName),'filled');
                    hold on
                    bar(iCat,meanRunsToPlot,'FaceColor',options.plotColor.(catName),'BarWidth', 0.5);
                    errorbar(iCat,meanRunsToPlot,stdRunsToPlot,'.','Color',[0.3 0.3 0.3],'LineWidth', 1.2, 'CapSize', 14)

                    elseif strcmp(statsMode,'median')
                        
                    scatter(iCat*ones(1,size(allDataToPlot,2)),allDataToPlot, 20, options.plotColorOverlay.(catName),'filled');
                    hold on
                    bar(iCat,medianRunsToPlot,'FaceColor',options.plotColor.(catName),'BarWidth', 0.5);
                       
                    elseif strcmp(statsMode,'mode')
                        
                    scatter(iCat*ones(1,size(allDataToPlot,2)),allDataToPlot, 20, options.plotColorOverlay.(catName),'filled');
                    hold on
                    bar(iCat,modeRunsToPlot,'FaceColor',options.plotColor.(catName),'BarWidth', 0.5);
                        
                    elseif strcmp(statsMode,'prctile')
                        
                    bar([1,2,3],prctileRunsToPlot,'FaceColor',options.plotColor.(catName),'BarWidth', 0.5);
                              
                    end
                    
                end
                
            else
                
                positionShift = 0.2*(-1)^(iCat);
                scatter(ones(size(allDataToPlot,1),1)*[1:nDataPerRun] + positionShift,allDataToPlot, 20, options.plotColorOverlay.(catName),'filled');
                hold on
                bar([1:nDataPerRun] + positionShift,meanRunsToPlot,'FaceColor',options.plotColor.(catName),'BarWidth', 0.3);
                errorbar([1:nDataPerRun] + positionShift,meanRunsToPlot,stdRunsToPlot,'.','Color',[0.3 0.3 0.3],'LineWidth', 1.2, 'CapSize', 14)

            end
                
        end
        
        hScatter = findobj(gca, 'Type', 'scatter');
        uistack(hScatter, 'top');
        
        if ~isempty(options.yLimits)
            
            yLimMax = max(options.yLimits(end),0);
            yLimMin = min(options.yLimits(1),0);
            yRange = abs(yLimMax - yLimMin);
            ylim([yLimMin,yLimMax]);
            
        end

        if nDataPerRun == 1 && ~strcmp(statsMode,'prctile'); xticks([1:nCat]);
        elseif nDataPerRun == 1 && strcmp(statsMode,'prctile'); xticks([1:3]);
        else; xticks([1:nDataPerRun]); end
        
        if ~isempty(options.xLabelText); xLabelText{end+1} = [options.xLabelText.(catName)]; xticklabels(xLabelText); end
        if ~isempty(options.xAxisText); xlabel(options.xAxisText,'FontSize', 10); end
        
        ylabel(options.yLabelText,'FontSize', 10)
        box on;
        set(gca, 'LineWidth', 0.5, 'FontSize', 10, 'FontName', 'Arial');
        title(options.plotTitle,'FontSize',12)
    
    end

end