function generatePlotTracesCC(plotData,options)

    nCat = numel(fieldnames(plotData.selectedRunGroupsAnalysisTable));
    legendText = {};

    for iCat = 1:nCat
        
        catName = ['cat' num2str(iCat)];
        selectedRunGroupsAnalysisTable = plotData.selectedRunGroupsAnalysisTable.(catName);
        
        if height(selectedRunGroupsAnalysisTable) == 0; continue; end
        
        allTracesToPlot = cell2mat(selectedRunGroupsAnalysisTable.trace);
        
        if options.isOverallExperimentPlot == 1; allCellNames = cell2mat(selectedRunGroupsAnalysisTable.overallCellName);
        else; allCellNames = cell2mat(selectedRunGroupsAnalysisTable.cellName); end
        
        uniqueCellNames = unique(allCellNames);
        
        if options.groupByCell == 1
            
            nCellsToPlot = numel(uniqueCellNames);
            meansCellsToPlot = {};
            
            for iCell = 1:nCellsToPlot
                
                cellName = uniqueCellNames(iCell);
                rowsCellTraces = find(allCellNames == cellName);
                cellTracesToPlot = allTracesToPlot(rowsCellTraces,:);                
                meanCellToPlot = nanmean(cellTracesToPlot,1);
                meansCellsToPlot{end+1} = meanCellToPlot;
            
            end
            
            meansCellsToPlot = cell2mat(meansCellsToPlot.');
            meanCellsToPlot = nanmean(meansCellsToPlot,1);
            stdCellsToPlot = nanstd(meansCellsToPlot,1);
            
        else
            
            meanRunsToPlot = nanmean(allTracesToPlot,1);
            stdRunsToPlot = nanstd(allTracesToPlot,1);
            
        end
        
        if ~isempty(options.legendText) 
            legendText{end+1} = strjoin([options.legendText.(catName)," (n = ", num2str(numel(uniqueCellNames)) ")"],'');
        end
        
        % Plotting

        if options.groupByCell == 1

            plot(squeeze(meanCellsToPlot),'Color',options.plotColor.(catName),'LineWidth',1.2);
            hold on

            for iCell = 1:nCellsToPlot

                plotSingleTraces = plot(squeeze(meansCellsToPlot(iCell,:)),'Color',options.plotColor.(catName),'LineWidth',0.2);
                set(plotSingleTraces, 'HandleVisibility', 'off');
                hold on

            end
            
        else

            plot(meanRunsToPlot,'Color',options.plotColor.(catName),'LineWidth',1.2);
            hold on

            for iTrace = 1:size(allTracesToPlot,1)

                plotSingleTraces = plot(allTracesToPlot(iTrace,:),'Color',options.plotColor.(catName),'LineWidth',0.2);
                set(plotSingleTraces, 'HandleVisibility', 'off');
                hold on

            end
            
        end
        
        % Plot light

        if isempty(options.xEvent)

        plotLine = xline(5000,'--k','LineWidth',1);
        set(plotLine, 'HandleVisibility', 'off');

        else 

        plotLine = xline(options.xEvent,'--k','LineWidth',1);
        set(plotLine, 'HandleVisibility', 'off');

        end

        if isempty(options.xLimits)

        xlim([4000,7000])
        xticks([4000,5000,6000,7000]);
        xticklabels([0,0.1,0.2,0.3]);

        else

          xLimMax = options.xLimits(end);
          xLimMin = options.xLimits(1);
          xlim([xLimMin, xLimMax]);
        %               xticks([0,4000,8000,12000,16000,20000]);
        %               xticklabels([0,0.4,0.8,1.2,1.6,2.0]);
          xticks([18000,22000,26000,30000,34000,38000]);
          xticklabels([0,0.4,0.8,1.2,1.6,2.0]);
%           xticks([20000,22000,24000,26000,28000,30000]);
%           xticklabels([0,0.2,0.4,0.6,0.8,1.0]);

        end
           
        
        if ~isempty(options.yLimits)
            
            yLimMax = max(options.yLimits(end),0);
            yLimMin = min(options.yLimits(1),0);
            yRange = abs(yLimMax - yLimMin);
            ylim([yLimMin - 0.3*yRange,yLimMax + 0.3*yRange]);

        end
        
        xlabel('Time [s]','FontSize', 10)
        
        if isempty(options.yLabel); ylabel('Current [pA]','FontSize', 10);
        else; ylabel(options.yLabel,'FontSize', 10); end
            
        box on;
        set(gca, 'LineWidth', 0.5, 'FontSize', 10, 'FontName', 'Arial');
        title(options.plotTitle,'FontSize',12, 'Interpreter', 'none')
        
        yLimits = ylim;
        if abs(yLimits(2)) >= abs(yLimits(1)); legendLocation = 'northeast';
        else; legendLocation = 'southeast'; end
        
        if ~isempty(options.legendText); legend(legendText,'FontSize',10,'Location',legendLocation); end
 
    end

end