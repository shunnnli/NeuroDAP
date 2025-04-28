function generateDotplot(plotData,options)

    nCat = numel(fieldnames(plotData.selectedRunGroupsAnalysisTable));
    legendText = {};

    for iCat = 1:nCat
        
        catName = ['cat' num2str(iCat)];
        queryRunFeature = plotData.queryRunFeature.(catName);
        selectedRunGroupsAnalysisTable = plotData.selectedRunGroupsAnalysisTable.(catName);
        
        if height(selectedRunGroupsAnalysisTable) == 0; continue; end
        
        allDataToPlot = cell2mat(selectedRunGroupsAnalysisTable.(queryRunFeature));
    
        % Plotting
        
        plot(1:numel(allDataToPlot),allDataToPlot, '-o','Color',options.plotColor.(catName),'MarkerFaceColor',options.plotColor.(catName));
        hold on;
        xlabel(options.xLabelText.(catName),'FontSize', 10)
        xlim([0,numel(allDataToPlot)+1]);
        ylim([0,max(allDataToPlot)*1.3+1]);
        
        if ~isempty(options.legendText); legendText{end+1} = [options.legendText.(catName)]; end

        xticks([]);
        xticklabels([]);
        ylabel(options.yLabelText,'FontSize', 10)
        box on;
        set(gca, 'LineWidth', 0.5, 'FontSize', 10, 'FontName', 'Arial');
        title(options.plotTitle,'FontSize',12)
        if ~isempty(options.legendText); legend(legendText,'FontSize',10,'Location','southeast'); end
    
    end

end