
% responseTypeBlue = 'excitatory';
% responseTypeRed = 'inhibitory';
% 
% experimentName = 'DCN_ZI_Ephys_2_2';
% 
% user = getenv('USER'); 
% matlabDirectory = ['/Users/' user '/Desktop/EPFL/Thesis/Matlab/'];
% mainDirectory = ['/Users/' user '/Desktop/EPFL/Thesis/Matlab/Experiments1/'];
% addpath(genpath(matlabDirectory));
% experimentDirectory = [mainDirectory experimentName];

experimentDirectory = expPath; % for Paolo's computer and Shun's main script
experimentNameSplit = split(experimentDirectory,'\');
experimentName = experimentNameSplit{end};

folders = dir(fullfile(experimentDirectory, '*Results*'));
folderPath = fullfile(experimentDirectory, folders(1).name); 

searchDataFileName = ['cells_DMD_', experimentName, '.mat'];
searchDataPath = fullfile(folderPath,searchDataFileName);
savePlotPath = folderPath;

searchData = load(searchDataPath);
searchData = searchData.cells;

nSearchCells = height(searchData);

columnNames = {'mouseName','mouseCell','mouseEpoch','runDepth','holdingVoltage','responseType','runProtocol','hotspotNumber','hotspotDistances','hotspotDistanceMean','hotspotDistanceStd','hotspotDistanceMedian','hotspotDistanceMode','hotspotDistancePrctile','radii','hotspotCount','hotspotCountMaxValue','hotspotCountMaxRadius','hotspotCountMeanRadius','hotspotCountStdRadius','hotspotCountMedianRadius','hotspotCountModeRadius','hotspotCountPrctile','hotspotCountNorm','hotspotCountNormMaxValue','hotspotCountNormMaxRadius','hotspotCountNormMeanRadius','hotspotCountNormStdRadius','hotspotCountNormMedianRadius','hotspotCountNormModeRadius','hotspotCountNormPrctile','radialCharge','radialChargeMaxValue','radialChargeMaxRadius','radialChargeNorm','radialChargeNormMaxValue','radialChargeNormMaxRadius'};
dataTypes = {'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell'};
searchFeaturesTable = table('Size', [0, numel(columnNames)], 'VariableNames', columnNames, 'VariableTypes', dataTypes);
rowTable = 1; warning('off', 'MATLAB:table:RowsAddedExistingVars');

isOverallExperimentPlot = 0;
groupByCell = 0;

for iSearchCell = 1:nSearchCells
    
    searchCell = searchData.Cell(iSearchCell);
    nSearchEpochs = size(searchData.Epochs{iSearchCell,1},1);
    
    for iSearchEpoch = 1:nSearchEpochs
        
        searchEpoch = str2double(regexp(searchData.Epochs{iSearchCell,1}{iSearchEpoch,1}, '(?<=epoch)\d+', 'match', 'once'));
        
        responseMapsEpoch = searchData.(8){iSearchCell,1}.responseMap{iSearchEpoch,1};
        isResponseMapsEpoch = searchData.(8){iSearchCell,1}.isResponseMap{iSearchEpoch,1};
        hotspotMapsEpoch = searchData.(8){iSearchCell,1}.hotspotMap{iSearchEpoch,1};
        
        runProtocol = searchData.Protocol{iSearchCell,1}{iSearchEpoch,1}{1,1};
        cellX = searchData.Protocol{iSearchCell,1}{iSearchEpoch,1}{1,1}.cellX;
        cellY = searchData.Protocol{iSearchCell,1}{iSearchEpoch,1}{1,1}.cellY;
        holdingVoltage = searchData.Vhold{iSearchCell,1}(iSearchEpoch,1);
        
        if holdingVoltage <= -65; responseType = 'excitatory';
        elseif holdingVoltage >= 0; responseType = 'inhibitory';
        else; responseType = 'excitatory'; end
        
        nSearchDepth = size(responseMapsEpoch,3);
        
        initializeFig(0.9,0.9);
        nMaxDepthsPlot = 5;
        masterLayout = tiledlayout(8,2*nMaxDepthsPlot);
        masterLayout.TileSpacing = 'tight';
        masterLayout.Padding = 'tight';
        
        colorBlue = [85, 161, 254]./255;
        colorRed = [255, 50, 58]./255;
        colorGrey = [.92, .92, .92];
        colorBlueOverlay = [85, 161, 254]./255 * 0.6;
        colorRedOverlay = [255, 50, 58]./255 * 0.6;
        colorGreyOverlay = [.92, .92, .92] * 0.6;
        cmapBlue = createcolormap(colorGrey,colorBlue);
        cmapRed = createcolormap(colorRed,colorGrey);
        
        for iSearchDepth = 1:nSearchDepth

            if iSearchDepth > 5; continue; end
            
            responseMapDepth = responseMapsEpoch(:,:,iSearchDepth);
            isResponseMapDepth = isResponseMapsEpoch(:,:,iSearchDepth);
            hotspotMapDepth = hotspotMapsEpoch(:,:,iSearchDepth);
            
            nBins = 40;

            if strcmp(responseType,'excitatory')
                color = colorBlue;
                colorOverlay = colorBlueOverlay;
                cmap = cmapBlue;
                climits = [min(responseMapDepth,[],'all'),0];
                if sum(abs(climits)) == 0; climits = [-eps 0]; end
                responseMapDepth(responseMapDepth > 0) = 0;
            elseif strcmp(responseType,'inhibitory')
                color = colorRed;
                colorOverlay = colorRedOverlay;
                cmap = cmapRed;
                climits = [0,max(responseMapDepth,[],'all')];
                if sum(abs(climits)) == 0; climits = [0 eps]; end
                responseMapDepth(responseMapDepth < 0) = 0;
            end

            hotspotDistanceStats = computeHotspotDistanceStats(hotspotMapDepth,iSearchDepth,cellX,cellY);
            [radialCharge,radialChargeNorm,~] = computeRadialCharge(responseMapDepth,nBins,cellX,cellY,responseType);
            [hotspotCount,hotspotCountNorm,radii] = computeRadialCount(hotspotMapDepth,nBins,cellX,cellY);
          
            if strcmp(responseType,'excitatory'); [radialChargeMaxValue,indMaxCharge] = min(radialCharge); [radialChargeNormMaxValue,indMaxChargeNorm] = min(radialChargeNorm);
            elseif strcmp(responseType,'inhibitory'); [radialChargeMaxValue,indMaxCharge] = max(radialCharge); [radialChargeNormMaxValue,indMaxChargeNorm] = max(radialChargeNorm); end
            
            radialChargeMaxRadius = radii(indMaxCharge);
            radialChargeNormMaxRadius = radii(indMaxChargeNorm);
            
            if sum(hotspotCount) == 0
                
                hotspotCountMaxValue = 0;
                hotspotCountMaxRadius = NaN;
                hotspotCountNormMaxValue = NaN;
                hotspotCountNormMaxRadius = NaN; 
                
                radiiHotspotCountDistribution = repelem(radii, hotspotCount);   
                radiiHotspotCountNormDistribution = NaN;
                
            else
            
                [hotspotCountMaxValue,indMaxCount] = max(hotspotCount);
                [hotspotCountNormMaxValue,indMaxCountNorm] = max(hotspotCountNorm);
                hotspotCountMaxRadius = radii(indMaxCount);
                hotspotCountNormMaxRadius = radii(indMaxCountNorm);
                
                radiiHotspotCountDistribution = repelem(radii, hotspotCount);
                weights = hotspotCountNorm/sum(hotspotCountNorm); nSamples = 10000;          
                radiiHotspotCountNormDistribution = randsample(radii,nSamples,true,weights);
                
            end
                
            hotspotCountMeanRadius = mean(radiiHotspotCountDistribution);
            hotspotCountStdRadius = std(radiiHotspotCountDistribution);
            hotspotCountMedianRadius = median(radiiHotspotCountDistribution);
            hotspotCountModeRadius = mode(radiiHotspotCountDistribution);
            hotspotCountPrctile = prctile(radiiHotspotCountDistribution,[25,50,75]);

            hotspotCountNormMeanRadius = mean(radiiHotspotCountNormDistribution);
            hotspotCountNormStdRadius = std(radiiHotspotCountNormDistribution);
            hotspotCountNormMedianRadius = median(radiiHotspotCountNormDistribution);
            hotspotCountNormModeRadius = mode(radiiHotspotCountNormDistribution);
            hotspotCountNormPrctile = prctile(radiiHotspotCountNormDistribution,[25,50,75]);
            
            % Plots
            
            nexttile(masterLayout,2*(iSearchDepth-1) + 1,[1 2])
            imshow(isResponseMapDepth)
            hold on
            rectangle('Position', [0.5, 0.5, size(isResponseMapDepth,2), size(isResponseMapDepth,1)], 'EdgeColor', colorGrey, 'LineWidth', 1);
            hold on
            scatter(cellX,cellY,50,'o','filled','MarkerFaceColor','g','MarkerEdgeColor','g');
            title(['Depth ' num2str(iSearchDepth)],'FontSize',12,'Interpreter','none')
            
            nexttile(masterLayout,2*(iSearchDepth-1) + 1 + 2*nMaxDepthsPlot,[1 2])
            imshow(hotspotMapDepth)
            hold on
            rectangle('Position', [0.5, 0.5, size(hotspotMapDepth,2), size(hotspotMapDepth,1)], 'EdgeColor', colorGrey, 'LineWidth', 1);
            hold on
            scatter(cellX,cellY,50,'o','filled','MarkerFaceColor','g','MarkerEdgeColor','g');
            %title(['Hotspot map - Depth ' num2str(iSearchDepth)],'FontSize',12,'Interpreter','none')
            
            nexttile(masterLayout,2*(iSearchDepth-1) + 1 + 4*nMaxDepthsPlot,[1 2])
            imshow(responseMapDepth)
            colormap(gca, flip(cmap))
            caxis(climits); cbar = colorbar; ylabel(cbar,'Total charge [pC]','FontSize',8);
            hold on
            scatter(cellX,cellY,50,'o','filled','MarkerFaceColor','k','MarkerEdgeColor','k');
            %title(['Charge map - Depth ' num2str(iSearchDepth)],'FontSize',12,'Interpreter','none')
            
            nexttile(masterLayout,2*(iSearchDepth-1) + 1 + 6*nMaxDepthsPlot,[1 2])
            plot(radii,hotspotCount,'-o','Color',color,'MarkerFaceColor',color,'MarkerSize',4);
            xlabel('Radial distance [pixel]','FontSize',10)
            ylabel('Hotspot count [pixel]','FontSize',10)
            hold on
            xline(hotspotCountMaxRadius,'-g','LineWidth',1)
            hold on
            
            nexttile(masterLayout,2*(iSearchDepth-1) + 1 + 8*nMaxDepthsPlot,[1 2])
            plot(radii,hotspotCountNorm,'-o','Color',color,'MarkerFaceColor',color,'MarkerSize',4);
            xlabel('Radial distance [pixel]','FontSize',10)
            ylabel(["Normalized" + newline + "hotspot count [1]"],'FontSize',10)
            hold on
            xline(hotspotCountNormMaxRadius,'-g','LineWidth',1)
            hold on
            
            nexttile(masterLayout,2*(iSearchDepth-1) + 1 + 10*nMaxDepthsPlot,[1 2])
            plot(radii,radialCharge,'-o','Color',color,'MarkerFaceColor',color,'MarkerSize',4);
            xlabel('Radial distance [pixel]','FontSize',10)
            ylabel(["Normalized" + newline + "mean charge [1]"],'FontSize',10)
            hold on
            xline(radialChargeMaxRadius,'-g','LineWidth',1)
            hold on
            
            % Collect features

            searchFeaturesTable.mouseName(rowTable) = {experimentName};
            searchFeaturesTable.mouseCell(rowTable) = {searchCell};
            searchFeaturesTable.mouseEpoch(rowTable) = {searchEpoch};
            searchFeaturesTable.runDepth(rowTable) = {iSearchDepth};
            searchFeaturesTable.holdingVoltage(rowTable) = {holdingVoltage};
            searchFeaturesTable.responseType(rowTable) = {responseType};
            searchFeaturesTable.runProtocol(rowTable) = {runProtocol};
            searchFeaturesTable.hotspotNumber(rowTable) = {hotspotDistanceStats.hotspotNumber};
            searchFeaturesTable.hotspotDistances(rowTable) = {hotspotDistanceStats.hotspotDistances};
            searchFeaturesTable.hotspotDistanceMean(rowTable) = {hotspotDistanceStats.hotspotDistanceMean};
            searchFeaturesTable.hotspotDistanceStd(rowTable) = {hotspotDistanceStats.hotspotDistanceStd};
            searchFeaturesTable.hotspotDistanceMedian(rowTable) = {hotspotDistanceStats.hotspotDistanceMedian};
            searchFeaturesTable.hotspotDistanceMode(rowTable) = {hotspotDistanceStats.hotspotDistanceMode};
            searchFeaturesTable.hotspotDistancePrctile(rowTable) = {hotspotDistanceStats.hotspotDistancePrctile};
            searchFeaturesTable.radii(rowTable) = {radii};
            searchFeaturesTable.hotspotCount(rowTable) = {hotspotCount};
            searchFeaturesTable.hotspotCountMaxValue(rowTable) = {hotspotCountMaxValue};
            searchFeaturesTable.hotspotCountMaxRadius(rowTable) = {hotspotCountMaxRadius};
            searchFeaturesTable.hotspotCountMeanRadius(rowTable) = {hotspotCountMeanRadius};
            searchFeaturesTable.hotspotCountStdRadius(rowTable) = {hotspotCountStdRadius};
            searchFeaturesTable.hotspotCountMedianRadius(rowTable) = {hotspotCountMedianRadius};
            searchFeaturesTable.hotspotCountModeRadius(rowTable) = {hotspotCountModeRadius};
            searchFeaturesTable.hotspotCountPrctile(rowTable) = {hotspotCountPrctile};
            searchFeaturesTable.hotspotCountNorm(rowTable) = {hotspotCountNorm};
            searchFeaturesTable.hotspotCountNormMaxValue(rowTable) = {hotspotCountNormMaxValue};
            searchFeaturesTable.hotspotCountNormMaxRadius(rowTable) = {hotspotCountNormMaxRadius};
            searchFeaturesTable.hotspotCountNormMeanRadius(rowTable) = {hotspotCountNormMeanRadius};
            searchFeaturesTable.hotspotCountNormStdRadius(rowTable) = {hotspotCountNormStdRadius};
            searchFeaturesTable.hotspotCountNormMedianRadius(rowTable) = {hotspotCountNormMedianRadius};
            searchFeaturesTable.hotspotCountNormModeRadius(rowTable) = {hotspotCountNormModeRadius};
            searchFeaturesTable.hotspotCountNormPrctile(rowTable) = {hotspotCountNormPrctile};
            searchFeaturesTable.radialCharge(rowTable) = {radialCharge};
            searchFeaturesTable.radialChargeMaxValue(rowTable) = {radialChargeMaxValue};
            searchFeaturesTable.radialChargeMaxRadius(rowTable) = {radialChargeMaxRadius};
            searchFeaturesTable.radialChargeNorm(rowTable) = {radialChargeNorm};
            searchFeaturesTable.radialChargeNormMaxValue(rowTable) = {radialChargeNormMaxValue};
            searchFeaturesTable.radialChargeNormMaxRadius(rowTable) = {radialChargeNormMaxRadius};
            
            rowTable = rowTable + 1;
            
            % Plot stats
            
            nexttile(masterLayout,2*(iSearchDepth-1) + 1 + 12*nMaxDepthsPlot,[1 1])

            queryFeature = 'hotspotDistances';
            plotData = generatePlotDataStruct;
            plotData.queryRunFeature.cat1 = queryFeature;
            plotData.selectedSearchesTable.cat1 = searchFeaturesTable(end,:);

            options = generateOptionsStruct(groupByCell,isOverallExperimentPlot);
            options.plotColor.cat1 = color;
            options.plotColorOverlay.cat1 = colorOverlay;
            options.xLabelText.cat1 = [' '];
            options.yLabelText = ["Hotspot distance" + newline + "mean [pixel]"];
            options.yLimits = [0,600];
            
            generateBarplotSearch(plotData,options,'mean')
            
            nexttile(masterLayout,2*(iSearchDepth-1) + 2 + 12*nMaxDepthsPlot,[1 1])

            queryFeature = 'hotspotDistances';
            plotData = generatePlotDataStruct;
            plotData.queryRunFeature.cat1 = queryFeature;
            plotData.selectedSearchesTable.cat1 = searchFeaturesTable(end,:);

            options = generateOptionsStruct(groupByCell,isOverallExperimentPlot);
            options.plotColor.cat1 = color;
            options.plotColorOverlay.cat1 = colorOverlay;
            options.xLabelText.cat1 = [' '];
            options.yLabelText = ["Hotspot distance" + newline + "median [pixel]"];
            options.yLimits = [0,600];
            
            generateBarplotSearch(plotData,options,'median')
            
            nexttile(masterLayout,2*(iSearchDepth-1) + 1 + 14*nMaxDepthsPlot,[1 1])

            queryFeature = 'hotspotDistances';
            plotData = generatePlotDataStruct;
            plotData.queryRunFeature.cat1 = queryFeature;
            plotData.selectedSearchesTable.cat1 = searchFeaturesTable(end,:);

            options = generateOptionsStruct(groupByCell,isOverallExperimentPlot);
            options.plotColor.cat1 = color;
            options.plotColorOverlay.cat1 = colorOverlay;
            options.xLabelText.cat1 = [' '];
            options.yLabelText = ["Hotspot distance" + newline + "mode [pixel]"];
            options.yLimits = [0,600];
            
            generateBarplotSearch(plotData,options,'mode')
            
            nexttile(masterLayout,2*(iSearchDepth-1) + 2 + 14*nMaxDepthsPlot,[1 1])

            queryFeature = 'hotspotDistances';
            plotData = generatePlotDataStruct;
            plotData.queryRunFeature.cat1 = queryFeature;
            plotData.selectedSearchesTable.cat1 = searchFeaturesTable(end,:);

            options = generateOptionsStruct(groupByCell,isOverallExperimentPlot);
            options.plotColor.cat1 = color;
            options.plotColorOverlay.cat1 = colorOverlay;
            options.xLabelText.cat1 = [' '];
            options.yLabelText = ["Hotspot distance" + newline + "prctile [pixel]"];
            options.yLimits = [0,600];
            
            generateBarplotSearch(plotData,options,'prctile')
            xticklabels({'25%','50%','75%'})
               
        end
        
        % Title and saving

        plotMainTitle = ['Cell ', num2str(searchCell),' - ','Epoch ', num2str(searchEpoch),' - ',experimentName];
        sgtitle(plotMainTitle,'FontSize',14,'Interpreter', 'none')

        plotName = ['searchSummary_','Cell', num2str(searchCell), '_', 'Epoch', num2str(searchEpoch)];
        plotFullPath = fullfile(savePlotPath, plotName);

        if savePlots == 1
            exportgraphics(gcf, fullfile(strcat(plotFullPath,'.pdf')), 'ContentType', 'vector');
            saveas(gcf, plotFullPath, 'fig')
        end

         saveTableName = 'SearchFeaturesTable.mat';
         saveTablePath = fullfile(savePlotPath, saveTableName);
         save(saveTablePath, 'searchFeaturesTable');

    end
    
end
