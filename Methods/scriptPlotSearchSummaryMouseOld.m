% experimentName = 'DCN_ZI_Ephys_3_1';
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

columnNames = {'mouseName','mouseCell','mouseEpoch','runDepth','holdingVoltage','responseType','runProtocol','radii','hotspotCount','hotspotCountMaxValue','hotspotCountMaxRadius','hotspotCountMeanRadius','hotspotCountStdRadius','hotspotCountMedianRadius','hotspotCountModeRadius','hotspotCountPrctile','hotspotCountNorm','hotspotCountNormMaxValue','hotspotCountNormMaxRadius','hotspotCountNormMeanRadius','hotspotCountNormStdRadius','hotspotCountNormMedianRadius','hotspotCountNormModeRadius','hotspotCountNormPrctile','radialCharge','radialChargeMaxValue','radialChargeMaxRadius','radialChargeNorm','radialChargeNormMaxValue','radialChargeNormMaxRadius'};
dataTypes = {'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell'};
searchFeaturesTable = table('Size', [0, numel(columnNames)], 'VariableNames', columnNames, 'VariableTypes', dataTypes);
rowTable = 1; warning('off', 'MATLAB:table:RowsAddedExistingVars');

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
        else; responseType = 'inhibitory'; end
        
        nSearchDepth = size(responseMapsEpoch,3);
        
        initializeFig(0.9,0.9);
        nMaxDepthsPlot = 5;
        masterLayout = tiledlayout(6,2*nMaxDepthsPlot);
        masterLayout.TileSpacing = 'tight';
        masterLayout.Padding = 'tight';

        colorBlue = [85, 161, 254]./255;
        colorRed = [255, 50, 58]./255;
        colorGrey = [.92, .92, .92];
        cmapBlue = createcolormap(colorGrey,colorBlue);
        cmapRed = createcolormap(colorRed,colorGrey);
        
        for iSearchDepth = 1:nSearchDepth
            
            responseMapDepth = responseMapsEpoch(:,:,iSearchDepth);
            isResponseMapDepth = isResponseMapsEpoch(:,:,iSearchDepth);
            hotspotMapDepth = hotspotMapsEpoch(:,:,iSearchDepth);
            
            nBins = 40;
            spotSizeX = 684/(2^iSearchDepth);
            spotSizeY = 608/(2^iSearchDepth);
            spotSize = spotSizeX * spotSizeY;

            if strcmp(responseType,'excitatory')
                color = colorBlue;
                cmap = cmapBlue;
                climits = [min(responseMapDepth,[],'all'),0];
                responseMapDepth(responseMapDepth > 0) = 0;
            elseif strcmp(responseType,'inhibitory')
                color = colorRed;
                cmap = cmapRed;
                climits = [0,max(responseMapDepth,[],'all')];
                responseMapDepth(responseMapDepth < 0) = 0;
            end

            [radialCharge,radialChargeNorm,~] = radialavg(responseMapDepth,nBins,cellX,cellY,responseType);
            [hotspotCount,hotspotCountNorm,radii] = radialcounts(hotspotMapDepth,nBins,cellX,cellY);
          
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
            title(['Response map - Depth ' num2str(iSearchDepth)],'FontSize',12,'Interpreter','none')
            
            nexttile(masterLayout,2*(iSearchDepth-1) + 1 + 2*nMaxDepthsPlot,[1 2])
            imshow(hotspotMapDepth)
            hold on
            rectangle('Position', [0.5, 0.5, size(hotspotMapDepth,2), size(hotspotMapDepth,1)], 'EdgeColor', colorGrey, 'LineWidth', 1);
            hold on
            scatter(cellX,cellY,50,'o','filled','MarkerFaceColor','g','MarkerEdgeColor','g');
            title(['Hotspot map - Depth ' num2str(iSearchDepth)],'FontSize',12,'Interpreter','none')

            nexttile(masterLayout,2*(iSearchDepth-1) + 1 + 4*nMaxDepthsPlot,[1 2])
            imshow(responseMapDepth)
            colormap(gca, flip(cmap))
            caxis(climits); cbar = colorbar; ylabel(cbar,'Total charge [pC]','FontSize',8);
            hold on
            scatter(cellX,cellY,50,'o','filled','MarkerFaceColor','k','MarkerEdgeColor','k');
            title(['Charge map - Depth ' num2str(iSearchDepth)],'FontSize',12,'Interpreter','none')
            
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

function [meanRadialResponses, meanRadialResponsesNorm, radialDistances] = radialavg(responseMap,nBins,cellX,cellY,responseType)

dmdWidth = 684;
dmdHeight = 608;
dmdX = 1:dmdWidth;
dmdY = 1:dmdHeight;
[X,Y] = meshgrid(dmdX,dmdY);

X = X - cellX;
Y = Y - cellY;
spotDistance = sqrt(X.^2+Y.^2);
maxDistance = max(spotDistance,[],'all');

binWidth = maxDistance/(nBins);
radialDistances = (0:binWidth:maxDistance);
distanceBins = round(spotDistance/binWidth)+1;

distanceBinsVector = reshape(distanceBins,dmdWidth*dmdHeight,1);
responseMapVector = reshape(responseMap,dmdWidth*dmdHeight,1);

meanRadialResponses = accumarray(distanceBinsVector,responseMapVector,[],@mean);

    if strcmp(responseType,'excitatory')
        meanRadialResponsesNorm = meanRadialResponses'/abs(min(meanRadialResponses));
    elseif strcmp(responseType,'inhibitory')
        meanRadialResponsesNorm = meanRadialResponses'/abs(max(meanRadialResponses));
    end
    
    meanRadialResponses = meanRadialResponses';

end

function [meanRadialResponses, meanRadialResponsesNormalized, radialDistances] = radialcounts(isResponseMap,nBins,cellX,cellY)

dmdWidth = 684;
dmdHeight = 608;
dmdX = 1:dmdWidth;
dmdY = 1:dmdHeight;
[X,Y] = meshgrid(dmdX,dmdY);

X = X - cellX;
Y = Y - cellY;
spotDistance = sqrt(X.^2+Y.^2);
maxDistance = max(spotDistance,[],'all');

binWidth = maxDistance/(nBins);
radialDistances = (0:binWidth:maxDistance);
distanceBins = round(spotDistance/binWidth)+1;

distanceBinsVector = reshape(distanceBins,dmdWidth*dmdHeight,1);
responseMapVector = reshape(isResponseMap,dmdWidth*dmdHeight,1);

radialSurfaces = accumarray(distanceBinsVector, 1);
meanRadialResponses = accumarray(distanceBinsVector,responseMapVector);
meanRadialResponses = (meanRadialResponses');
meanRadialResponsesNormalized = (meanRadialResponses./radialSurfaces');

end