
% experimentName = 'DCN_ZI_Ephys_3_2';

% user = getenv('USER'); 
% matlabDirectory = ['/Users/' user '/Desktop/EPFL/Thesis/Matlab/'];
% mainDirectory = ['/Users/' user '/Desktop/EPFL/Thesis/Matlab/Experiments1/'];
% addpath(genpath(matlabDirectory));
% experimentDirectory = [mainDirectory experimentName];

% for Paolo's computer and Shun's main script
experimentDirectory = expPath;
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

for iSearchCell = 1:nSearchCells
    
    searchCell = searchData.Cell(iSearchCell);
    nSearchEpochs = size(searchData.Epochs{iSearchCell,1},1);
    
    for iSearchEpoch = 1:nSearchEpochs
        
        searchEpoch = str2double(regexp(searchData.Epochs{iSearchCell,1}{iSearchEpoch,1}, '(?<=epoch)\d+', 'match', 'once'));
        
        responseMapsEpoch = searchData.(8){iSearchCell,1}.responseMap{iSearchEpoch,1};
        isResponseMapsEpoch = searchData.(8){iSearchCell,1}.isResponseMap{iSearchEpoch,1};
        hotspotMapsEpoch = searchData.(8){iSearchCell,1}.hotspotMap{iSearchEpoch,1};
        
        cellX = searchData.Protocol{iSearchCell,1}{1,1}{1,1}.cellX;
        cellY = searchData.Protocol{iSearchCell,1}{1,1}{1,1}.cellY;
        holdingVoltage = searchData.Vhold{iSearchCell,1}(iSearchEpoch,1);
        
        if holdingVoltage <= -65; responseType = 'excitatory';
        elseif holdingVoltage >= 0; responseType = 'inhibitory';
        else; responseType = 'inhibitory'; end
        
        nSearchDepth = size(responseMapsEpoch,3);
        
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
        markers = {'-o', '-*', '-+', '-x', '-s', '-d', '-^', '-v', '->', '-<', '-p', '-h'};
        
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

            [meanRadialResponseAUC,radialDistancesAUC] = radialavg(responseMapDepth,nBins,cellX,cellY,responseType);
            [meanRadialResponseCount,meanRadialResponseCountNorm,radialDistancesCount] = radialcounts(isResponseMapDepth,nBins,cellX,cellY);

            [maxAUC,indMaxAUC] = max(meanRadialResponseAUC);
            [maxCount,indMaxCount] = max(meanRadialResponseCount);
            [maxCountNorm,indMaxCountNorm] = max(meanRadialResponseCountNorm);
            radiusMaxAUC = radialDistancesAUC(indMaxAUC);
            radiusMaxCount = radialDistancesCount(indMaxCount);
            radiusMaxCountNorm = radialDistancesCount(indMaxCountNorm);
            
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
            caxis([climits(1),climits(end)+0.001]); cbar = colorbar; ylabel(cbar,'Total charge [pC]','FontSize',8);
            hold on
            scatter(cellX,cellY,50,'o','filled','MarkerFaceColor','k','MarkerEdgeColor','k');
            title(['Charge map - Depth ' num2str(iSearchDepth)],'FontSize',12,'Interpreter','none')
            
            nexttile(masterLayout,2*(iSearchDepth-1) + 1 + 6*nMaxDepthsPlot,[1 2])
            plot(radialDistancesCount,meanRadialResponseCount,'-o','Color',color,'MarkerFaceColor',color);
            xlabel('Radial distance [pixels]','FontSize',10)
            ylabel('Hotspot area [pixel^2]','FontSize',10)
            hold on
            xline(radiusMaxCount,'-g','LineWidth',1)
            hold on
            
            nexttile(masterLayout,2*(iSearchDepth-1) + 1 + 8*nMaxDepthsPlot,[1 2])
            plot(radialDistancesCount,meanRadialResponseCountNorm,'-o','Color',color,'MarkerFaceColor',color);
            xlabel('Radial distance [pixels]','FontSize',10)
            ylabel(["Normalized" + newline + "hotspot area [1]"],'FontSize',10)
            hold on
            xline(radiusMaxCountNorm,'-g','LineWidth',1)
            hold on
            
            nexttile(masterLayout,2*(iSearchDepth-1) + 1 + 10*nMaxDepthsPlot,[1 2])
            plot(radialDistancesAUC,meanRadialResponseAUC,'o','Color',color,'MarkerFaceColor',color);
            xlabel('Radial distance [pixels]','FontSize',10)
            ylabel(["Normalized" + newline + "mean response"],'FontSize',10)
            hold on
            xline(radiusMaxAUC,'-g','LineWidth',1)
            hold on
   
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

    end
    
end

function [meanRadialResponses, radialDistances] = radialavg(responseMap,nBins,cellX,cellY,responseType)

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
        meanRadialResponses = meanRadialResponses'/abs(min(meanRadialResponses));
    elseif strcmp(responseType,'inhibitory')
        meanRadialResponses = meanRadialResponses'/abs(max(meanRadialResponses));
    end

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
meanRadialResponsesNormalized = (meanRadialResponses'./radialSurfaces);

end