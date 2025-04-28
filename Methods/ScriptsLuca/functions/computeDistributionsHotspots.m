depthBlue = 4;
depthRed = 3;
nBins = 40;
nSearches = 5;

initializeFig(0.9,0.9);
masterLayout = tiledlayout(3,nSearches);
masterLayout.TileSpacing = 'compact';
masterLayout.Padding = 'compact';

colorBlue = [85, 161, 254]./255;
colorRed = [255, 50, 58]./255;
colorGrey = [.92, .92, .92];
cmapBlue = createcolormap(colorGrey,colorBlue);
cmapRed = createcolormap(colorRed,colorGrey);
markers = {'-o', '-*', '-+', '-x', '-s', '-d', '-^', '-v', '->', '-<', '-p', '-h'};
legendText = {};

for iSearch = 1:nSearches
    
    if iSearch == 1
        responseType = 'excitatory';
        depth = depthBlue;
        cellX = cells1.Protocol{4,1}{1,1}{depth,1}.cellX;
        cellY = cells1.Protocol{4,1}{1,1}{depth,1}.cellY;
        responseMap = cells1.(8){4,1}.responseMap{1,1}(:,:,depth); 
        isResponseMap = cells1.(8){4,1}.isResponseMap{1,1}(:,:,depth);
        splitExpPath = split(cells1.(1){1},'\');
        expName = splitExpPath{end-1,1};
    elseif iSearch == 2
        responseType = 'excitatory';
        depth = depthBlue;
        cellX = cells2.Protocol{1,1}{1,1}{depth,1}.cellX;
        cellY = cells2.Protocol{1,1}{1,1}{depth,1}.cellY;
        responseMap = cells2.(8){1,1}.responseMap{1,1}(:,:,depth); 
        isResponseMap = cells2.(8){1,1}.isResponseMap{1,1}(:,:,depth);
        splitExpPath = split(cells2.(1){1},'\');
        expName = splitExpPath{end-1,1};
    elseif iSearch == 3
        responseType = 'excitatory';
        depth = depthBlue;
        cellX = cells3.Protocol{4,1}{1,1}{depth,1}.cellX;
        cellY = cells3.Protocol{4,1}{1,1}{depth,1}.cellY;
        responseMap = cells3.(8){4,1}.responseMap{1,1}(:,:,depth); 
        isResponseMap = cells3.(8){4,1}.isResponseMap{1,1}(:,:,depth);
        splitExpPath = split(cells3.(1){1},'\');
        expName = splitExpPath{end-1,1};
    elseif iSearch == 4
        responseType = 'inhibitory';
        depth = depthRed;
        cellX = cells4.Protocol{5,1}{3,1}{depth,1}.cellX;
        cellY = cells4.Protocol{5,1}{3,1}{depth,1}.cellY;
        responseMap = cells4.(8){5,1}.responseMap{3,1}(:,:,depth); 
        isResponseMap = cells4.(8){5,1}.isResponseMap{3,1}(:,:,depth);  
        splitExpPath = split(cells4.(1){1},'\');
        expName = splitExpPath{end-1,1};
    elseif iSearch == 5
        responseType = 'inhibitory';
        depth = depthRed;
        cellX = cells1.Protocol{3,1}{1,1}{depth,1}.cellX;
        cellY = cells1.Protocol{3,1}{1,1}{depth,1}.cellY;
        responseMap = cells1.(8){3,1}.responseMap{1,1}(:,:,depth); 
        isResponseMap = cells1.(8){3,1}.isResponseMap{1,1}(:,:,depth); 
        splitExpPath = split(cells1.(1){1},'\');
        expName = splitExpPath{end-1,1};
    end
    
    spotSizeX = 684/(2^depth);
    spotSizeY = 608/(2^depth);
    spotSize = spotSizeX * spotSizeY;
    legendText{end+1,1} = expName;
    
    if strcmp(responseType,'excitatory')
        color = colorBlue;
        cmap = cmapBlue;
        climits = [min(responseMap,[],'all'),0];
        responseMap(responseMap > 0) = 0;
    elseif strcmp(responseType,'inhibitory')
        color = colorRed;
        cmap = cmapRed;
        climits = [0,max(responseMap,[],'all')];
        responseMap(responseMap < 0) = 0;
    end
    
    [meanRadialResponseAUC,radialDistancesAUC] = radialavg(responseMap,nBins,cellX,cellY,responseType);
    [meanRadialResponseCount,radialDistancesCount] = radialcounts(isResponseMap,nBins,cellX,cellY);

    nexttile(masterLayout,iSearch,[1 1])
    imshow(isResponseMap)
    hold on
    scatter(cellX,cellY,50,'o','filled','MarkerFaceColor','g','MarkerEdgeColor','g');
    title(expName,'FontSize',14,'Interpreter','none')
    
    nexttile(masterLayout,iSearch + nSearches,[1 1])
    imshow(responseMap)
    colormap(gca, flip(cmap))
    caxis(climits); cbar = colorbar; ylabel(cbar,'Total charge [pC]','FontSize',8);
    hold on
    scatter(cellX,cellY,50,'o','filled','MarkerFaceColor','k','MarkerEdgeColor','k');
    
    nexttile(masterLayout,2*nSearches + 1,[1 floor(nSearches/2)])
    
    plot(radialDistancesAUC,meanRadialResponseAUC,markers{iSearch},'Color',color,'MarkerFaceColor',color);
    xlabel('Radial distance [pixels]','FontSize',14)
    ylabel('Normalized mean response','FontSize',14)
    legend(legendText,'FontSize',10,'Location','northeast','Interpreter','none')
    hold on
    
    nexttile(masterLayout,2*nSearches + ceil(nSearches/2) + 1,[1 floor(nSearches/2)])
    plot(radialDistancesCount,meanRadialResponseCount,markers{iSearch},'Color',color,'MarkerFaceColor',color);
    xlabel('Radial distance [pixels]','FontSize',14)
    ylabel('Normalized hotspot area [1]','FontSize',14)
    legend(legendText,'FontSize',10,'Location','northeast','Interpreter','none')
    hold on
    
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

function [meanRadialResponses, radialDistances] = radialcounts(isResponseMap,nBins,cellX,cellY)

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

end